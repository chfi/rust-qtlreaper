use ndarray::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::ops::Range;
use std::path::PathBuf;

// `Metadata` is really only used for parsing; it's the data above the
// header line in the genome data.
#[derive(Debug, PartialEq, Clone)]
struct Metadata {
    name: String,
    maternal: String,
    paternal: String,
    dataset_type: String, // "intercross" or "riset"
    heterozygous: String, // defaults to "H"
    unknown: String,      // defaults to "U"
}

impl Metadata {
    fn parse_genoprob(&self, geno: &str) -> f64 {
        if geno == self.maternal {
            -1.0
        } else if geno == self.paternal {
            1.0
        } else if geno == self.heterozygous {
            0.0
        } else if geno == self.unknown {
            99.0
        } else {
            panic!("Failed to parse genotype: {}\n{:?}", geno, self);
        }
    }

    fn parse_genotype(&self, geno: &str) -> Genotype {
        if geno == self.maternal {
            Genotype::Mat
        } else if geno == self.paternal {
            Genotype::Pat
        } else if geno == self.heterozygous {
            Genotype::Het
        } else if geno == self.unknown {
            Genotype::Unk
        } else {
            panic!("Failed to parse genotype: {}\n{:?}", geno, self);
        }
    }

    fn parse_dominance(&self, geno: &str) -> f64 {
        if geno == self.maternal || geno == self.paternal {
            0.0
        } else if geno == self.heterozygous || geno == self.unknown {
            1.0
        } else {
            panic!("Failed to parse genotype: {}\n{:?}", geno, self);
        }
    }

    fn parse_line(line: &str) -> Option<(&str, &str)> {
        let line = line.trim();
        if line.starts_with('#') {
            return None;
        }

        if line.starts_with('@') {
            let sep = line.find(':').unwrap();
            let name = &line[1..sep];
            let val = &line[sep + 1..];

            return Some((name, val));
        }

        None
    }

    // panic!s if the provided lines do not contain @name, @mat, and @pat fields
    fn from_lines(lines: Vec<&str>) -> Metadata {
        let mut name = None;
        // the type should be either `riset` or `intercross`; fix later
        let mut typ = None;

        let mut mat = None;
        let mut pat = None;
        let mut het = "H".into();
        let mut unk = "U".into();

        for line in lines.iter() {
            if let Some((n, v)) = Metadata::parse_line(line) {
                let val = v.into();
                match n {
                    "name" => name = Some(val),
                    "mat" => mat = Some(val),
                    "pat" => pat = Some(val),
                    "type" => typ = Some(val),
                    "het" => het = val,
                    "unk" => unk = val,
                    _ => (),
                }
            }
        }

        if name == None || mat == None || pat == None || typ == None {
            panic!(
                "Required metadata was not provided!\nname = {:?}\nmat = {:?}\npat = {:?}\ntype = {:?}",
                name, mat, pat, typ
            );
        }

        Metadata {
            name: name.unwrap(),
            maternal: mat.unwrap(),
            paternal: pat.unwrap(),
            dataset_type: typ.unwrap(),
            heterozygous: het,
            unknown: unk,
        }
    }
}

#[derive(Debug, PartialEq, Clone, Serialize)]
pub struct Marker {
    pub name: String,
    pub centi_morgan: f64,
    pub mega_basepair: Option<f64>,
    pub chromosome: String,
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
pub enum Genotype {
    Mat,
    Pat,
    Het,
    Unk,
}

#[derive(Debug, PartialEq, Clone)]
pub struct Locus {
    dominance: Option<Array1<f64>>,
    pub genotype: Array1<Genotype>,
    genoprob: Array1<f64>,
    pub marker: Marker,
}

/// UnknownIntervals holds a list of ranges of unknown genotypes, per strain
struct UnknownIntervals(Vec<Vec<Range<usize>>>);

impl Locus {
    // corresponds to lines 950-1044 in dataset.c
    fn parse_line(
        metadata: &Metadata,
        has_mb: bool,
        // header: &DatasetHeader,
        dominance: bool,
        line: &str,
    ) -> (String, Locus) {
        // Example locus is: "1	D1Mit1	8.3	B6	B6	D	D"
        // where the first three columns are chromosome, name, cM;
        // remaining columns are the genotypes

        let words: Vec<_> = line.split_terminator('\t').collect();

        let chromosome = words[0].to_string();
        let name = words[1].into();
        let centi_morgan = words[2].parse::<f64>().unwrap();
        let mega_basepair = if has_mb {
            words[3].parse::<f64>().ok()
        } else {
            None
        };

        let marker = Marker {
            name,
            centi_morgan,
            mega_basepair,
            chromosome: chromosome.clone(),
        };

        let range = if has_mb { 4.. } else { 3.. };

        let genotype = words[range.clone()]
            .iter()
            .map(|g| metadata.parse_genotype(g))
            .collect();

        let genoprob = words[range.clone()]
            .iter()
            .map(|g| metadata.parse_genoprob(g))
            .collect();

        let dominance = if dominance {
            Some(
                words[range]
                    .iter()
                    .map(|g| metadata.parse_dominance(g))
                    .collect(),
            )
        } else {
            None
        };

        (
            chromosome,
            Locus {
                genotype,
                genoprob,
                dominance,
                marker,
            },
        )
    }

    /// Steps through a list of genotypes per strain, building up a list of ranges of missing data for each strain
    fn step_many_unknown_intervals_mut(
        state: &mut (Vec<Option<usize>>, Vec<Vec<Range<usize>>>),
        next: (usize, &ArrayView1<Genotype>),
    ) {
        let (ix, genotype) = next;

        for (strain_ix, geno) in genotype.iter().enumerate() {
            if let Genotype::Unk = geno {
                match state.0[strain_ix] {
                    None => state.0[strain_ix] = Some(ix),
                    Some(start) => state.0[strain_ix] = Some(start),
                }
            } else if let Some(start) = state.0[strain_ix] {
                state.1[strain_ix].push(start..ix);
                state.0[strain_ix] = None;
            }
        }
    }

    fn find_unknown_intervals(loci: &[Locus]) -> UnknownIntervals {
        let n_strains = loci.first().unwrap().genotype.len();
        let mut state = (vec![None; n_strains], vec![Vec::new(); n_strains]);

        for (locus_ix, locus) in loci.iter().enumerate() {
            Self::step_many_unknown_intervals_mut(
                &mut state,
                (locus_ix, &locus.genotype.slice(s![..])),
            )
        }

        UnknownIntervals(state.1)
    }

    fn estimate_unknown_locus(
        strain_ix: usize,
        locus: &mut Locus,
        prev: &Locus,
        next: &Locus,
    ) {
        let rec_1 = (locus.cm() - prev.cm()) / 100.0;
        let rec_2 = (next.cm() - locus.cm()) / 100.0;
        let rec_0 = (next.cm() - prev.cm()) / 100.0;

        let f1 = (1.0 - f64::exp(-2.0 * rec_1)) / 2.0;
        let f2 = (1.0 - f64::exp(-2.0 * rec_2)) / 2.0;
        let f0 = (1.0 - f64::exp(-2.0 * rec_0)) / 2.0;

        let r_0 = (1.0 - f1) * (1.0 - f2) / (1.0 - f0);
        let r_1 = f1 * (1.0 - f2) / f0;
        let r_2 = f2 * (1.0 - f1) / f0;
        let r_3 = f1 * f2 / (1.0 - f0);

        let prev_geno = prev.genotype[strain_ix];
        let next_geno = next.genotype[strain_ix];

        let new_genoprob = match (prev_geno, next_geno) {
            (Genotype::Mat, Genotype::Mat) => 1.0 - 2.0 * r_0,
            (Genotype::Het, Genotype::Mat) => 1.0 - r_0 - r_1,
            (Genotype::Pat, Genotype::Mat) => 1.0 - 2.0 * r_1,
            (Genotype::Mat, Genotype::Het) => r_1 - r_0,
            (Genotype::Het, Genotype::Het) => 0.0,
            (Genotype::Pat, Genotype::Het) => r_0 - r_1,
            (Genotype::Mat, Genotype::Pat) => 2.0 * r_1 - 1.0,
            (Genotype::Het, Genotype::Pat) => r_0 + r_1 - 1.0,
            (Genotype::Pat, Genotype::Pat) => 2.0 * r_0 - 1.0,
            _ => panic!("Genotype was unknown when it shouldn't be!"),
        };

        locus.genoprob[strain_ix] = new_genoprob;

        if let Some(d) = &mut locus.dominance {
            // if dominance {
            let new_dominance = match (prev_geno, next_geno) {
                (Genotype::Mat, Genotype::Mat) => 2.0 * r_0 * r_3,
                (Genotype::Het, Genotype::Mat) => r_1 * (r_2 + r_3),
                (Genotype::Pat, Genotype::Mat) => 2.0 * r_1 * r_2,
                (Genotype::Mat, Genotype::Het) => r_1 * r_0 + r_2 * r_3,
                (Genotype::Het, Genotype::Het) => {
                    let w = ((1.0 - f0) * (1.0 - f0))
                        / (1.0 - 2.0 * f0 * (1.0 - f0));
                    1.0 - 2.0 * w * r_0 * r_3 - 2.0 * (1.0 - w) * r_1 * r_2
                }
                (Genotype::Pat, Genotype::Het) => r_0 * r_1 + r_2 * r_3,
                (Genotype::Mat, Genotype::Pat) => 2.0 * r_1 * r_2,
                (Genotype::Het, Genotype::Pat) => r_1 * (r_2 + r_3),
                (Genotype::Pat, Genotype::Pat) => 2.0 * r_1 * r_3,
                _ => panic!("Genotype was unknown when it shouldn't be!"),
            };

            d[strain_ix] = new_dominance;
        }
    }

    fn estimate_unknown_genotypes(
        loci: &mut [Locus],
        intervals: UnknownIntervals,
    ) {
        for (strain_ix, strain) in intervals.0.iter().enumerate() {
            for range in strain {
                for locus_ix in range.clone() {
                    // wasteful to be cloning here, really
                    let prev = &loci[range.start - 1].clone();
                    let next = &loci[range.end].clone();
                    let locus = loci.get_mut(locus_ix).unwrap();
                    Locus::estimate_unknown_locus(strain_ix, locus, prev, next);
                }
            }
        }
    }

    pub fn cm(&self) -> f64 {
        self.marker.centi_morgan
    }

    // allocating this every step is probably slowing things down (it was twice as fast without)
    pub fn genotypes_subset(&self, strain_ixs: &[usize]) -> Vec<f64> {
        strain_ixs.iter().map(|ix| self.genoprob[*ix]).collect()
    }

    pub fn genotypes_subindices(
        &self,
        indices: &[usize],
        subset: &mut Vec<f64>,
    ) {
        for (data_ix, ix) in indices.iter().enumerate() {
            subset[data_ix] = self.genoprob[*ix];
        }
    }

    pub fn dominance_subset(&self, strain_ixs: &[usize]) -> Vec<f64> {
        match &self.dominance {
            None => {
                panic!("Attempted to extract dominance subset, but dataset has no dominance data")
            }
            Some(d) => strain_ixs.iter().map(|ix| d[*ix]).collect(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Genome {
    chr_order: Vec<String>,
    pub chromosomes: BTreeMap<String, Vec<Locus>>, // chromosomes: Vec<(String, Vec<Locus>)>
}

/// Iterator that steps through the genome in order
pub struct GenomeIter<'a> {
    keys: Vec<String>,
    chromosomes: &'a BTreeMap<String, Vec<Locus>>,
}

impl<'a> Iterator for GenomeIter<'a> {
    type Item = &'a Vec<Locus>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.keys.is_empty() {
            None
        } else {
            let chr = self.keys.remove(0);
            self.chromosomes.get(&chr)
        }
    }
}

impl Genome {
    fn new() -> Genome {
        Genome {
            chr_order: Vec::new(),
            chromosomes: BTreeMap::new(),
        }
    }

    fn or_push_chromosome(&mut self, chr: String) -> &mut Vec<Locus> {
        if self.chr_order.iter().find(|&c| c == &chr).is_none() {
            self.chr_order.push(chr.clone());
        }

        self.chromosomes.entry(chr).or_insert_with(Vec::new)
    }

    fn push_locus(&mut self, chr: String, locus: Locus) {
        self.or_push_chromosome(chr).push(locus);
    }

    /// Iterates through the chromosomes in the order they were added to the genotype
    pub fn iter(&self) -> GenomeIter<'_> {
        GenomeIter {
            keys: self.chr_order.clone(),
            chromosomes: &self.chromosomes,
        }
    }

    /// Mutably iterates through the chromosomes, using the sorted order from BTreeMap
    fn iter_mut(
        &mut self,
    ) -> impl Iterator<Item = (&'_ String, &'_ mut Vec<Locus>)> {
        self.chromosomes.iter_mut()
    }

    pub fn find_locus(&self, name: &str) -> Option<&Locus> {
        let mut result = None;
        for loci in self.iter() {
            if let Some(l) = loci.iter().find(|locus| locus.marker.name == name)
            {
                result = Some(l)
            }
        }

        result
    }

    fn chromosome_interval(chromosome: &[Locus], interval: f64) -> Vec<Locus> {
        let find_adj_known = |geno_ix: usize, ix: usize| {
            let (lhs, rhs) = chromosome.split_at(ix + 1);
            let pred =
                |locus: &&Locus| locus.genotype[geno_ix] != Genotype::Unk;
            let err = "Error when searching for adjacent loci";
            (
                lhs.iter().rev().find(pred).expect(err),
                rhs.iter().find(pred).expect(err),
            )
        };

        let interval = interval.min(1.0);
        let mut interval_chromosome = Vec::new();

        chromosome.iter().enumerate().for_each(|(ix, locus)| {
            let mut cur_cm = locus.cm();
            let next_locus = &chromosome[(ix + 1).min(chromosome.len() - 1)];
            let mut first = true;

            loop {
                if first {
                    interval_chromosome.push(locus.clone());
                } else {
                    let mut new_locus = locus.clone();
                    new_locus.marker.name = String::from(" - ");
                    new_locus.marker.centi_morgan = cur_cm;
                    for (geno_ix, _geno) in locus.genotype.iter().enumerate() {
                        let (prev, next) = find_adj_known(geno_ix, ix);
                        Locus::estimate_unknown_locus(
                            geno_ix,
                            &mut new_locus,
                            prev,
                            next,
                        );
                    }
                    interval_chromosome.push(new_locus);
                }

                cur_cm += interval;
                first = false;

                if cur_cm >= next_locus.cm() {
                    break;
                }
            }
        });

        interval_chromosome
    }

    fn interval_mapped(&self, interval: f64) -> Genome {
        let mut chromosomes = BTreeMap::new();
        for (chr, loci) in self.chromosomes.iter() {
            let new_loci = Self::chromosome_interval(&loci, interval);
            chromosomes.insert(chr.clone(), new_loci);
        }

        Genome {
            chr_order: self.chr_order.clone(),
            chromosomes,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Dataset {
    metadata: Metadata,
    pub genome: Genome,
    strains: Vec<String>,
    pub dominance: bool, // true if dataset type is "intercross"
    has_mb: bool,
}

impl Dataset {
    fn new(metadata: Metadata, strains: Vec<String>, has_mb: bool) -> Dataset {
        let dominance = metadata.dataset_type == "intercross";
        Dataset {
            metadata,
            strains,
            genome: Genome::new(),
            dominance,
            has_mb,
        }
    }

    pub fn has_mb(&self) -> bool {
        self.has_mb
    }

    /// Corresponds to `addintervals` in C implementation
    pub fn interval_mapped_clone(&self, interval: f64) -> Dataset {
        let genome = self.genome.interval_mapped(interval);
        Dataset {
            genome,
            metadata: self.metadata.clone(),
            strains: self.strains.clone(),
            dominance: self.dominance,
            has_mb: self.has_mb,
        }
    }

    pub fn strains(&self) -> &[String] {
        &self.strains
    }

    pub fn strain_indices(&self, strains: &[String]) -> Vec<usize> {
        strains
            .iter()
            .map(|s| self.strains.iter().position(|p| p == s).unwrap())
            .collect()
    }

    pub fn n_loci(&self) -> usize {
        self.genome
            .chromosomes
            .iter()
            .map(|(_, loci)| loci.len())
            .sum()
    }

    fn parse_dataset_header(line: &str) -> (bool, Vec<String>) {
        let header_words: Vec<_> = line.split_terminator('\t').collect();

        let has_mb = match header_words.get(3) {
            None => panic!(
                "Dataset header had less than four elements; no strains!"
            ),
            Some(w) => *w == "Mb",
        };

        let skip_n = if has_mb { 4 } else { 3 };

        let strains = header_words
            .into_iter()
            .skip(skip_n)
            .map(String::from)
            .collect();

        (has_mb, strains)
    }

    pub fn read_file(path: &PathBuf) -> Dataset {
        let f = File::open(path)
            .unwrap_or_else(|_| panic!("Error opening file {:?}", path));

        let reader = BufReader::new(f);
        let mut lines = reader.lines();

        let has_mb;
        let strains;

        let mut metadata_lines = vec![];

        loop {
            match lines.next() {
                None => {
                    panic!("Reached end of file before parsing dataset header")
                }
                Some(l) => {
                    let ll = l.unwrap();
                    if ll.starts_with("Chr	Locus	cM") {
                        let header = Dataset::parse_dataset_header(&ll);
                        has_mb = header.0;
                        strains = header.1;
                        break;
                    } else {
                        metadata_lines.push(ll);
                    }
                }
            }
        }

        let metadata = Metadata::from_lines(
            metadata_lines.iter().map(String::as_str).collect(),
        );

        let mut dataset = Dataset::new(metadata, strains, has_mb);

        for line in lines {
            let (chr, locus) = Locus::parse_line(
                &dataset.metadata,
                has_mb,
                dataset.dominance,
                &line.unwrap(),
            );
            dataset.genome.push_locus(chr, locus);
        }
        dataset.estimate_unknown();

        dataset
    }

    // Corresponds to lines 1071-1152 in dataset.c
    fn estimate_unknown(&mut self) {
        // first replace any cases of "Unknown" in the first and last loci of each chromosome

        for (_chr, loci) in self.genome.iter_mut() {
            let replace_genotype = |locus: Option<&mut Locus>| {
                let l = locus.unwrap();

                let gt = &mut l.genotype;
                let gp = &mut l.genoprob;

                azip!(mut gt (gt), mut gp (gp) in {
                    if let Genotype::Unk = *gt {
                        *gt = Genotype::Het;
                        *gp = 0.0;
                    }
                });
            };

            replace_genotype(loci.first_mut());
            replace_genotype(loci.last_mut());
        }

        for (_chr, loci) in self.genome.iter_mut() {
            // then, for each chromosome, construct the intervals of
            // unknown genotypes
            let unk = Locus::find_unknown_intervals(loci);

            // ... and use those intervals to estimate the
            // missing genotypes
            Locus::estimate_unknown_genotypes(loci, unk);
        }
    }
}

#[derive(Debug, Serialize)]
pub struct QTL {
    pub lrs: f64,
    pub additive: f64,
    pub dominance: Option<f64>,
    pub marker: Marker,
}

impl QTL {
    pub fn new(
        marker: Marker,
        lrs: f64,
        additive: f64,
        dominance: Option<f64>,
    ) -> QTL {
        QTL {
            lrs,
            additive,
            dominance,
            marker,
        }
    }
}

// formatter for QTL, used to print the regression output, tab-delimited
impl fmt::Display for QTL {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{:.*}",
            self.marker.name,
            self.marker.chromosome,
            3,
            self.marker.centi_morgan
        )?;

        if let Some(mb) = self.marker.mega_basepair {
            write!(f, "\t{:.*}", 3, mb)?;
        }

        write!(f, "\t{:.*}\t{:.*}", 3, self.lrs, 3, self.additive)?;

        if let Some(d) = self.dominance {
            write!(f, "\t{:.*}", 3, d)?;
        }

        Ok(())
    }
}

pub struct Traits {
    pub strains: Vec<String>,
    pub traits: Vec<(String, Vec<f64>)>,
}

impl Traits {
    pub fn read_file(path: &PathBuf) -> Traits {
        let f = File::open(path)
            .unwrap_or_else(|_| panic!("Error opening traits file {:?}", path));

        let reader = BufReader::new(f);
        let mut lines = reader.lines();

        let strains = match lines.next() {
            None => panic!("Reached end of file before parsing traits header"),
            Some(l) => {
                let ll = l.unwrap();
                if ll.starts_with("Trait") {
                    ll.split_terminator('\t')
                        .skip(1)
                        .map(ToString::to_string)
                        .collect()
                } else {
                    panic!(
                        "Traits file did not begin with \"Trait\", aborting"
                    );
                }
            }
        };

        let mut traits = Vec::new();

        for line in lines {
            let ll = line.unwrap();
            let mut words = ll.split_terminator('\t');
            let key = words.next().unwrap().to_string();
            let values = words.map(|s| s.parse::<f64>().unwrap()).collect();
            traits.push((key, values));
        }

        Traits { strains, traits }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_can_parse_header() {
        let header1 = "Chr	Locus	cM	BXD1	BXD2	BXD5	BXD6";
        let (has_mb_1, strains_1) = Dataset::parse_dataset_header(header1);

        assert_eq!(false, has_mb_1);
        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], strains_1);

        let header2 = "Chr	Locus	cM	Mb	BXD1	BXD2	BXD5	BXD6";
        let (has_mb_2, strains_2) = Dataset::parse_dataset_header(header2);

        assert_eq!(true, has_mb_2);
        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], strains_2);
    }

    #[test]
    fn it_can_parse_metadata() {
        let lines = vec![
        "#comment line always start with a '#'",
        "#type riset or intercross",
        "@type:riset",
        "@name:BXD",
        "#abbreviation of maternal or paternal parents",
        "@mat:B6",
        "@pat:D",
        "#heterozygous , optional, default is \"H\"",
        "@het:H",
        "#Unknown , optional, default is \"U\"",
        "@unk:U",
        "Chr	Locus	cM	BXD1	BXD2	BXD5	BXD6	BXD8	BXD9	BXD11	BXD12	BXD13	BXD14	BXD15	BXD16	BXD18	BXD19	BXD20	BXD21	BXD22	BXD23	BXD24	BXD25	BXD27	BXD28	BXD29	BXD30	BXD31	BXD32	BXD33	BXD34	BXD35	BXD36	BXD37	BXD38	BXD39	BXD40	BXD42",
            ];

        assert_eq!(
            Metadata::from_lines(lines),
            Metadata {
                name: "BXD".into(),
                maternal: "B6".into(),
                paternal: "D".into(),
                dataset_type: "riset".into(),
                heterozygous: "H".into(),
                unknown: "U".into(),
            }
        );
    }

    #[test]
    fn it_can_estimate_unknown_genotypes() {
        let genotypes = vec![
            array![Genotype::Mat, Genotype::Mat, Genotype::Pat],
            array![Genotype::Unk, Genotype::Pat, Genotype::Unk],
            array![Genotype::Unk, Genotype::Unk, Genotype::Pat],
            array![Genotype::Unk, Genotype::Unk, Genotype::Unk],
            array![Genotype::Pat, Genotype::Mat, Genotype::Unk],
            array![Genotype::Pat, Genotype::Mat, Genotype::Mat],
        ];

        let genoprobs = vec![
            array![-1.0, -1.0, 1.0],
            array![99.0, 1.0, 99.0],
            array![99.0, 99.0, 1.0],
            array![99.0, 99.0, 99.0],
            array![1.0, -1.0, 99.0],
            array![1.0, -1.0, -1.0],
        ];

        let mk_locus = |name, cm, genotype, genoprob| Locus {
            marker: Marker {
                name: String::from(name),
                centi_morgan: cm,
                mega_basepair: None,
                chromosome: String::from("1"),
            },
            dominance: None,
            genotype,
            genoprob,
        };

        let loci_new = vec![
            mk_locus(
                "Mk1",
                10.0,
                array![Genotype::Mat, Genotype::Mat, Genotype::Pat],
                array![-1.0, -1.0, 1.0],
            ),
            mk_locus(
                "Mk2",
                30.3,
                array![Genotype::Unk, Genotype::Pat, Genotype::Unk],
                array![-0.18523506128077272, 1.0, 0.9616255554798838],
            ),
            mk_locus(
                "Mk3",
                40.1,
                array![Genotype::Unk, Genotype::Unk, Genotype::Pat],
                array![0.18906668494617707, 0.3421367343627405, 1.0],
            ),
            mk_locus(
                "Mk4",
                50.2,
                array![Genotype::Unk, Genotype::Unk, Genotype::Unk],
                array![
                    0.5826065314914579,
                    -0.3223330030526561,
                    0.3223330030526552,
                ],
            ),
            mk_locus(
                "Mk5",
                60.3,
                array![Genotype::Pat, Genotype::Mat, Genotype::Unk],
                array![1.0, -1.0, -0.34213673436274084],
            ),
            mk_locus(
                "Mk6",
                70.1,
                array![Genotype::Pat, Genotype::Mat, Genotype::Mat],
                array![1.0, -1.0, -1.0],
            ),
        ];

        let mut loci = vec![
            mk_locus("Mk1", 10.0, genotypes[0].clone(), genoprobs[0].clone()),
            mk_locus("Mk2", 30.3, genotypes[1].clone(), genoprobs[1].clone()),
            mk_locus("Mk3", 40.1, genotypes[2].clone(), genoprobs[2].clone()),
            mk_locus("Mk4", 50.2, genotypes[3].clone(), genoprobs[3].clone()),
            mk_locus("Mk5", 60.3, genotypes[4].clone(), genoprobs[4].clone()),
            mk_locus("Mk6", 70.1, genotypes[5].clone(), genoprobs[5].clone()),
        ];

        let unk = Locus::find_unknown_intervals(&loci);

        Locus::estimate_unknown_genotypes(&mut loci, unk);

        assert_eq!(loci, loci_new);
    }

}
