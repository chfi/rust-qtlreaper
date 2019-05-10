use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;

#[derive(Debug)]
pub struct Locus {
    pub name: String,
    dominance: Option<Vec<f64>>,
    pub genotype: Vec<(Genotype, f64)>,
    centi_morgan: f64,
    mega_basepair: Option<f64>,
}

/// UnknownIntervals holds a list of ranges of unknown genotypes, per strain
struct UnknownIntervals(Vec<Vec<Range<usize>>>);

impl UnknownIntervals {
    fn empty(n: usize) -> UnknownIntervals {
        UnknownIntervals(Vec::with_capacity(n))
    }
}

/// Steps through a list of genotypes (corresponding to contiguous
/// loci for one strain), building up a list of ranges of missing data
fn step_unknown_intervals_mut(
    state: &mut (Option<usize>, Vec<Range<usize>>),
    next: (usize, Genotype),
) {
    let (ix, geno) = next;

    if let Genotype::Unk = geno {
        match state.0 {
            None => state.0 = Some(ix),
            Some(start) => state.0 = Some(start),
        }
    } else {
        if let Some(start) = state.0 {
            state.1.push(start..ix);
            state.0 = None;
        }
    }
}

/// Steps through a list of genotypes per strain, building up a list of ranges of missing data for each strain
fn step_many_unknown_intervals_mut(
    state: &mut (Vec<Option<usize>>, Vec<Vec<Range<usize>>>),
    next: (usize, &[Genotype]),
) {
    let (ix, genotype) = next;

    for (strain_ix, geno) in genotype.iter().enumerate() {
        if let Genotype::Unk = geno {
            match state.0[strain_ix] {
                None => state.0[strain_ix] = Some(ix),
                Some(start) => state.0[strain_ix] = Some(start),
            }
        } else {
            if let Some(start) = state.0[strain_ix] {
                state.1[strain_ix].push(start..ix);
                state.0[strain_ix] = None;
            }
        }
    }
}

impl Locus {
    // corresponds to lines 950-1044 in dataset.c
    fn parse_line(metadata: &Metadata, line: &str) -> (String, Locus) {
        // Example locus is: "1	D1Mit1	8.3	B6	B6	D	D"
        // where the first three columns are chromosome, name, cM;
        // remaining columns are the genotypes

        let words: Vec<_> = line.split_terminator('\t').collect();

        // let chromosome = String::from(words.next().unwrap());
        let chromosome = String::from(words[0]);
        let name = String::from(words[1]);
        let centi_morgan = words[2].parse::<f64>().unwrap();

        let genotype = words[3..]
            .iter()
            .map(|g| metadata.parse_genotype(g))
            .collect();

        (
            chromosome,
            Locus {
                name,
                centi_morgan,
                genotype,
                dominance: None,
                mega_basepair: None,
            },
        )
    }

    pub fn cm(&self) -> f64 {
        self.centi_morgan
    }

    fn unknown_data_intervals(strains: &[String], loci: &[Locus]) -> UnknownIntervals {
        // let n_strains = strains.len();
        // the current state tracks the current left-hand index for each strain (if the start of an unknown interval has been found), and the unknown intervals so far.
        let mut state: (Vec<Option<usize>>, UnknownIntervals) =
            (Vec::new(), UnknownIntervals::empty(strains.len()));

        state.1
    }
}

#[derive(Debug)]
pub struct Chromosome {
    pub name: String,
    pub loci: Vec<Locus>,
    pub size: usize,
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
pub enum Genotype {
    Mat,
    Pat,
    Het,
    Unk,
}

#[derive(Debug, PartialEq)]
pub struct Metadata {
    pub name: String,
    pub maternal: String,
    pub paternal: String,
    pub dataset_type: String,
    pub heterozygous: String, // defaults to "H"
    pub unknown: String,      // defaults to "U"
}

impl Metadata {
    fn new(name: &str, maternal: &str, paternal: &str, dataset_type: &str) -> Metadata {
        Metadata {
            name: String::from(name),
            maternal: String::from(maternal),
            paternal: String::from(paternal),
            dataset_type: String::from(dataset_type),
            heterozygous: String::from("H"),
            unknown: String::from("U"),
        }
    }

    fn parse_genotype(&self, geno: &str) -> (Genotype, f64) {
        if geno == self.maternal.as_str() {
            (Genotype::Mat, -1.0)
        } else if geno == self.paternal.as_str() {
            (Genotype::Pat, 1.0)
        } else if geno == self.heterozygous.as_str() {
            (Genotype::Het, 0.0)
        } else if geno == self.unknown.as_str() {
            (Genotype::Unk, 99.0)
        } else {
            panic!("Failed to parse genotype: {}\n{:?}", geno, self);
        }
    }

    fn parse_line(line: &str) -> Option<(&str, &str)> {
        if line.starts_with("#") {
            return None;
        }

        if line.starts_with("@") {
            let sep = line.find(':').unwrap();
            let name = &line[1..sep];
            let val = &line[sep + 1..];

            return Some((name, val));
        }

        None
    }

    // panic!s if the provided lines do not contain @name, @mat, and @pat fields
    pub fn from_lines(lines: Vec<&str>) -> Metadata {
        let mut name: Option<String> = None;
        let mut mat: Option<String> = None;
        let mut pat: Option<String> = None;

        // the type should be either `riset` or `intercross`; fix later
        let mut typ: Option<String> = None;
        let mut het = String::from("H");
        let mut unk = String::from("U");

        for line in lines.iter() {
            if let Some((n, v)) = Metadata::parse_line(line) {
                match n {
                    "name" => name = Some(String::from(v)),
                    "mat" => mat = Some(String::from(v)),
                    "pat" => pat = Some(String::from(v)),
                    "type" => typ = Some(String::from(v)),
                    "het" => het = String::from(v),
                    "unk" => unk = String::from(v),
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

#[derive(Debug)]
pub struct Dataset {
    metadata: Metadata,
    header: DatasetHeader,
    chromosomes: HashMap<String, Vec<Locus>>,
    strains: Vec<String>,
}

impl Dataset {
    // Corresponds to the outer loop in lines 941-1045,
    // and the loop at lines 1050-1069, in dataset.c
    // fn parse_dataset(metadata: &Metadata, header: &DatasetHeader) -> Dataset {

    // partition the dataset
    // }

    pub fn new(metadata: Metadata, header: DatasetHeader, strains: Vec<String>) -> Dataset {
        Dataset {
            metadata,
            header,
            strains,
            chromosomes: HashMap::new(),
        }
    }

    pub fn chromosomes(&self) -> &HashMap<String, Vec<Locus>> {
        &self.chromosomes
    }

    fn strain_ix(&self, strain: &str) -> Option<usize> {
        self.strains.iter().position(|s| s == &strain)
    }

    pub fn show_self(&self) {
        let ks = self.chromosomes.keys();
        for k in ks {
            println!("{}", k);
        }
    }

    // fn strain_genotype(&mut self, strain: &str) -> Vec<(String, String, (Genotype, f64))> {
    fn strain_genotype(&mut self, strain: &str) -> HashMap<String, Vec<(Genotype, f64)>> {
        let strain_ix = match self.strains.iter().position(|s| s == &strain) {
            None => panic!("Tried to find genotype for nonexistent strain {}", strain),
            Some(s) => s,
        };

        let mut strain_chromosomes = HashMap::new();

        for (chr, loci) in self.chromosomes.iter() {
            // let (chr, locus) = Locus::parse_line(&self.metadata, line);
            let strain_loci = loci.iter().map(|locus| locus.genotype[strain_ix]).collect();

            strain_chromosomes.insert(chr.to_string(), strain_loci);
        }

        strain_chromosomes
    }

    /*
    // returns the previous and next known loci's genotype and position in centimorgan
    fn adj_known_loci(
        &self,
        chromosome: &str,
        locus_ix: usize,
        strain_ix: usize,
        // strain: &str,
    ) -> ((Genotype, f64), (Genotype, f64)) {
        // let (strain_ix, _) = self
        //     .strains
        //     .iter()
        //     .enumerate()
        //     .find(|(_, s)| s == &strain)
        //     .unwrap();

        let loci = self.chromosomes.get(chromosome).unwrap();

        let find_adj_known = |step: i32| {
            let mut ix = (locus_ix as i32) + step;
            loop {
                let (geno, _val) = &loci[ix as usize].genotype[strain_ix];
                if *geno != Genotype::Unk {
                    return ix as usize;
                } else {
                    ix += step;
                }
            }
        };

        let find_adj_known_2 = |step: i32| {
            let mut ix = (locus_ix as i32) + step;
            loop {
                let (geno, _val) = &loci[ix as usize].genotype[strain_ix];
                if *geno != Genotype::Unk {
                    let loc = &loci[ix as usize];
                    let (adj_geno, _) = loc.genotype[strain_ix];
                    let cm = loc.cm();
                    return (adj_geno, cm);
                } else {
                    ix += step;
                }
            }
        };

        let prev = &loci[find_adj_known(-1)];
        let (prev_geno, _) = prev.genotype[strain_ix];

        let next = &loci[find_adj_known(1)];
        let (next_geno, _) = next.genotype[strain_ix];

        ((prev_geno, prev.cm()), (next_geno, next.cm()))
    }
    */

    // parse a line with a locus, adding it to the corresponding
    // chromosome in the dataset. NOTE: we assume the input data is
    // sorted by chromosome and marker position!
    pub fn parse_locus(&mut self, line: &str) {
        let (chr, locus) = Locus::parse_line(&self.metadata, line);

        let loci = self.chromosomes.entry(chr).or_insert_with(|| Vec::new());

        loci.push(locus);

        // match self.chromosomes.get_mut(&chr) {
        //     None => self.chromosomes.insert(chr, vec![locus]),
        //     Some(loci) => loci.push(locus),
        // }
        // match self.chromosomes.get_mut(&chr) {
        //     None =>
        //     Some(loci) => loci.push(locus);,
        // }
        // self.chromosomes
        //     .entry(chr)
        //     .
        //     .or_insert_with(|loci| loci.push(locus));
    }

    // pub fn parse_loci(&mut self, lines: Vec<&str>) {

    // }

    // Corresponds to lines 1071-1152 in dataset.c
    pub fn estimate_unknown(&mut self) {
        // first replace any cases of "Unknown" in the first and last loci of each chromosome

        let replace = |geno: &mut Genotype, val: &mut f64| {
            if let Genotype::Unk = *geno {
                *geno = Genotype::Het;
                *val = 0.0;
            }
        };

        for (_chr, loci) in self.chromosomes.iter_mut() {
            let replace_genotype = |locus: Option<&mut Locus>| {
                locus
                    .unwrap()
                    .genotype
                    .iter_mut()
                    .for_each(|(geno, val)| replace(geno, val))
            };

            replace_genotype(loci.first_mut());
            replace_genotype(loci.last_mut());
        }

        // then, for each strain,
        for (strain_ix, strain) in self.strains.iter().enumerate() {
            // walk through the loci, e

            for (chr, loci) in self.chromosomes.iter() {}
        }

        ///////// OLD

        // then walk through each locus, estimating the value for each remaining "Unknown" case

        /*
                    m = k-1;
                    while ((m>=0) && (((Locus *)(cptr->loci[m]))->txtstr[j] == GENOSYMBOL[3]))
                      m--;
                    n = k+1;
                    while ( (n<=cptr->size) && (((Locus *)(cptr->loci[n]))->txtstr[j] == GENOSYMBOL[3]) )
                      n++;
        */
        // these closures correspond to the while-loops in lines 1077-1082 (above block comment)
        // `locus_ix` is the initial locus index
        // `geno_strain` is the index into the `genotype` vec; it corresponds to a given strain
        let find_adj_known = |step: i32, loci: &Vec<Locus>, locus_ix: usize, geno_strain: usize| {
            let mut ix = (locus_ix as i32) + step;
            loop {
                let (geno, _val) = loci[ix as usize].genotype[geno_strain];
                if geno != Genotype::Unk {
                    return ix as usize;
                } else {
                    ix += step;
                }
            }
        };

        /*
           Rewrite this whole loop (and probably the struct, really)

           Might need to use Rc or RefCell, though I have a feeling
           that just thinking it through will be enough.
        */
        /*
        for (chr, loci) in self.chromosomes.iter() {
            for (locus_ix, locus) in loci.iter_mut().enumerate() {
                for (geno_strain, (geno, val)) in locus.genotype.clone().iter().enumerate() {
                    if *geno == Genotype::Unk {
                        let prev_ix = find_adj_known(-1, &loci, locus_ix, geno_strain);
                        let next_ix = find_adj_known(1, &loci, locus_ix, geno_strain);

                        let prev = &loci[prev_ix];
                        let next = &loci[next_ix];

                        /*
                        println!("For {}, {}", locus.name, self.strains[geno_strain]);
                        println!("Prev\n{:?}\n", prev);
                        println!("Next\n{:?}\n\n", next);
                        */
        let ((prev_geno, prev_cm), (next_geno, next_cm)) =
        self.adj_known_loci(chr, locus_ix, geno_strain);

        let rec_1 = (locus.cm() - prev.cm()) / 100.0;
        let rec_2 = (next.cm() - locus.cm()) / 100.0;
        let rec_0 = (next.cm() - prev.cm()) / 100.0;

        let f1 = (1.0 - f64::exp(-2.0 * rec_1)) / 2.0;
        let f2 = (1.0 - f64::exp(-2.0 * rec_2)) / 2.0;
        let f0 = (1.0 - f64::exp(-2.0 * rec_0)) / 2.0;

        // NOTE make sure the parens act the same as the C version!!
        let r_0 = (1.0 - f1) * (1.0 - f2) / (1.0 - f0);
        let r_1 = f1 * (1.0 - f2) / f0;
        let r_2 = f2 * (1.0 - f1) / f0;
        let r_3 = f1 * f2 / (1.0 - f0);

        let (prev_geno, _) = prev.genotype[geno_strain];
        let (next_geno, _) = next.genotype[geno_strain];

        // use Genotype::*;

        let new_genotype = match (prev_geno, next_geno) {
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

        // locus.genotype[geno_strain] = (*geno, new_genotype);
        // *val = new_genotype;

        // match (prev.genotype[geno_strain]._1, next.genotype[geno_strain]._1) {}
        }
        }
        }
        }
         */
    }
}

#[derive(Debug)]
pub struct QTL {
    pub lrs: f64,
    pub additive: f64,
    pub dominance: f64,
    pub locus: Locus,
}

// this could work better as a hashmap, maybe
// or use a hashmap-like interface
// or maybe look like { name: String, values: Vec<(String, f64)>
pub struct Trait {
    pub name: String,
    pub values: Vec<(String, f64)>,
    // pub strains: Vec<String>,
    // pub values: Vec<f64>,
}

// pub struct Traits {

// impl Locus {
//     pub fn new() -> Locus {
//     };
// }

impl QTL {
    pub fn new(locus: Locus, lrs: f64, additive: f64, dominance: f64) -> QTL {
        QTL {
            lrs,
            additive,
            dominance,
            locus,
        }
    }
}

// the PartialEq and PartialOrd definitions are taken from the C version,
// but they don't actually make sense -- two QTLs are the *same* iff they
// have the same LRS? yeah no.
impl PartialEq<QTL> for QTL {
    fn eq(&self, other: &QTL) -> bool {
        self.lrs == other.lrs
    }
}

impl PartialOrd for QTL {
    fn partial_cmp(&self, other: &QTL) -> Option<Ordering> {
        self.lrs.partial_cmp(&other.lrs)
    }
}

fn parse_tab_delim_line(line: &str) -> Vec<String> {
    line.split_terminator('\t')
        .map(|s| String::from(s))
        .collect()
}

#[derive(Debug, Clone)]
pub struct DatasetHeader {
    has_cm: bool,
    pub strains: Vec<String>,
}

impl DatasetHeader {
    pub fn from_line(line: &str) -> Option<DatasetHeader> {
        let header_words = parse_tab_delim_line(&line);

        let has_cm = match header_words.get(3) {
            None => panic!("Dataset header had less than four elements; no strains!"),
            Some(w) => w == "Mb",
        };

        let skip_n = if has_cm { 4 } else { 3 };

        let strains = header_words
            .into_iter()
            .skip(skip_n)
            .map(|s| String::from(s))
            .collect();

        Some(DatasetHeader { has_cm, strains })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header_line() -> String {
        String::from("Chr	Locus	cM	BXD1	BXD2	BXD5	BXD6")
    }

    fn header_line_2() -> String {
        String::from("Chr	Locus	cM	Mb	BXD1	BXD2	BXD5	BXD6")
    }

    fn data_line1() -> String {
        String::from("1	D1Mit1	8.3	B6	B6	D	D")
    }

    #[test]
    fn it_can_parse_header() {
        let header = header_line();
        let result = parse_tab_delim_line(&header);

        assert_eq!(
            vec!["Chr", "Locus", "cM", "BXD1", "BXD2", "BXD5", "BXD6"],
            result
        );

        let parsed = DatasetHeader::from_line(&header).unwrap();

        assert_eq!(false, parsed.has_cm);
        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], parsed.strains);

        let parsed_2 = DatasetHeader::from_line(&header_line_2()).unwrap();

        assert_eq!(true, parsed_2.has_cm);
        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], parsed_2.strains);
    }

    /*
    #[test]
    fn it_can_parse_data_lines() {
        let line = data_line1();
        let result = parse_tab_delim_line(&line);

        assert_eq!(vec!["1", "D1Mit1", "8.3", "B6", "B6", "D", "D"], result);

        let header = DatasetHeader::from_line(&header_line()).unwrap();
        let metadata = Metadata::new("BXD", "B6", "D", "riset");
        let parsed = DatasetLine::from_line(&header, &metadata, &line).unwrap();

        assert_eq!(
            DatasetLine {
                chromosome: String::from("1"),
                locus: String::from("D1Mit1"),
                centi_morgan: 8.3,
                mega_basepair: None,
                genotypes: vec!["B6", "B6", "D", "D"],
            },
            parsed
        );
    }
    */

    #[test]
    fn it_can_parse_metadata_lines() {
        let lines = vec![
            "#@type:intercross",
            "@name:BXD",
            "#abbreviation of maternal or paternal parents",
            "@mat:B6",
        ];

        assert_eq!(Metadata::parse_line(lines[0]), None);
        assert_eq!(Metadata::parse_line(lines[1]), Some(("name", "BXD")));
        assert_eq!(Metadata::parse_line(lines[2]), None);
        assert_eq!(Metadata::parse_line(lines[3]), Some(("mat", "B6")));
    }

    #[test]
    fn it_can_parse_metadata() {
        let lines = vec![
            "#the first three/four columns should be \"Chr, Locus, cM, [Mb]\" (case sensitive)",
        "#please save as Unix format text file.",
        "#comment line always start with a '#'",
        "#type riset or intercross",
        "@type:riset",
        "#@type:intercross",
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
            Metadata::new("BXD", "B6", "D", "riset")
        );
    }

    #[test]
    fn it_can_find_unknown_intervals_in_one_strain() {
        let geno1 = vec![
            (Genotype::Mat, -1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Unk, 99.0),
            (Genotype::Pat, 1.0),
            (Genotype::Pat, 1.0),
        ];
        let geno2 = vec![
            (Genotype::Mat, -1.0),
            (Genotype::Pat, 1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Unk, 99.0),
            (Genotype::Mat, -1.0),
        ];
        let geno3 = vec![
            (Genotype::Pat, 1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Pat, 1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Unk, 99.0),
            (Genotype::Mat, -1.0),
        ];

        let mut state = (None, Vec::new());
        geno1.iter().enumerate().for_each(|(ix, (geno, _))| {
            step_unknown_intervals_mut(&mut state, (ix, *geno));
            // println!("{} - {:?}", ix, state);
        });

        assert_eq!(state.1, vec![1..3]);

        state = (None, Vec::new());
        geno2.iter().enumerate().for_each(|(ix, (geno, _))| {
            step_unknown_intervals_mut(&mut state, (ix, *geno));
            // println!("{} - {:?}", ix, state);
        });

        assert_eq!(state.1, vec![2..4]);

        state = (None, Vec::new());
        geno3.iter().enumerate().for_each(|(ix, (geno, _))| {
            step_unknown_intervals_mut(&mut state, (ix, *geno));
            // println!("{} - {:?}", ix, state);
        });

        assert_eq!(state.1, vec![1..2, 3..5]);
    }

    #[test]
    fn it_can_find_unknown_intervals_in_many_strains() {
        let genos = vec![
            vec![Genotype::Mat, Genotype::Mat, Genotype::Pat],
            vec![Genotype::Unk, Genotype::Pat, Genotype::Unk],
            vec![Genotype::Unk, Genotype::Unk, Genotype::Pat],
            vec![Genotype::Unk, Genotype::Unk, Genotype::Unk],
            vec![Genotype::Pat, Genotype::Mat, Genotype::Unk],
            vec![Genotype::Pat, Genotype::Mat, Genotype::Mat],
        ];

        let strains = 3;
        let mut state = (vec![None; strains], vec![Vec::new(); strains]);

        for (geno_ix, genos_line) in genos.iter().enumerate() {
            step_many_unknown_intervals_mut(&mut state, (geno_ix, &genos_line));
        }

        assert_eq!(state.1, vec![vec![1..4], vec![2..4], vec![1..2, 3..5]]);
    }

    /*
    #[test]
    fn it_can_estimate_unknown_genotypes() {
        // let mut chromosomes = HashMap::new();
        let strains = vec!["S1".to_string(), "S2".to_string(), "S3".to_string()];
        let geno1 = vec![
            (Genotype::Mat, -1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Unk, 99.0),
            (Genotype::Pat, 1.0),
            (Genotype::Pat, 1.0),
        ];
        let geno2 = vec![
            (Genotype::Mat, -1.0),
            (Genotype::Pat, 1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Unk, 99.0),
            (Genotype::Mat, -1.0),
        ];
        let geno3 = vec![
            (Genotype::Pat, 1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Pat, 1.0),
            (Genotype::Unk, 99.0),
            (Genotype::Mat, -1.0),
        ];

        let mk_locus = |name, cm, genotype| Locus {
            name: String::from(name),
            dominance: None,
            genotype,
            centi_morgan: cm,
            mega_basepair: None,
        };

        let loci = vec![
            mk_locus("Mk1", 1.0, geno1),
            mk_locus("Mk2", 3.3, geno2),
            mk_locus("Mk3", 4.1, geno3),
        ];

        // chromosomes.

        // let dataset = Dataset {
        //     metadata: Metadata::new("test", "B", "D", "test"),
        //     header: DatasetHeader {
        //         has_cm: false,
        //         strains: strains.clone(),
        //     },
        //     chromosomes,
        //     strains: strains.clone(),
        // };
    }
    */

}

/*
impl Dataset {
    pub fn read_from_file(self, filename: &str) -> Dataset {
        let path = Path::new(filename);

        let mut file = match File::open(&path) {
            Err(why) => panic!(
                "Couldn't open dataset file {}: {}",
                path.display(),
                why.description()
            ),
            Ok(file) => file,
        };

        println!("reaper: parsing {}", filename);

        let mut buf_reader = BufReader::new(file);
        // let name =

        /*
        The `name`, `maternal`, `paternal`, and `dataset_type` fields are
        read from the dataset.

        The `chromosome` field is populated by creating the chromosomes,
        based on the data read from the file.
         */
// pub name: String,
// pub maternal: String,
// pub paternal: String,
// pub dataset_type: String,
// pub chromosome: Vec<Chromosome>,

// pub size: usize,
// pub nprgy: usize,
// pub prgy: Vec<String>,
// pub parentsf1: bool,
// pub dominance: bool,
// pub mega_basepair: usize,
// pub interval: i32,
//
}
}

*/
