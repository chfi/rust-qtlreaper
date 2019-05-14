use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::ops::Range;

#[derive(Debug, PartialEq)]
struct Metadata {
    name: String,
    maternal: String,
    paternal: String,
    dataset_type: String, // "intercross" or "riset"
    heterozygous: String, // defaults to "H"
    unknown: String,      // defaults to "U"
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

    fn parse_dominance(&self, geno: &str) -> f64 {
        if geno == self.maternal.as_str() {
            0.0
        } else if geno == self.paternal.as_str() {
            0.0
        } else if geno == self.heterozygous.as_str() {
            1.0
        } else if geno == self.unknown.as_str() {
            1.0
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
    fn from_lines(lines: Vec<&str>) -> Metadata {
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

#[derive(Debug, Clone)]
pub struct DatasetHeader {
    has_mb: bool,
    pub strains: Vec<String>,
}

impl DatasetHeader {
    fn from_line(line: &str) -> Option<DatasetHeader> {
        let header_words = parse_tab_delim_line(&line);

        let has_mb = match header_words.get(3) {
            None => panic!("Dataset header had less than four elements; no strains!"),
            Some(w) => w == "Mb",
        };

        let skip_n = if has_mb { 4 } else { 3 };

        let strains = header_words
            .into_iter()
            .skip(skip_n)
            .map(|s| String::from(s))
            .collect();

        Some(DatasetHeader { has_mb, strains })
    }
}

#[derive(Debug, PartialEq)]
pub struct Locus {
    pub name: String,
    dominance: Option<Vec<f64>>,
    pub genotype: Vec<(Genotype, f64)>,
    centi_morgan: f64,
    mega_basepair: Option<f64>,
}

/// UnknownIntervals holds a list of ranges of unknown genotypes, per strain
struct UnknownIntervals(Vec<Vec<Range<usize>>>);

impl Locus {
    // corresponds to lines 950-1044 in dataset.c
    fn parse_line(
        metadata: &Metadata,
        header: &DatasetHeader,
        dominance: bool,
        line: &str,
    ) -> (String, Locus) {
        // Example locus is: "1	D1Mit1	8.3	B6	B6	D	D"
        // where the first three columns are chromosome, name, cM;
        // remaining columns are the genotypes

        let words: Vec<_> = line.split_terminator('\t').collect();

        // let chromosome = String::from(words.next().unwrap());
        let chromosome = String::from(words[0]);
        let name = String::from(words[1]);
        let centi_morgan = words[2].parse::<f64>().unwrap();

        let range = if header.has_mb { 4.. } else { 3.. };

        let genotype = words[range.clone()]
            .iter()
            .map(|g| metadata.parse_genotype(g))
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
                name,
                centi_morgan,
                genotype,
                dominance,
                mega_basepair: None,
            },
        )
    }

    /// Steps through a list of genotypes per strain, building up a list of ranges of missing data for each strain
    fn step_many_unknown_intervals_mut(
        state: &mut (Vec<Option<usize>>, Vec<Vec<Range<usize>>>),
        next: (usize, &[(Genotype, f64)]),
    ) {
        let (ix, genotype) = next;

        for (strain_ix, (geno, _)) in genotype.iter().enumerate() {
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

    fn find_unknown_intervals(loci: &[Locus]) -> UnknownIntervals {
        let n_strains = loci.first().unwrap().genotype.len();
        let mut state = (vec![None; n_strains], vec![Vec::new(); n_strains]);

        for (locus_ix, locus) in loci.iter().enumerate() {
            Self::step_many_unknown_intervals_mut(&mut state, (locus_ix, &locus.genotype))
        }

        UnknownIntervals(state.1)
    }

    fn estimate_unknown_genotypes(
        dominance: bool,
        loci: &mut [Locus],
        intervals: UnknownIntervals,
    ) {
        for (strain_ix, strain) in intervals.0.iter().enumerate() {
            for range in strain {
                for locus_ix in range.clone() {
                    let prev = &loci[range.start - 1];
                    let next = &loci[range.end];
                    let locus = &loci[locus_ix];
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

                    let (prev_geno, _) = prev.genotype[strain_ix];
                    let (next_geno, _) = next.genotype[strain_ix];

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

                    if dominance {
                        let new_dominance = match (prev_geno, next_geno) {
                            (Genotype::Mat, Genotype::Mat) => 2.0 * r_0 * r_3,
                            (Genotype::Het, Genotype::Mat) => r_1 * (r_2 + r_3),
                            (Genotype::Pat, Genotype::Mat) => 2.0 * r_1 * r_2,
                            (Genotype::Mat, Genotype::Het) => r_1 * r_0 + r_2 * r_3,
                            (Genotype::Het, Genotype::Het) => {
                                let w = ((1.0 - f0) * (1.0 - f0)) / (1.0 - 2.0 * f0 * (1.0 - f0));
                                1.0 - 2.0 * w * r_0 * r_3 - 2.0 * (1.0 - w) * r_1 * r_2
                            }
                            (Genotype::Pat, Genotype::Het) => r_0 * r_1 + r_2 * r_3,
                            (Genotype::Mat, Genotype::Pat) => 2.0 * r_1 * r_2,
                            (Genotype::Het, Genotype::Pat) => r_1 * (r_2 + r_3),
                            (Genotype::Pat, Genotype::Pat) => 2.0 * r_1 * r_3,
                            _ => panic!("Genotype was unknown when it shouldn't be!"),
                        };

                        if let Some(d) = &mut loci[locus_ix].dominance {
                            d[strain_ix] = new_dominance;
                        }
                    }

                    loci[locus_ix].genotype[strain_ix].1 = new_genotype
                }
            }
        }
    }

    pub fn cm(&self) -> f64 {
        self.centi_morgan
    }
}

#[derive(Debug)]
pub struct Chromosome {
    pub name: String,
    pub loci: Vec<Locus>,
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
pub enum Genotype {
    Mat,
    Pat,
    Het,
    Unk,
}

#[derive(Debug)]
pub struct Dataset {
    metadata: Metadata,
    header: DatasetHeader,
    chromosomes: HashMap<String, Vec<Locus>>,
    strains: Vec<String>,
    dominance: bool, // true if dataset type is "intercross"
}

impl Dataset {
    fn new(metadata: Metadata, header: DatasetHeader, strains: Vec<String>) -> Dataset {
        let dominance = metadata.dataset_type == String::from("intercross");
        Dataset {
            metadata,
            header,
            strains,
            chromosomes: HashMap::new(),
            dominance,
        }
    }

    pub fn strains(&self) -> &[String] {
        &self.strains
    }

    pub fn read_file(path: &str) -> Dataset {
        let f = File::open(path).expect(&format!("Error opening file {}", path));

        let reader = BufReader::new(f);
        let mut lines = reader.lines();
        let mut header = None;

        let mut metadata_lines = vec![];

        loop {
            match lines.next() {
                None => panic!("Reached end of file before parsing dataset header"),
                Some(l) => {
                    let ll = l.unwrap();
                    if ll.starts_with("Chr	Locus	cM") {
                        header = DatasetHeader::from_line(&ll);
                        break;
                    } else {
                        metadata_lines.push(ll);
                    }
                }
            }
        }

        let metadata = Metadata::from_lines(metadata_lines.iter().map(|s| s.as_str()).collect());

        let strains = header.clone().unwrap().strains.clone();

        let mut dataset = Dataset::new(metadata, header.unwrap().clone(), strains);

        for line in lines {
            dataset.parse_locus(&line.unwrap());
        }
        dataset.estimate_unknown();

        dataset
    }

    pub fn chromosomes(&self) -> &HashMap<String, Vec<Locus>> {
        &self.chromosomes
    }

    // parse a line with a locus, adding it to the corresponding
    // chromosome in the dataset. NOTE: we assume the input data is
    // sorted by chromosome and marker position!
    fn parse_locus(&mut self, line: &str) {
        let (chr, locus) = Locus::parse_line(&self.metadata, &self.header, self.dominance, line);

        let loci = self.chromosomes.entry(chr).or_insert_with(|| Vec::new());

        loci.push(locus);
    }

    // Corresponds to lines 1071-1152 in dataset.c
    fn estimate_unknown(&mut self) {
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

        for (_chr, loci) in self.chromosomes.iter_mut() {
            // then, for each chromosome, construct the intervals of
            // unknown genotypes
            let unk = Locus::find_unknown_intervals(loci);

            // ... and use those intervals to estimate the
            // missing genotypes
            Locus::estimate_unknown_genotypes(self.dominance, loci, unk);
        }
    }
}

#[derive(Debug)]
pub struct QTL {
    pub lrs: f64,
    pub additive: f64,
    pub dominance: Option<f64>,
    pub locus: Locus,
}

impl QTL {
    pub fn new(locus: Locus, lrs: f64, additive: f64, dominance: Option<f64>) -> QTL {
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

pub struct Traits {
    strains: Vec<String>,
    traits: HashMap<String, Vec<f64>>,
}

impl Traits {
    pub fn read_file(path: &str) -> Traits {
        let f = File::open(path).expect(&format!("Error opening traits file {}", path));

        let reader = BufReader::new(f);
        let mut lines = reader.lines();

        let strains = match lines.next() {
            None => panic!("Reached end of file before parsing traits header"),
            Some(l) => {
                let ll = l.unwrap();
                if ll.starts_with("Trait") {
                    ll.split_terminator('\t')
                        .skip(1)
                        .map(|s| s.to_string())
                        .collect()
                } else {
                    panic!("Traits file did not begin with \"Trait\", aborting");
                }
            }
        };

        let mut traits = HashMap::new();

        for line in lines {
            let ll = line.unwrap();
            let mut words = ll.split_terminator('\t');
            let key = words.next().unwrap().to_string();
            let values = words.map(|s| s.parse::<f64>().unwrap()).collect();
            traits.insert(key, values);
        }

        println!("parsed strains: {:?}", strains);
        println!("parsed traits: {:?}", traits);

        Traits { strains, traits }
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

    #[test]
    fn it_can_parse_header() {
        let header = header_line();
        let result = parse_tab_delim_line(&header);

        assert_eq!(
            vec!["Chr", "Locus", "cM", "BXD1", "BXD2", "BXD5", "BXD6"],
            result
        );

        let parsed = DatasetHeader::from_line(&header).unwrap();

        assert_eq!(false, parsed.has_mb);
        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], parsed.strains);

        let parsed_2 = DatasetHeader::from_line(&header_line_2()).unwrap();

        assert_eq!(true, parsed_2.has_mb);
        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], parsed_2.strains);
    }

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
    fn it_can_find_unknown_intervals_in_many_strains() {
        let genos = vec![
            vec![
                (Genotype::Mat, -1.0),
                (Genotype::Mat, -1.0),
                (Genotype::Pat, 1.0),
            ],
            vec![
                (Genotype::Unk, 99.0),
                (Genotype::Pat, 1.0),
                (Genotype::Unk, 99.0),
            ],
            vec![
                (Genotype::Unk, 99.0),
                (Genotype::Unk, 99.0),
                (Genotype::Pat, 1.0),
            ],
            vec![
                (Genotype::Unk, 99.0),
                (Genotype::Unk, 99.0),
                (Genotype::Unk, 99.0),
            ],
            vec![
                (Genotype::Pat, 1.0),
                (Genotype::Mat, -1.0),
                (Genotype::Unk, 99.0),
            ],
            vec![
                (Genotype::Pat, 1.0),
                (Genotype::Mat, -1.0),
                (Genotype::Mat, -1.0),
            ],
        ];

        let strains = 3;
        let mut state = (vec![None; strains], vec![Vec::new(); strains]);

        for (geno_ix, genos_line) in genos.iter().enumerate() {
            Locus::step_many_unknown_intervals_mut(&mut state, (geno_ix, &genos_line));
        }

        assert_eq!(state.1, vec![vec![1..4], vec![2..4], vec![1..2, 3..5]]);
    }

    #[test]
    fn it_can_estimate_unknown_genotypes() {
        // let mut chromosomes = HashMap::new();
        let strains = vec!["S1".to_string(), "S2".to_string(), "S3".to_string()];

        let genos = vec![
            vec![
                (Genotype::Mat, -1.0),
                (Genotype::Mat, -1.0),
                (Genotype::Pat, 1.0),
            ],
            vec![
                (Genotype::Unk, 99.0),
                (Genotype::Pat, 1.0),
                (Genotype::Unk, 99.0),
            ],
            vec![
                (Genotype::Unk, 99.0),
                (Genotype::Unk, 99.0),
                (Genotype::Pat, 1.0),
            ],
            vec![
                (Genotype::Unk, 99.0),
                (Genotype::Unk, 99.0),
                (Genotype::Unk, 99.0),
            ],
            vec![
                (Genotype::Pat, 1.0),
                (Genotype::Mat, -1.0),
                (Genotype::Unk, 99.0),
            ],
            vec![
                (Genotype::Pat, 1.0),
                (Genotype::Mat, -1.0),
                (Genotype::Mat, -1.0),
            ],
        ];

        let mk_locus = |name, cm, genotype| Locus {
            name: String::from(name),
            dominance: None,
            genotype,
            centi_morgan: cm,
            mega_basepair: None,
        };

        let loci_new = vec![
            mk_locus(
                "Mk1",
                10.0,
                vec![
                    (Genotype::Mat, -1.0),
                    (Genotype::Mat, -1.0),
                    (Genotype::Pat, 1.0),
                ],
            ),
            mk_locus(
                "Mk2",
                30.3,
                vec![
                    (Genotype::Unk, -0.18523506128077272),
                    (Genotype::Pat, 1.0),
                    (Genotype::Unk, 0.9616255554798838),
                ],
            ),
            mk_locus(
                "Mk3",
                40.1,
                vec![
                    (Genotype::Unk, 0.18906668494617707),
                    (Genotype::Unk, 0.3421367343627405),
                    (Genotype::Pat, 1.0),
                ],
            ),
            mk_locus(
                "Mk4",
                50.2,
                vec![
                    (Genotype::Unk, 0.5826065314914579),
                    (Genotype::Unk, -0.3223330030526561),
                    (Genotype::Unk, 0.3223330030526552),
                ],
            ),
            mk_locus(
                "Mk5",
                60.3,
                vec![
                    (Genotype::Pat, 1.0),
                    (Genotype::Mat, -1.0),
                    (Genotype::Unk, -0.34213673436274084),
                ],
            ),
            mk_locus(
                "Mk6",
                70.1,
                vec![
                    (Genotype::Pat, 1.0),
                    (Genotype::Mat, -1.0),
                    (Genotype::Mat, -1.0),
                ],
            ),
        ];

        let mut loci = vec![
            mk_locus("Mk1", 10.0, genos[0].clone()),
            mk_locus("Mk2", 30.3, genos[1].clone()),
            mk_locus("Mk3", 40.1, genos[2].clone()),
            mk_locus("Mk4", 50.2, genos[3].clone()),
            mk_locus("Mk5", 60.3, genos[4].clone()),
            mk_locus("Mk6", 70.1, genos[5].clone()),
        ];

        let unk = Locus::find_unknown_intervals(&loci);

        Locus::estimate_unknown_genotypes(false, &mut loci, unk);

        assert_eq!(loci, loci_new);
    }

}
