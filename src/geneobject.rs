use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug)]
pub struct Locus {
    pub name: String,
    pub chromosome: Chromosome,
    pub genotype: f64,
    pub dominance: f64,
    pub text: String,
    pub size: usize,
    pub centi_morgan: f64,
    pub mega_basepair: f64,
}

#[derive(Debug)]
pub struct Chromosome {
    pub name: String,
    pub loci: Vec<Locus>,
    pub size: usize,
}

#[derive(Debug)]
pub struct Metadata {
    pub name: String,
    pub maternal: String,
    pub paternal: String,
    pub heterozygous: String, // defaults to "H"
    pub unknown: String,      // defaults to "U"
}

impl Metadata {
    fn new(name: &str, maternal: &str, paternal: &str) -> Metadata {
        Metadata {
            name: String::from(name),
            maternal: String::from(maternal),
            paternal: String::from(paternal),
            heterozygous: String::from("H"),
            unknown: String::from("U"),
        }
    }
}

#[derive(Debug)]
pub struct Dataset {
    pub name: String,
    pub maternal: String,
    pub paternal: String,
    pub dataset_type: String,
    pub chromosome: Vec<Chromosome>,
    // pub size: usize,
    // pub nprgy: usize,
    // pub prgy: Vec<String>,
    // pub parentsf1: bool,
    // pub dominance: bool,
    // pub mega_basepair: usize,
    // pub interval: i32,
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

// fn parse_metadata_line(line: &str) -> Option<(String, String)> {
fn parse_metadata_line(line: &str) -> Option<(&str, &str)> {
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
fn parse_metadata_lines(lines: Vec<&str>) -> Metadata {
    let mut name: Option<String> = None;
    let mut mat: Option<String> = None;
    let mut pat: Option<String> = None;
    let mut het = String::from("H");
    let mut unk = String::from("U");

    for line in lines.iter() {
        if let Some((n, v)) = parse_metadata_line(line) {
            match n {
                "name" => name = Some(String::from(v)),
                "mat" => mat = Some(String::from(v)),
                "pat" => pat = Some(String::from(v)),
                "het" => het = String::from(v),
                "unk" => unk = String::from(v),
                _ => (),
            }
        }
    }

    if name == None || mat == None || pat == None {
        panic!(
            "Required metadata was not provided!\nname = {:?}\nmat = {:?}\npat = {:?}",
            name, mat, pat
        );
    }

    Metadata {
        name: name.unwrap(),
        maternal: mat.unwrap(),
        paternal: pat.unwrap(),
        heterozygous: het,
        unknown: unk,
    }
}

#[derive(Debug)]
struct DatasetHeader {
    has_cm: bool,
    strains: Vec<String>,
}

fn parse_header_line(line: &str) -> Option<DatasetHeader> {
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

#[derive(Debug, PartialEq)]
struct DatasetLine<'a> {
    chromosome: String,
    locus: String,
    centi_morgan: f64,
    mega_basepair: Option<f64>,
    genotypes: Vec<&'a str>,
}

fn parse_data_line<'a>(
    header: &DatasetHeader,
    metadata: &'a Metadata,
    line: &str,
) -> Option<DatasetLine<'a>> {
    let words = parse_tab_delim_line(&line);

    let chromosome = words.get(0).unwrap().clone();
    let locus = words.get(1).unwrap().clone();
    let centi_morgan = words.get(2).unwrap().parse::<f64>().unwrap();

    let skip_n = if header.has_cm { 4 } else { 3 };
    let mega_basepair = if header.has_cm {
        words.get(3).map(|s| s.parse::<f64>().unwrap())
    } else {
        None
    };

    let geno_ref = |s| {
        let mat = metadata.maternal.as_str();
        let pat = metadata.paternal.as_str();
        if s == mat {
            mat
        } else if s == pat {
            pat
        } else {
            panic!(
                "Locus {} has a genotype that's neither maternal nor paternal!",
                locus
            );
        }
    };

    let genotypes = words.into_iter().skip(skip_n).map(geno_ref).collect();

    Some(DatasetLine {
        chromosome,
        locus,
        centi_morgan,
        mega_basepair,
        genotypes,
    })
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

        let parsed = parse_header_line(&header).unwrap();

        assert_eq!(false, parsed.has_cm);

        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], parsed.strains);

        let header_2 = header_line_2();
        let parsed_2 = parse_header_line(&header_2).unwrap();

        assert_eq!(true, parsed_2.has_cm);

        assert_eq!(vec!["BXD1", "BXD2", "BXD5", "BXD6"], parsed_2.strains);
    }

    #[test]
    fn it_can_parse_data_lines() {
        let line = data_line1();
        let result = parse_tab_delim_line(&line);

        assert_eq!(vec!["1", "D1Mit1", "8.3", "B6", "B6", "D", "D"], result);

        let header = parse_header_line(&header_line()).unwrap();
        let metadata = Metadata {
            name: String::from("BXD"),
            maternal: String::from("B6"),
            paternal: String::from("D"),
            heterozygous: String::from("H"),
            unknown: String::from("U"),
        };
        let parsed = parse_data_line(&header, &metadata, &line).unwrap();

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

    #[test]
    fn it_can_parse_metadata_lines() {
        let lines = vec![
            "#@type:intercross",
            "@name:BXD",
            "#abbreviation of maternal or paternal parents",
            "@mat:B6",
        ];

        assert_eq!(parse_metadata_line(lines[0]), None);
        assert_eq!(parse_metadata_line(lines[1]), Some(("name", "BXD")));
        assert_eq!(parse_metadata_line(lines[2]), None);
        assert_eq!(parse_metadata_line(lines[3]), Some(("mat", "B6")));
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

        println!("{:?}", parse_metadata_lines(lines));
    }

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
