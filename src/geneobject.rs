use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

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

pub struct Chromosome {
    pub name: String,
    pub loci: Vec<Locus>,
    pub size: usize,
}

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

fn parse_tab_delim_line_2(line: &str) -> Vec<&str> {
    line.split_terminator('\t')
        // .map(|s| String::from(s))
        .collect()
}

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
    // match header_iter.next() {
    //     None => return None;

    // }
    // None
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

    /*
    let data_line2 = "1	D1Mit231	12	D	B6	B6	D";
    let data_line3 = "1	D1Mit169	14.5	D	B6	B6	D";
    let data_line3 = "1	D1Mit373	17	D	B6	B6	D";
    let data_line3 = "1	D1Mit318	18.5	D	B6	B6	D"; */

    #[test]
    fn parse_header() {
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
    fn parse_data_lines() {
        let line = data_line1();
        let result = parse_tab_delim_line(&line);

        assert_eq!(vec!["1", "D1Mit1", "8.3", "B6", "B6", "D", "D"], result);
    }

}

// fn parse_tab_delimited(filename: &str) ->

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
