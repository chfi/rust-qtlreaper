use super::geneobject::{Genotype, Metadata};
use ndarray::prelude::*;
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::ops::Range;
use std::path::PathBuf;

type GenotypesRaw = Array1<Genotype>;

// struct Chromosome<G, L> {
//     genotypes: G,
//     loci: L,
// }

type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Debug, Clone)]
struct Marker {
    name: String,
    c_m: f64,
    m_bp: Option<f64>,
}

impl Marker {
    pub fn show(&self) -> String {
        let mbp = match self.m_bp {
            None => String::from(""),
            Some(m) => format!("\t-{}", m),
        };
        format!("{}\t- {}{}", self.name, self.c_m, mbp)
    }
}

#[derive(Debug, Clone)]
struct Genome<T> {
    pub chromosomes: Vec<(String, T)>,
}

// struct ChrName(String)

// Comparison of chromosome names (1 < 10, 2 < 10, and 10 < x)
fn chr_cmp(l: &str, r: &str) -> Ordering {
    match (
        l.chars().all(|s| s.is_numeric()),
        r.chars().all(|s| s.is_numeric()),
    ) {
        // both l and r are numbers
        (true, true) => {
            // compare by length then value
            match l.len().cmp(&r.len()) {
                Ordering::Equal => l.cmp(r),
                o => o,
            }
        }
        // only one of l and r is a number, other is text
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        // both l and r are text
        (false, false) => l.cmp(r),
    }
}

impl<T> Genome<T> {
    fn new() -> Self {
        Genome {
            chromosomes: Vec::new(),
        }
    }

    fn get_or_new_mut(&mut self, chr: &str, default: T) -> &mut T {
        if self.chromosomes.iter().find(|(c, _)| c == &chr).is_none() {
            self.chromosomes.push((chr.to_string(), default));
        }

        let (_, t) = self
            .chromosomes
            .iter_mut()
            .find(|(c, _)| c == &chr)
            .expect("impossible");

        t
    }

    fn get(&self, chr: &str) -> Option<&T> {
        self.chromosomes
            .iter()
            .find(|(c, _)| c == &chr)
            .map(|(_, t)| t)
    }

    fn get_mut(&mut self, chr: &str) -> Option<&mut T> {
        self.chromosomes
            .iter_mut()
            .find(|(c, _)| c == &chr)
            .map(|(_, t)| t)
    }

    pub fn map_chrs<U, F>(&self, f: F) -> Genome<U>
    where
        F: Fn(&T) -> U,
    {
        let mut new_chrs = Vec::new();

        for (chr, ls) in self.chromosomes.iter() {
            new_chrs.push((chr.clone(), f(ls)))
        }

        Genome {
            chromosomes: new_chrs,
        }
    }
}

impl<'a, T> Genome<ArrayView2<'a, T>> {
    pub fn get_column(&self, col: usize) -> Vec<ArrayView1<T>> {
        self.chromosomes
            .iter()
            .map(|(_, c)| c.column(col))
            .collect()
    }
}

struct Chromosome {
    genotypes: Array3<f64>,
    loci: Array1<Marker>,
}

fn parse_genotype(meta: &Metadata, geno: &str) -> Result<Genotype> {
    if geno == meta.maternal {
        Ok(Genotype::Mat)
    } else if geno == meta.paternal {
        Ok(Genotype::Pat)
    } else if geno == meta.heterozygous {
        Ok(Genotype::Het)
    } else if geno == meta.unknown {
        Ok(Genotype::Unk)
    } else {
        Err(format!("Failed to parse genotype: {}\n", geno))?
    }
}

struct Dataset<T> {
    pub metadata: Metadata,
    pub genome: Genome<T>,
    pub gmap: Genome<Array1<Marker>>,
    pub strains: Vec<String>,
    pub has_mb: bool,
}

impl<T> Dataset<T> {
    pub fn strain_ix(&self, strain: &str) -> Option<usize> {
        self.strains
            .iter()
            .enumerate()
            .find(|(_, s)| *s == strain)
            .map(|(ix, _)| ix)
    }
}

// impl<'a, T> Dataset<ArrayView2<'a, T>> {
//     pub fn get_column(&self, col: usize) -> Option<Genome<ArrayView1<'a, T>>> {
//         if self.genome.chromosomes[0].1.cols() <= col {
//             Some(Genome {
//                 chromosomes: self
//                     .genome
//                     .chromosomes
//                     .iter()
//                     .map(|(_, chr)| chr.column(col))
//                     .collect(),
//             })
//         // Some(self.genome.column(col))
//         } else {
//             None
//         }
//     }
// }

impl Dataset<Array2<Genotype>> {
    fn to_genoprob(&self) -> Dataset<Array3<f64>> {
        // find the unknown genotype ranges for each strain

        // find the closest known loci for each unknown range for each strain

        // generate the new genome with probabilities

        let genome = Genome::new();
        Dataset {
            // self..
            genome,
            metadata: self.metadata.clone(),
            strains: self.strains.clone(),
            gmap: self.gmap.clone(),
            has_mb: self.has_mb,
        }
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

    fn read_line(
        metadata: &Metadata,
        has_mb: bool,
        line: String,
    ) -> Result<(String, Marker, Vec<Genotype>)> {
        // let line_ = line?;
        let words: Vec<_> = line.split_terminator('\t').collect();

        let chr = words[0].to_string();
        let name = words[1].to_string();
        let c_m = words[2].parse::<f64>()?;

        let m_bp = match (has_mb, words[3].parse::<f64>()) {
            (true, Ok(m)) => Some(m),
            (_, _) => None,
        };

        let marker = Marker { name, c_m, m_bp };

        let range = if has_mb { 4.. } else { 3.. };

        let genotype = words[range.clone()]
            .iter()
            .map(|g| metadata.parse_genotype(g))
            .collect();

        Ok((chr, marker, genotype))
    }

    pub fn read_file(path: &PathBuf) -> Result<Self> {
        let f =
            File::open(path).expect(&format!("Error opening file {:?}", path));

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
                    let ll = l.expect("Error parsing dataset");
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

        // let mut chrs: BTreeMap<String, Vec< = BTreeMap::new();

        let mut genome_vec: Genome<Vec<(Marker, Vec<Genotype>)>> =
            Genome::new();

        for line_ in lines {
            let line = line_?;
            let (chr, marker, geno) = Self::read_line(&metadata, has_mb, line)?;

            let a_chr = genome_vec.get_or_new_mut(&chr, Vec::new());
            a_chr.push((marker, geno))
        }

        let mut gmap: Genome<Array1<Marker>> = Genome::new();

        let mut genome: Genome<Array2<Genotype>> = Genome::new();

        for (chr, loci) in genome_vec.chromosomes.into_iter() {
            let rows = loci.len();
            let (m, g): (Vec<_>, Vec<_>) = loci.into_iter().unzip();

            let chr_gmap = Array::from_vec(m);
            let geno_cols = g[0].len();
            let chr_geno = Array::from_shape_vec((rows, geno_cols), g.concat())
                .expect("Error when reshaping genotype matrix");

            gmap.get_or_new_mut(&chr, chr_gmap);
            genome.get_or_new_mut(&chr, chr_geno);
        }

        Ok(Dataset {
            metadata,
            strains,
            genome,
            gmap,
            has_mb,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_genotypes() {
        let path = "./tests/data/input/BXD.txt";

        let ds = Dataset::read_file(&PathBuf::from(path));

        let dataset = match ds {
            Err(err) => panic!("Error: {}", err),
            Ok(d) => d,
        };

        println!("Metadata\n{:?}\n", dataset.metadata);

        let gmap_rows: usize =
            dataset.gmap.chromosomes.iter().map(|(_, l)| l.len()).sum();

        println!("gmap rows: {}\n", gmap_rows);

        println!("strains: {:?}\n", dataset.strains);

        assert_eq!(true, false);
    }
}
