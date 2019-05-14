use std::cmp::Ordering;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use std::process;

use qtlreaper::geneobject;

/*
The original code has a lot of reference cycles and stuff,
with each locus referring back to its chromosome, and
a chromosome having a list of loci; then the QTLs also
have a locus, of course.

the locus doesn't need to know about the chromosome, I bet. the chr
can be passed, if and when needed, since each chromosome contains its
loci.

buuuut I still don't know how the actual strain data is stored/represented

the Locus has the genotype and dominance as *floats*, not vectors.

there's also the "text" field... wthf

in the c code it's the `txtstr` field, which is filled from an array of
chars containing "BDHU", probably referring to "B" and "D" as in "BXD",
and "H", "U" as "Heterozygous" and "Unknown"

soooo... the `txtstr` field is filled, using the metadata which says
what string corresponds to the maternal genotype, etc., with those
known bytes.


looking at the parsing and regression functions in dataset.c, it looks like
the "genotype" field, which is... an array of floats? is set, for each
strain, to some number based on the "previous" and the "next" genotypes.
I don't know if the "previous" and the "next" depend on the strain, or
if the dependency is on the next and previous loci.

it's the next and previous loci, which makes sense given what I know
about the QTL mapping algorithm.

*/

/*
usage:
qtlreaper <genotype> <traits> <strains>
*/

pub struct Config {
    pub genotype_file: String,
    // pub traits_file: String,
    // pub strains_file: String,
}

// pub enum Verbosity {
//     Verbal,
//     Normal,
//     Quiet
// }

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 2 {
            return Err("Not enough arguments");
        }

        let genotype_file = args[1].clone();
        // let traits_file = args[2].clone();
        // let strains_file = args[3].clone();

        Ok(Config {
            genotype_file,
            // traits_file,
            // strains_file,
        })
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let config = Config::new(&args).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    let dataset = geneobject::Dataset::read_file(&config.genotype_file);

    println!("Parsed genotype file {}", config.genotype_file);

    println!("Dataset has strains: ");

    dataset.strains().iter().for_each(|s| println!("{}", s));

    geneobject::Traits::read_file("examples/data/input/trait.txt");
    // for strain in dataset.strains.iter() {

    // }

    /*
    println!("--------------");
    for (chr, loci) in dataset.chromosomes() {
        println!("Chromosome {}", chr);
        for locus in loci {
            println!("name: {}\n{:?}\n-------", locus.name, locus.genotype);
        }
        println!("----------------");
    }
    */
}
