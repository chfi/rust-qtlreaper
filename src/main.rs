use std::cmp::Ordering;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use std::process;

use qtlreaper::geneobject;
use qtlreaper::regression;

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

    for (chr, loci) in dataset.chromosomes().iter() {
        for locus in loci.iter() {
            let genotypes = locus.genotypes();
            println!("Locus(\"{}\", {:?})", locus.marker.name, genotypes);
        }
    }

    println!("Parsed genotype file {}", config.genotype_file);

    println!("Dataset has strains: ");

    dataset.strains().iter().for_each(|s| println!("{}", s));

    let traits = geneobject::Traits::read_file("examples/data/input/trait_one.txt");
    // for strain in dataset.strains.iter() {

    // }

    /*
    for qtl in qtlresults:
          fout_qtl.write("%s\t%s\t%s\t%2.3f\t%2.3f\t%2.3f\t\n" % (traitName,
              qtl.locus.name, qtl.locus.chr, qtl.locus.cM, qtl.lrs, qtl.additive))
    */

    let mut fout = File::create("output.txt").unwrap();

    fout.write(b"ID\tLocus\tChr\tcM\tLRS\tAdditive\tpValue\n");

    for (name, values) in traits.traits.iter() {
        let qtls = regression::regression(&dataset, values, &traits.strains);
        let permu = regression::permutation(&dataset, values);

        // for p in permu.iter() {
        //     println!("{}", p);
        // }

        for qtl in qtls.iter() {
            // println!(
            //     "{}\t{}\t{}\t{}\t{}\t{}",
            //     name,
            //     qtl.marker.name,
            //     qtl.marker.chromosome,
            //     qtl.marker.centi_morgan,
            //     qtl.lrs,
            //     qtl.additive
            // )

            let pvalue = regression::pvalue(qtl.lrs, &permu);
            let line = format!(
                "{}\t{}\t{}\t{:.*}\t{:.*}\t{:.*}\t{:.*}\n",
                name,
                qtl.marker.name,
                qtl.marker.chromosome,
                3,
                qtl.marker.centi_morgan,
                3,
                qtl.lrs,
                3,
                qtl.additive,
                3,
                pvalue
            );

            fout.write(line.as_bytes()).unwrap();
        }
    }
    // let
    // regression::regression(dataset,

    // let vec_in = vec![1, 2, 3, 4, 5, 6];
    // let vec_perm = regression::permuted(&vec_in);

    // println!("before:");
    // println!("{:?}\n", vec_in);

    // println!("after:");
    // println!("{:?}\n", vec_perm);

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
