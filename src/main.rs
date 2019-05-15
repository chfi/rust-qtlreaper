use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::process;

use qtlreaper::geneobject;
use qtlreaper::regression;

pub struct Config {
    pub genotype_file: String,
    pub traits_file: String,
    // pub strains_file: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 3 {
            return Err("Not enough arguments");
        }

        let genotype_file = args[1].clone();
        let traits_file = args[2].clone();
        // let strains_file = args[3].clone();

        Ok(Config {
            genotype_file,
            traits_file,
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

    println!("dataset has strains: ");

    dataset.strains().iter().for_each(|s| println!("{}", s));

    println!("dominance: {}", dataset.dominance);

    let traits = geneobject::Traits::read_file(&config.traits_file);

    let mut fout = File::create("output.txt").unwrap();

    fout.write(b"ID\tLocus\tChr\tcM\tLRS\tAdditive\tpValue\n")
        .unwrap();

    // println!("{:?}", dataset.genome.chromosomes);

    for (name, values) in traits.traits.iter() {
        let qtls = regression::regression(&dataset, values, &traits.strains);
        let permu = regression::permutation(&dataset, values, &traits.strains);

        for qtl in qtls.iter() {
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
                // 3,
                // qtl.dominance.unwrap(),
                3,
                pvalue
            );

            fout.write(line.as_bytes()).unwrap();
        }
    }
}
