extern crate structopt;

use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;
use structopt::StructOpt;

use qtlreaper::geneobject;
use qtlreaper::regression;

#[derive(StructOpt, Debug)]
#[structopt(name = "qtlreaper")]
struct Opt {
    #[structopt(long = "geno")]
    genotype_file: PathBuf,

    #[structopt(long = "traits")]
    traits_file: PathBuf,

    #[structopt(short = "c", long = "control", long_help = r"control marker name")]
    control: Option<String>,

    #[structopt(
        short = "o",
        long = "output",
        long_help = r"output file",
        default_value = "output.txt"
    )]
    output_file: PathBuf,

    #[structopt(
        short = "n",
        long = "n_permutations",
        long_help = r"number of permutations",
        default_value = "1000"
    )]
    n_permutations: usize,

    #[structopt(
        short = "t",
        long = "threads",
        long_help = r"number of threads to use",
        default_value = "1"
    )]
    threads: usize,
}

fn main() {
    let dataset = geneobject::Dataset::read_file(&PathBuf::from("tests/data/input/BXD_Test.txt"));

    let chromosome = dataset.genome.chromosomes.get("1").unwrap();

    for locus in chromosome {
        let geno: Vec<_> = locus.genotype.iter().map(|(_, g)| g).collect();
        println!("{:?}", geno);
    }

    let new_chr = geneobject::Genome::chromosome_interval(chromosome, 1.0);

    println!("{}", new_chr.len());
    for locus in new_chr {
        let geno: Vec<_> = locus.genotype.iter().map(|(_, g)| g).collect();
        println!("{:?}", geno);
    }
    // println!("{:?}", chromosome);
}

fn main_() {
    let opt = Opt::from_args();

    let dataset = geneobject::Dataset::read_file(&opt.genotype_file);

    let traits = geneobject::Traits::read_file(&opt.traits_file);

    let mut fout = File::create(opt.output_file).unwrap();

    if dataset.dominance {
        fout.write_all(b"ID\tLocus\tChr\tcM\tLRS\tAdditive\tDominance\tpValue\n")
            .unwrap();
    } else {
        fout.write_all(b"ID\tLocus\tChr\tcM\tLRS\tAdditive\tpValue\n")
            .unwrap();
    }

    for (name, values) in traits.traits.iter() {
        let qtls = regression::regression(
            &dataset,
            values,
            &traits.strains,
            opt.control.as_ref().map(|s| &**s),
        );
        let permu = regression::permutation(
            &dataset,
            values,
            &traits.strains,
            opt.n_permutations,
            opt.threads,
        );

        for qtl in qtls.iter() {
            let pvalue = regression::pvalue(qtl.lrs, &permu);

            let line = if dataset.dominance {
                format!(
                    "{}\t{}\t{}\t{:.*}\t{:.*}\t{:.*}\t{:.*}\t{:.*}\n",
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
                    qtl.dominance.unwrap(),
                    3,
                    pvalue
                )
            } else {
                format!(
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
                )
            };

            fout.write_all(line.as_bytes()).unwrap();
        }
    }
}
