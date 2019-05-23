extern crate structopt;

use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;
use structopt::StructOpt;

use qtlreaper::geneobject::{Dataset, Traits};
use qtlreaper::regression;

use serde_json;

#[derive(StructOpt, Debug)]
#[structopt(name = "qtlreaper")]
struct Opt {
    #[structopt(long = "geno")]
    genotype_file: PathBuf,

    #[structopt(long = "traits")]
    traits_file: PathBuf,

    #[structopt(
        short = "c",
        long = "control",
        long_help = r"control marker name"
    )]
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

    #[structopt(
        long = "json",
        long_help = r"output in JSON instead of tab-delimited"
    )]
    output_json: bool,
}

fn format_header(dataset: &Dataset) -> String {
    let mut start = String::from("ID\tLocus\tChr\tcM");

    if dataset.has_mb() {
        start += "\tMb";
    }
    start += "\tLRS\tAdditive";
    if dataset.dominance {
        start += "\tDominance";
    }

    start + "\tpValue\n"
}

fn main() {
    let opt = Opt::from_args();

    let dataset = Dataset::read_file(&opt.genotype_file);

    let traits = Traits::read_file(&opt.traits_file);

    let mut fout = File::create(opt.output_file).unwrap();

    if !opt.output_json {
        fout.write_all(format_header(&dataset).as_bytes()).unwrap();
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

        /*
        let bootstrap = regression::bootstrap(
            &dataset,
            values,
            &traits.strains,
            None,
            1000,
        );
        */

        if opt.output_json {
            for qtl in qtls.iter() {
                fout.write_all(serde_json::to_string(qtl).unwrap().as_bytes())
                    .unwrap();
            }
        } else {
            for qtl in qtls.iter() {
                let pvalue = regression::pvalue(qtl.lrs, &permu);

                let line = format!("{}\t{}\t{:.*}\n", name, qtl, 3, pvalue);

                fout.write_all(line.as_bytes()).unwrap();
            }
        }
    }
}
