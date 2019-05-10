use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

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

fn main() {
    let dataset = geneobject::Dataset::read_file("examples/data/input/BXD_Test.txt");

    println!("--------------");
    for (chr, loci) in dataset.chromosomes() {
        println!("Chromosome {}", chr);
        for locus in loci {
            println!("name: {}\n{:?}\n-------", locus.name, locus.genotype);
        }
        println!("----------------");
    }
}
