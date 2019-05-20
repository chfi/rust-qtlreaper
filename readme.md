#### rust-qtlreaper


A reimplementation of [genenetwork/QTLReaper] in Rust.

Build with `cargo build --release`, output in `target/release`.

##### Usage

```
qtlreaper 0.1.1
Christian Fischer <christian@chfi.se>

USAGE:
    qtlreaper [OPTIONS] --geno <genotype_file> --traits <traits_file>

FLAGS:
    -h, --help
            Prints help information

    -V, --version
            Prints version information


OPTIONS:
    -c, --control <control>
            control marker name

        --geno <genotype_file>


    -n, --n_permutations <n_permutations>
            number of permutations [default: 1000]

    -o, --output <output_file>
            output file [default: output.txt]

    -t, --threads <threads>
            number of threads to use [default: 1]

        --traits <traits_file>
```


###### Example

```
qtlreaper --geno tests/data/input/BXD.txt --traits tests/data/input/trait.txt -o output
```


####### TODO

- [ ] `addinterval`
- [ ] `addparentsf1`
- [ ] `bootstrap`
- [X] configurable number of permutations
- [ ] configurable LRS threshold on permutations
- [ ] abort permutations when LRS threshold reached
