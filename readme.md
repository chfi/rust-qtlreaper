#### rust-qtlreaper


A reimplementation of [genenetwork/QTLReaper] in Rust.

Build with `cargo build --release`, output in `target/release`.

##### Usage

```
qtlreaper 0.1.2
Christian Fischer <christian@chfi.se>

USAGE:
    qtlreaper [FLAGS] [OPTIONS] --geno <genotype_file> --traits <traits_file>

FLAGS:
    -b, --bootstrap
            run bootstrap

    -h, --help
            Prints help information

        --json
            output in JSON instead of tab-delimited

    -V, --version
            Prints version information


OPTIONS:
        --bootstrap_output <bootstrap_output>
            bootstrap output file [default: bootstrap.txt]

    -c, --control <control>
            control marker name

        --geno <genotype_file>


        --n_bootstrap <n_bootstrap>
            bootstrap count [default: 1000]

    -n, --n_permutations <n_permutations>
            number of permutations [default: 1000]

    -o, --main_output <output_file>
            p-values output file [default: output.txt]

    -t, --threads <threads>
            number of threads to use [default: 1]

        --traits <traits_file>
```


###### Example

```
qtlreaper --geno tests/data/input/BXD.txt --traits tests/data/input/trait.txt -o output
```


####### TODO

- [X] `addinterval`
- [ ] `addparentsf1`
- [X] `bootstrap`
- [X] configurable number of permutations
- [ ] configurable LRS threshold on permutations
- [ ] abort permutations when LRS threshold reached
