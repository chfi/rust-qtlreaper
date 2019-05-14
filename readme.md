#### rust-qtlreaper

A reimplementation of [genenetwork/QTLReaper] in Rust.

##### TODO

- [X] Parse input files
- [ ] Implement the simplest type of regression used
- [ ] Test against QTLReaper output


###### notes

the `Dataset_regression` and `Dataset_permutation` functions in the
original implementation are essentially the same; the latter simply
permutes the trait values until the requisite number of significant
QTLs are found (or a limit is reached)

permutation only uses the `_2nRegression` function, too
