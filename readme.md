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



`genotype.regression` outputs a list of QTLs, like:

```
QTL (Locus: "D18Mit4", Chr: 18, LRS: 0.717, Additive: 0.022)
QTL (Locus: "D19Mit109", Chr: 19, LRS: 5.504, Additive: -0.064)
QTL (Locus: "D19Mit44", Chr: 19, LRS: 2.342, Additive: -0.041)
QTL (Locus: "D19Mit42", Chr: 19, LRS: 2.342, Additive: -0.041)
QTL (Locus: "D19Mit22", Chr: 19, LRS: 5.169, Additive: -0.062)
QTL (Locus: "D19Mit127", Chr: 19, LRS: 4.641, Additive: -0.057)
```

in the example, it's 300 entries long, which is also the number
of markers present in BXD.txt


`genotype.permutation` outputs a list of numbers; in the example1.py case,
it's 1000 numbers long. it is sorted, asc, and the values range roughly
from 4.0 to 23.0
