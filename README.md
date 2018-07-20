# RNAMotif
RNA Motif generator and finder


## Dependencies

- the [SeqAn library](https://github.com/seqan/seqan) (tested with release 2.3.0) and its usual dependencies like zlib and OpenMP for parallel folding.  
- the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) headers and some compiled library like `libRNA.a`

The rest is self-contained in the repository.

## Known bugs

The ViennaRNA C library sometimes seems to crash when using multiple threads for folding. The cause is unclear since I don't think that I share state between threads, but maybe something leaks internally in the ViennaRNA library.