# athanor

**Last Updated:** 2026-04-19

The goal of `athanor` is to integrate gene expression (GEX) and B cell
receptor (BCR) data to improve B cell subtyping. The package provides
functions for data preprocessing, feature extraction, concatenation,
weighted nearest neighbors, and plotting.

## Installation

You can install the development version of `athanor` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("eba28/athanor")
```

## Example

``` r
library(athanor)

## basic example code
```

## Details

- [`automated_annotation()`](https://eba28.github.io/athanor/reference/automated_annotation.md)
  used to support Azimuth, but Azimuth is not currently available on
  CRAN and causes installation issues due to a hidden requirement of a
  Signac function that is no longer a part of Signac (`RunChromVAR()`)
- The functions were designed to operate on v5 of Seurat. If you are
  using an older version, you may have to change `layer` to `slot` in
  some of the code.
