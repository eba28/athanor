# athanor

**Last Updated:** 2026-06-17

The goal of `athanor` is to integrate gene expression (GEX) and B cell
receptor (BCR) data. The package provides functions for data
preprocessing, feature extraction, concatenation, weighted nearest
neighbors, and plotting.

## Installation

You can install the development version of `athanor` from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("eba28/athanor")
```

*Please note that this package is still under active development.*

## Details

- [`automated_annotation()`](https://eba28.github.io/athanor/reference/automated_annotation.md)
  used to support Azimuth, but Azimuth is not currently available on
  CRAN and causes installation issues due to a hidden requirement of a
  Signac function that is no longer a part of Signac (`RunChromVAR()`)
- The functions were designed to operate on v5 of Seurat. If you are
  using an older version, you may have to change `layer` to `slot` in
  some of the code.
- Since the focus is on having multiple modalities, the GEX object
  created by
  [`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
  will have “rpca” and “rna.umap” as reductions instead of the usual
  “pca” and “umap” respectively. Similarly,
  [`bcr_embeddings_pipeline()`](https://eba28.github.io/athanor/reference/bcr_embeddings_pipeline.md)
  will return an object with “bpca” and “bcr.umap” reductions instead of
  “pca” and “umap”. This is to avoid confusion when the two objects are
  merged together in
  [`merge_gex_bcr()`](https://eba28.github.io/athanor/reference/merge_gex_bcr.md).
