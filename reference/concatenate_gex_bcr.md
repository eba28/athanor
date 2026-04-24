# Concatenate GEX and BCR data in a Seurat object

Creates a new assay combining gene expression (GEX) and B-cell receptor
(BCR) data by concatenating processed BCR features with RNA data. Runs
the standard Seurat pipeline on the combined data.

## Usage

``` r
concatenate_gex_bcr(
  seurat_obj,
  pca_stage = c("Before", "After"),
  cols_to_include,
  var_features = FALSE,
  normalize = TRUE,
  num_pcs = 50,
  num_dims = 20,
  k_param = 20,
  filter_genes,
  ensembl_version = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing RNA assay and BCR metadata.

- pca_stage:

  Add BCR information before PCA or after PCA.

- cols_to_include:

  Character vector of BCR metadata column names to include in the
  concatenated assay e.g. c("mu_freq", "isotype") or embedded
  dimensions.

- var_features:

  If TRUE, run FindVariableFeatures on the combined assay. If FALSE,
  concatenate BCR features onto existing variable features.

- normalize:

  If TRUE, normalize the combined assay using LogNormalize. If FALSE,
  skip normalization.

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of PCA dimensions to use for neighbor finding and UMAP.

- k_param:

  Number of nearest neighbors.

- filter_genes:

  If specified, filter out genes from this category (e.g. "IG" and/or
  "TR")

- ensembl_version:

  If filtering genes, specify the Ensembl version to use for gene
  annotations (e.g. "GRCh38.104"). If NULL, uses the default version in
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).

## Value

A Seurat object with:

- New `RNA_BCR` assay containing concatenated GEX and BCR features

- PCA reduction (`rna_bcr.pca`)

- UMAP reduction (`rna_bcr.umap`)

- Neighbor graphs computed on the combined data

## Details

This would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
and
[`gex_add_airr()`](https://eba28.github.io/athanor/reference/gex_add_airr.md).
The function:

1.  Extracts specified BCR metadata columns from the Seurat object

2.  Processes BCR features using
    [`process_bcr_features()`](https://eba28.github.io/athanor/reference/process_bcr_features.md)

3.  Creates a new assay by row-binding RNA and BCR data

4.  Optionally filters out IG/TR genes from variable features

5.  Runs standard Seurat workflow: normalize, scale, PCA, neighbors,
    UMAP

## Note

Currently assumes BCR data is already integrated into object metadata.
