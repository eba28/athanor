# Concatenate GEX and BCR data in a Seurat object

Creates a new assay combining gene expression (GEX) and B-cell receptor
(BCR) data by concatenating processed BCR features with RNA data. Runs
the standard Seurat pipeline on the combined data.

## Usage

``` r
concatenate_gex_bcr(
  seurat_obj,
  pca_stage = "Before",
  cols_to_include,
  var_features = FALSE,
  normalize = TRUE,
  num_dims = 20,
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

- num_dims:

  number of dimensions to use for PCA and neighbor finding.

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
See `integrate_gex_airr()` for adding AIRR data to a Seurat object.
