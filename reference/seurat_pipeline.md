# Run Seurat's standard pipeline

This function takes a Seurat object containing combined GEX data and
runs the standard Seurat pipeline for normalization, scaling,
dimensionality reduction, clustering, and UMAP visualization. It
includes optional filtering of IG/TR genes and allows for customization
of various parameters such as the number of variable features, principal
components, and clustering resolution.

## Usage

``` r
seurat_pipeline(
  seurat_obj,
  nfeatures_RNA = 200,
  perc_mt = 15,
  num_features = 2000,
  num_pcs = 30,
  num_dims = 20,
  cluster_res = 0.4,
  filter_genes = TRUE,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  The Seurat object containing combined GEX data.

- nfeatures_RNA:

  Minimum number of RNA features.

- perc_mt:

  Maximum percentage of mitochondrial genes to retain.

- num_features:

  Desired number of variable features.

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of dimensions to use for neighbor finding and UMAP.

- cluster_res:

  Clustering resolution parameter.

- filter_genes:

  Whether to filter out IG/TR genes.

- verbose:

  Print out Seurat's progress messages.

## Value

A processed Seurat object with normalization, scaling, PCA, clustering,
and UMAP.

## Details

It is highly recommended to save the resulting object as an RDS or qs
file. For `filter_genes`, we assume that `features_meta` is already
loaded and `remove_genes` defined. This pipeline is loosely based on
[Seurat's
pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial). Unlike
previous analyses, `features = rownames(obj)` was removed from the
`ScaleData` step since the data is too large and only the top variable
features are needed to do `RunPCA`.
