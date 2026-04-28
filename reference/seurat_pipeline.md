# Run Seurat's standard pipeline

Runs the standard Seurat pipeline (normalize, scale, PCA, neighbors,
UMAP) on a Seurat object. Optionally filters cells by QC metrics first,
filters IG/TR genes from variable features, and clusters.

## Usage

``` r
seurat_pipeline(
  seurat_obj,
  nfeatures_RNA,
  perc_mt,
  num_features = 2000,
  num_pcs = 50,
  num_dims = 30,
  k_param = 20,
  cluster_res = NULL,
  filter_genes,
  ensembl_version = NULL,
  cache_file = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  The Seurat object containing combined GEX data.

- nfeatures_RNA:

  Minimum number of RNA features. If omitted, cell filtering is skipped.

- perc_mt:

  Maximum percentage of mitochondrial genes to retain. If omitted, cell
  filtering is skipped.

- num_features:

  Desired number of variable features.

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of PCA dimensions to use for neighbor finding.

- k_param:

  Number of nearest neighbors.

- cluster_res:

  Clustering resolution(s). If NULL, clustering is skipped.

- filter_genes:

  If specified, filter out genes from this category (e.g. "IG" and/or
  "TR").

- ensembl_version:

  Ensembl version for gene annotations (e.g. "GRCh38.104"). If NULL,
  uses the default in
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).

- cache_file:

  Passed to
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).
  Path to a cached RDS result to use instead of querying Ensembl.

- verbose:

  Print out Seurat's progress messages.

## Value

A processed Seurat object with PCA (`rpca`), neighbor graphs (`RNA_nn`,
`RNA_snn`, `RNA.nn`), optional clusters, and UMAP (`rna.umap`).

## Details

It is highly recommended to save the resulting object as an RDS or qs
file. This pipeline is loosely based on [Seurat's
pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial). Unlike
previous analyses, `features = rownames(obj)` was removed from the
`ScaleData` step since the data is too large and only the top variable
features are needed to do `RunPCA`.
