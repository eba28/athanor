# Run Seurat's standard pipeline

Runs the standard Seurat pipeline (normalize, scale, PCA, neighbors,
UMAP) on a Seurat object. Optionally filters cells by QC metrics first,
filters IG/TR genes from variable features, and clusters. Supports any
assay (e.g. `"RNA"`, `"RNA_BCR"`), with reduction names derived
automatically from the assay name.

## Usage

``` r
seurat_pipeline(
  seurat_obj,
  assay = "RNA",
  pca_name = NULL,
  nfeatures_RNA,
  perc_mt,
  num_features = 2000,
  num_pcs = 50,
  num_dims = 20,
  k_param = 20,
  normalize = TRUE,
  find_var_features = TRUE,
  cluster_res = NULL,
  filter_genes,
  ensembl_version = NULL,
  cache_file = NULL,
  post_scale = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- assay:

  Name of the assay to run the pipeline on. Defaults to `"RNA"`.

- pca_name:

  Name to give the PCA reduction. If `NULL`, defaults to `"rpca"` for
  the RNA assay, or `paste0(tolower(assay), ".pca")` for others (e.g.
  `"rna_bcr.pca"` for `assay = "RNA_BCR"`).

- nfeatures_RNA:

  Minimum number of RNA features to retain per cell. If omitted, cell
  filtering is skipped. Only applies when `assay = "RNA"`.

- perc_mt:

  Maximum percentage of mitochondrial genes to retain. If omitted, cell
  filtering is skipped. Only applies when `assay = "RNA"`.

- num_features:

  Desired number of variable features.

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of PCA dimensions to use for neighbor finding.

- k_param:

  Number of nearest neighbors.

- normalize:

  If `TRUE`, normalize the assay using LogNormalize before scaling. Set
  to `FALSE` if the assay has already been normalized.

- find_var_features:

  If `TRUE`, run
  [`Seurat::FindVariableFeatures()`](https://satijalab.org/seurat/reference/FindVariableFeatures.html)
  before scaling. Set to `FALSE` if variable features are already set on
  the assay (e.g. when calling this after manually setting them in
  [`concatenate_gex_bcr()`](https://eba28.github.io/athanor/reference/concatenate_gex_bcr.md)).

- cluster_res:

  Clustering resolution(s). If `NULL`, clustering is skipped.

- filter_genes:

  If specified, filter out genes from this category (e.g. `"IG"` and/or
  `"TR"`).

- ensembl_version:

  Ensembl version for gene annotations (e.g. `"GRCh38.104"`). If `NULL`,
  uses the default in
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).

- cache_file:

  Passed to
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).
  Path to a cached RDS result to use instead of querying Ensembl.

- post_scale:

  Optional function applied to the `scale.data` matrix after
  [`Seurat::ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
  but before
  [`Seurat::RunPCA()`](https://satijalab.org/seurat/reference/RunPCA.html)
  (e.g. to reweight blocks of rows relative to one another before the
  joint PCA in
  [`concatenate_gex_bcr()`](https://eba28.github.io/athanor/reference/concatenate_gex_bcr.md)).
  Takes and returns a features-by-cells matrix.

- verbose:

  Logical indicating whether or not to print messages.

## Value

A processed Seurat object with PCA, neighbor graphs, optional clusters,
and UMAP. Reduction names are derived from `assay` and `pca_name`.

## Details

Cell filtering (`nfeatures_RNA`, `perc_mt`) and ADT normalization are
only applied when `assay = "RNA"`. For other assays, pass
`normalize = FALSE` if the data has already been normalized.
