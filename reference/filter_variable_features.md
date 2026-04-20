# Filter AIRR genes from variable features

Removes IG and/or TR genes from a Seurat object's variable features
list. Optionally reports how many remaining features are GEX-only when
BCR features are present.

## Usage

``` r
filter_variable_features(
  seurat_obj,
  filter_genes,
  ensembl_version = NULL,
  bcr_features = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- filter_genes:

  Category of genes to remove (e.g. "IG" and/or "TR").

- ensembl_version:

  Ensembl version for gene annotations (e.g. "v114"). If NULL, uses the
  default in
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).

- bcr_features:

  Optional matrix of BCR features (rows = features). If provided, the
  log message also reports the number of GEX-only features.

## Value

The Seurat object with filtered variable features and Ensembl version
saved to `Misc(seurat_obj, "ensembl_version")`.

## Details

This will usually just be used as part of
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md).
