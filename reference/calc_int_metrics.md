# This function calculates internal clustering metrics for the given Seurat object.

This function calculates internal clustering metrics for the given
Seurat object.

## Usage

``` r
calc_int_metrics(
  seurat_obj,
  reduction_name,
  criteria = "Silhouette",
  labels_true = "annotated_clusters",
  labels_name,
  return_full = FALSE
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- reduction_name:

  The name of the reduction to use for distances.

- criteria:

  One of: DB, Dunn, Intra_Complete, Silhouette

- labels_true:

  The name of the column in the metadata that contains the cluster
  labels to evaluate.

- labels_name:

  A more descriptive name for the labels to use in plotting (optional).

- return_full:

  If TRUE, return the full silhouette object instead of just the mean
  silhouette width.

## Value

A data.frame with a row per metric containing the combined score.

## Details

The Satija lab used `cluster` for their analyses. The
[`as.factor()`](https://rdrr.io/r/base/factor.html) is needed in case
you give a categorical cluster col. Use the `bluster` package for speed
if desired e.g. bluster::approxSilhouette(x = embeddings, clusters =
clusters) Assumes you're using the embeddings approach.
