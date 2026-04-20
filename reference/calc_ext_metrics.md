# Calculate external clustering metrics

This function calculates various external clustering metrics using the
`mclust` and `clevr` packages, which provide a comprehensive set of
clustering evaluation statistics based on true cluster labels and
predicted cluster assignments.

## Usage

``` r
calc_ext_metrics(
  seurat_obj,
  reduction_name,
  criteria = c("Completeness", "Homogeneity"),
  labels_true = "annotated_clusters",
  labels_pred = "seurat_clusters"
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- reduction_name:

  The name of the reduction to use for distances.

- criteria:

  One of: ARI, Completeness, Homogeneity

- labels_true:

  The name of the column in the metadata that contains the true cluster
  labels to evaluate.

- labels_pred:

  The name of the column in the metadata that contains the predicted
  cluster labels to evaluate.

## Value

A data.frame with a row per metric containing the combined score.

## Details

The Satija lab used `cluster` for their analyses. The
[`as.factor()`](https://rdrr.io/r/base/factor.html) is needed in case
you give a categorical cluster col. `sklearn.metrics.cluster` has
`completeness_score` and `homogeneity_score` Assumes you're using the
embeddings approach.
