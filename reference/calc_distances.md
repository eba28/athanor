# This function calculates cluster distances for the given Seurat object.

This function calculates cluster distances for the given Seurat object.

## Usage

``` r
calc_distances(
  seurat_obj,
  reduction_name,
  criteria = "Within_Max",
  labels_true = "annotated_clusters",
  labels_name
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- reduction_name:

  The name of the reduction to use for distances.

- criteria:

  One of: Between_Mean, Between_Min, Calinski_Harabasz, Within_Between,
  Within_Max, Within_Mean, Within_Median.

- labels_true:

  The name of the column in the metadata that contains the cluster
  labels to evaluate.

- labels_name:

  A more descriptive name for the labels to use in plotting (optional).

## Value

A data.frame with a row per metric and cluster containing the score.

## Details

The [`as.factor()`](https://rdrr.io/r/base/factor.html) is needed in
case you give a categorical cluster col. Assumes you're using the
embeddings approach.
