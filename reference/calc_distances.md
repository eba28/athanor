# Calculate cluster distances for a Seurat object

This function calculates various cluster distance metrics using the
`fpc::cluster.stats` function, which provides a comprehensive set of
clustering statistics based on a distance matrix and cluster
assignments. The function computes both single-value metrics (e.g.,
Calinski-Harabasz index) and cluster-wise metrics (e.g., average
within-cluster distance) depending on the specified criteria. The
results are returned in a tidy data frame format for easy plotting and
comparison across different embeddings and reductions.

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
