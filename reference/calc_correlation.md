# Calculate the correlation between each cell's expression and the mean of its neighbors' expression

Calculate the correlation between each cell's expression and the mean of
its neighbors' expression

## Usage

``` r
calc_correlation(seurat_obj, features_adt, cor_method = "spearman")
```

## Arguments

- seurat_obj:

  The Seurat object.

- features_adt:

  Name of the ADT features to evaluate on (e.g. "CD27.1").

- cor_method:

  Correlation method to use (e.g. "pearson", "spearman").

## Value

A data frame with columns: Graph, Feature, Score (correlation value).

## Details

The Seurat object must have `FindNeighbors()` already run at least one
time and an assay named "ADT".
