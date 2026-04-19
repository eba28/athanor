# Calculate the correlation between each cell's expression and the mean of its neighbors' expression

Calculate the correlation between each cell's expression and the mean of
its neighbors' expression

## Usage

``` r
calc_correlation(
  seurat_obj,
  features_adt,
  adt_assay = "ADT",
  cor_method = "spearman"
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- features_adt:

  Name of the ADT features to evaluate (e.g. "CD27.1").

- adt_assay:

  Name of the assay containing ADT data.

- cor_method:

  Correlation method to use (e.g. "pearson", "spearman").

## Value

A data frame with columns: Graph, Feature, Score.

## Details

The Seurat object must have `FindNeighbors()` already run and an ADT
assay. Correlations are calculated across all neighbor graphs present in
the object.
