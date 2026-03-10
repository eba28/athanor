# Calculate homogeneity scores for binary ADT features across embeddings

This function calculates homogeneity scores for binary ADT features
across different embeddings and reductions in a Seurat object.

## Usage

``` r
calc_adt_scores(
  seurat_objs,
  meta_res,
  metric_type,
  metrics,
  adt_features = "CD27.1",
  adt_cutoff = 1
)
```

## Arguments

- seurat_objs:

  List of Seurat objects for each embedding type.

- meta_res:

  Named list of cluster columns for each reduction.

- metric_type:

  One of: Internal, External

- metrics:

  List of metric names.

- adt_features:

  ADT feature name (e.g. "CD27.1"). You can provide multiple.

- adt_cutoff:

  Numeric cutoff for binary classification.

## Value

Data frame of scores for each embedding/reduction.
