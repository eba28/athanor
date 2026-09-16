# Calculate the correlation between each cell's expression and the mean of its neighbors' expression

Calculate the correlation between each cell's expression and the mean of
its neighbors' expression

## Usage

``` r
calc_adt_correlation(
  seurat_obj,
  nn_names = names(seurat_obj@neighbors),
  features_adt,
  adt_assay = "ADT",
  cor_method = "spearman",
  drop_self = TRUE,
  verbose = FALSE
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- nn_names:

  Character vector of neighbor graph names to evaluate. Defaults to all
  graphs in the object.

- features_adt:

  Name of the ADT features to evaluate (e.g. "CD27.1").

- adt_assay:

  Name of the assay containing ADT data.

- cor_method:

  Correlation method to use (e.g. "pearson", "spearman").

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) so `k` is fixed
  across cells regardless of whether self-matches are present. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

- verbose:

  Logical indicating whether or not to print messages.

## Value

A data frame with columns: Graph, Assay, Feature, Method, Score (the
canonical schema shared by the `calc_*` evaluation functions). `Method`
records the correlation method (e.g. "Spearman").

## Details

The Seurat object must have `FindNeighbors()` already run and an ADT
assay. Correlations are calculated across all neighbor graphs present in
the object.
