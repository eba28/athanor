# Calculate the proportion of neighbors within an ADT marker's expression range

For each cell, calculates the proportion of its nearest neighbors whose
ADT expression falls within a symmetric relative threshold of the cell's
own expression. For example, with `range = 0.20`, neighbors within
+/-20% of the cell's expression are counted.

## Usage

``` r
calc_adt_nn_within_range(
  seurat_obj,
  nn_names = names(seurat_obj@neighbors),
  features_adt,
  adt_assay = "ADT",
  range = 0.2,
  return_mean = TRUE,
  drop_self = TRUE,
  verbose = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- nn_names:

  Character vector of neighbor graph names to evaluate. Defaults to all
  graphs in the object.

- features_adt:

  Character vector of ADT feature names to evaluate.

- adt_assay:

  Name of the ADT assay.

- range:

  Relative threshold. A value of 0.20 means neighbors within +/-20% of
  the cell's expression are counted as matches.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell
  values.

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) so `k` is fixed
  across cells regardless of whether self-matches are present. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

- verbose:

  If TRUE, print progress messages.

## Value

Data frame with columns: Graph, Assay, Feature, Method, Score. If
`return_mean = TRUE`, `Score` is the mean across all cells; else it
contains per-cell values with an additional `cell_id` column.
