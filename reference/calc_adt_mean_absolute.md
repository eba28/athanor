# Calculate mean absolute ADT distance to nearest neighbors

For each cell, calculates the mean absolute difference between its ADT
expression and that of its nearest neighbors across specified features.

## Usage

``` r
calc_adt_mean_absolute(
  seurat_obj,
  nn_names = names(seurat_obj@neighbors),
  features_adt,
  adt_assay = "ADT",
  return_mean = TRUE,
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

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell
  values.

- verbose:

  If TRUE, print progress messages.

## Value

Data frame with columns: Graph, Assay, Feature, Method, Score. If
`return_mean = TRUE`, `Score` is the mean across all cells; else it
contains per-cell values with an additional `cell_id` column.

## Details

Use
[`permute_adt()`](https://eba28.github.io/athanor/reference/permute_adt.md)
to compute a permuted baseline. For the range-based score, see
[`calc_adt_nn_within_range()`](https://eba28.github.io/athanor/reference/calc_adt_nn_within_range.md).
