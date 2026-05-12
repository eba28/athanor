# Compute a permuted baseline for ADT neighbor matching scores

Runs ADT neighbor matching `n_permutations` times with shuffled
expression values and returns the mean score across permutations as a
random baseline. Supports both the mean absolute and range methods.

## Usage

``` r
permute_adt(
  seurat_obj,
  nn_names = names(seurat_obj@neighbors),
  features_adt,
  adt_assay = "ADT",
  adt_range = 0.1,
  methods = c("mean_abs", "range"),
  n_permutations = 10,
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

- adt_range:

  Relative threshold for the range method.

- methods:

  Character vector of methods. One or both of `"mean_abs"` and
  `"range"`.

- n_permutations:

  Number of permutations to run.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell values
  (averaged across permutations).

- verbose:

  If TRUE, print progress messages.

## Value

Data frame in the same format as
[`calc_adt_mean_absolute()`](https://eba28.github.io/athanor/reference/calc_adt_mean_absolute.md)
and
[`calc_adt_nn_within_range()`](https://eba28.github.io/athanor/reference/calc_adt_nn_within_range.md).
