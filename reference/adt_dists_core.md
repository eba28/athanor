# Compute per-cell mean distances to nearest neighbors

Compute per-cell mean distances to nearest neighbors

## Usage

``` r
adt_dists_core(
  adt_mat,
  nn_idx,
  distance_metric = "euclidean",
  exclude_self = TRUE,
  return_mean = FALSE
)
```

## Arguments

- adt_mat:

  A cells by features numeric matrix.

- nn_idx:

  An integer matrix of nearest-neighbor indices (cells by k).

- distance_metric:

  One of `"euclidean"`, `"manhattan"`, `"mean_abs"`, or `"cosine"`.

- exclude_self:

  Logical; whether to drop a cell from its own neighbor list before
  computing distances.

- return_mean:

  Logical; if `TRUE` returns a single grand mean instead of a per-cell
  vector.

## Value

A named numeric vector of per-cell mean distances, or a scalar if
`return_mean = TRUE`.
