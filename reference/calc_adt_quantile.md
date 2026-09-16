# Calculate the proportion of neighbors within an ADT marker's quantile by expression

For each cell in a Seurat object, calculates how many of its k nearest
neighbors have ADT expression within the same quantile as the cell's own
ADT expression for a given feature. Note that the output of this
function is a vector and not a data.frame like most of the other
evaluation functions.

## Usage

``` r
calc_adt_quantile(
  seurat_obj,
  adt_assay = "ADT",
  features_adt,
  base_assay,
  k = 20,
  n_quantile = 10,
  method = c("quantile", "percentile_diff"),
  binning = c("frequency", "width"),
  drop_self = TRUE,
  return_counts = FALSE,
  return_mean = TRUE,
  verbose = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing ADT data and computed neighbor graphs.

- adt_assay:

  Name of the assay containing ADT data.

- features_adt:

  Name(s) of the ADT feature(s) to evaluate (e.g. "CD27.1"); one row
  group per feature is returned.

- base_assay:

  The assay used to compute neighbors. One of "RNA", "GEX", "BCR", or
  "WNN".

- k:

  Number of nearest neighbors to evaluate. Must match the k used when
  computing the neighbor graph.

- n_quantile:

  The number of bins to divide ADT expression into. Neighbors count as
  "within range" when they fall in the same bin as the cell.

- method:

  One of `c("quantile", "percentile_diff")`. `"quantile"` compares
  discrete bins. `"percentile_diff"` returns the mean absolute
  percentile rank difference per cell.

- binning:

  How the `n_quantile` bins are defined for `method = "quantile"`.
  `"frequency"` (the default) makes bins holding roughly equal numbers
  of cells (true quantiles, via
  [`ggplot2::cut_number()`](https://ggplot2.tidyverse.org/reference/cut_interval.html));
  `"width"` makes bins of equal expression range (via
  [`ggplot2::cut_interval()`](https://ggplot2.tidyverse.org/reference/cut_interval.html)).
  Ignored for `method = "percentile_diff"`.

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) so `k` is fixed
  across cells regardless of whether self-matches are present. When
  TRUE, `k` is reduced by one to match. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

- return_counts:

  If TRUE, `Score` is the count of neighbors in the same bin (only for
  `method = "quantile"`); if FALSE, the proportion (count / k). Ignored
  for `method = "percentile_diff"`.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell
  values.

- verbose:

  If TRUE, print progress messages.

## Value

Data frame with columns: Graph, Assay, Feature, Method, Score (the
canonical schema shared by the `calc_*` evaluation functions). `Method`
is "Quantile" or "Percentile_Diff". If `return_mean = TRUE`, `Score` is
the mean across all cells; else it contains per-cell values with an
additional `cell_id` column.
