# Calculate the proportion of neighbors within an ADT marker's quantile by expression

For each cell in a Seurat object, calculates how many of its k nearest
neighbors have ADT expression within the same quantile as the cell's own
ADT expression for a given feature.

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
  return_counts = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing ADT data and computed neighbor graphs.

- adt_assay:

  Name of the assay containing ADT data.

- features_adt:

  Name of the ADT feature to evaluate (e.g. "CD27.1").

- base_assay:

  The assay used to compute neighbors. One of "RNA", "GEX", "BCR", or
  "WNN".

- k:

  Number of nearest neighbors to evaluate. Must match the k used when
  computing the neighbor graph.

- n_quantile:

  The number of quantiles to divide ADT expression into. Neighbors are
  considered "within range" if they fall into the same quantile as the
  cell.

- method:

  One of `c("quantile", "percentile_diff")`. `"quantile"` compares
  discrete quantile bins. `"percentile_diff"` returns mean absolute
  percentile rank difference per cell.

- return_counts:

  If TRUE, returns the count of neighbors within the same quantile bin
  (only applies to `method = "quantile"`). If FALSE, returns the
  proportion (count/k). Ignored when `method = "percentile_diff"`.

## Value

A named numeric vector with one value per cell in the Seurat object. For
`method = "quantile"`: count or proportion of neighbors in the same
quantile bin. For `method = "percentile_diff"`: mean absolute percentile
rank difference per cell (0 to 1).
