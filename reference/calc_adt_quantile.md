# Calculate the proportion of neighbors within an ADT marker's quantile by expression

For each cell in a Seurat object, this function calculates how many of
its k nearest neighbors have ADT expression within the same quantile as
the cell's own ADT expression for a given feature.

## Usage

``` r
calc_adt_quantile(
  seurat_obj,
  adt_assay = "ADT",
  feature,
  base_assay,
  k = 20,
  use_k = TRUE,
  n_quantile = 10,
  return_counts = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing ADT data and computed neighbor graphs.

- adt_assay:

  Name of the assay containing ADT data.

- feature:

  Name of the ADT feature to evaluate (e.g. "CD27.1", "CD38").

- base_assay:

  The assay used to compute neighbors. One of "RNA", "GEX", "BCR", or
  "WNN".

- k:

  Numeric. Number of nearest neighbors to evaluate. Must match the k
  used when computing the neighbor graph.

- use_k:

  Logical. Whether to look for a neighbor slot specific to the provided
  k (e.g. "RNA.nn_20") or just use the generic one (e.g. "RNA.nn"). The
  former allows you to have multiple neighbor graphs with different k's,
  while the latter assumes you only have one neighbor graph per assay.

- n_quantile:

  Numeric. The number of quantiles to divide the ADT expression into.
  Neighbors are considered "within range" if they fall into the same
  quantile as the cell.

- return_counts:

  Logical. If TRUE, returns the count of neighbors within range. If
  FALSE, returns the proportion (count/k).

## Value

A named numeric vector with one value per cell in the Seurat object. If
`return_counts = TRUE`, returns the count of neighbors within the same
quantile If `return_counts = FALSE`, returns the proportion of neighbors
within the same quantile (ranging from 0 to 1). Vector names are cell
ids.
