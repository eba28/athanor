# calculate the proportion of neighbors within an ADT expression range

For each cell in a Seurat object, this function calculates how many of
its k nearest neighbors have ADT expression within a specified threshold
(default 20%) of the cell's own ADT expression for a given feature.

## Usage

``` r
calc_adt_nn_within_range(
  seurat_obj,
  adt_assay = "ADT",
  feature,
  base_assay,
  k = 20,
  use_k = TRUE,
  range = 0.2,
  return_counts = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing ADT data and computed neighbor graphs.

- adt_assay:

  Character. Name of the assay containing ADT data.

- feature:

  Character. Name of the ADT feature to evaluate (e.g., "CD27.1",
  "CD38.1").

- base_assay:

  Character. The assay used to compute neighbors. One of "RNA", "GEX",
  "BCR", or "WNN".

- k:

  Numeric. Number of nearest neighbors to evaluate. Must match the k
  used when computing the neighbor graph.

- use_k:

  Logical. Whether to look for a neighbor slot specific to the provided
  k (e.g., "RNA.nn_20") or just use the generic one (e.g., "RNA.nn").
  The former allows you to have multiple neighbor graphs with different
  k's, while the latter assumes you only have one neighbor graph per
  assay.

- range:

  Numeric. The relative threshold for considering neighbors similar. A
  value of 0.20 means neighbors within ±20% of the cell's expression are
  counted.

- return_counts:

  Logical. If TRUE, returns the count of neighbors within range. If
  FALSE, returns the proportion (count/k).

## Value

A named numeric vector with one value per cell in the Seurat object. If
`return_counts = TRUE`, returns the count of neighbors within range. If
`return_counts = FALSE`, returns the proportion of neighbors within
range (ranging from 0 to 1). Vector names are cell ids.

## Details

The range is symmetric around the cell's expression value. For example,
with range = 0.20:

- Lower bound = cell_expr \* 0.80

- Upper bound = cell_expr \* 1.20
