# Calculate the proportion of neighbors within an ADT marker's expression range

For each cell in a Seurat object, calculates how many of its k nearest
neighbors have ADT expression within a specified threshold (default 20%)
of the cell's own ADT expression for a given feature.

## Usage

``` r
calc_adt_nn_within_range(
  seurat_obj,
  adt_assay = "ADT",
  features_adt,
  base_assay,
  k = 20,
  range = 0.2,
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

- range:

  The relative threshold for considering neighbors similar. A value of
  0.20 means neighbors within +/-20% of the cell's expression are
  counted.

- return_counts:

  If TRUE, returns the count of neighbors within range. If FALSE,
  returns the proportion (count/k).

## Value

A named numeric vector with one value per cell in the Seurat object.

## Details

The range is symmetric around the cell's expression value. For example,
with range = 0.20:

- Lower bound = cell_expr \* 0.80

- Upper bound = cell_expr \* 1.20
