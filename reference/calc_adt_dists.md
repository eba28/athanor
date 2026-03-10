# Compute mean ADT distance to each cell's k nearest neighbors

This function calculates the mean distance (or similarity) of each
cell's ADT profile to its k nearest neighbors in a specified neighbor
space (e.g., RNA, BCR, WNN) using a chosen distance metric (e.g., mean
absolute difference, Manhattan distance, Euclidean distance, or cosine
similarity).

## Usage

``` r
calc_adt_dists(
  seurat_obj,
  base_assay,
  adt_assay = "ADT",
  layer = "data",
  feature,
  k,
  multiple_k = TRUE,
  distance_metric = "mean_abs",
  return_mean = TRUE,
  exclude_self = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object that contains neighbor slots for different k's.

- base_assay:

  Which neighbor space to use (e.g., "RNA", "BCR", "w").

- adt_assay:

  ADT assay name (ADT, ADTnorm).

- layer:

  Data layer to pull from (data, counts, or scale.data).

- feature:

  Optional vector of ADT features; if missing, use all features present.

- k:

  Number of nearest neighbors.

- multiple_k:

  Whether to look for a neighbor slot specific to the provided k (e.g.,
  "RNA.nn_20") or just use the generic one (e.g., "RNA.nn"). The former
  allows you to have multiple neighbor graphs with different k's, while
  the latter assumes you only have one neighbor graph per assay.

- distance_metric:

  One of `c(mean_abs, manhattan, euclidean, cosine)`.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell.

- exclude_self:

  Drop the cell itself from neighbors if present.

## Value

A single numeric value if `return_mean = TRUE`, or a named numeric
vector of per-cell distances if `return_mean = FALSE`.

## Details

Pearson correlation coefficient measures linear relationships. Euclidean
distance measures the "straight-line" distance. Cosine similarity
measures the cosine of the angle between two vectors, with values closer
to 1 meaning that they are more similar.
