# Compute mean ADT distance to each cell's k nearest neighbors

Calculates the mean distance of each cell's ADT profile to its k nearest
neighbors in a specified neighbor space (e.g. RNA, BCR, WNN) using a
chosen distance metric (mean absolute difference, Manhattan, Euclidean,
or cosine similarity).

## Usage

``` r
calc_adt_dists(
  seurat_obj,
  base_assay,
  adt_assay = "ADT",
  layer = "data",
  features_adt = NULL,
  k,
  distance_metric = "mean_abs",
  return_mean = TRUE,
  exclude_self = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object that contains neighbor slots.

- base_assay:

  Which neighbor space to use (e.g. "RNA", "BCR", "WNN").

- adt_assay:

  ADT assay name (e.g. "ADT", "ADTnorm").

- layer:

  Data layer to pull from (data, counts, or scale.data).

- features_adt:

  Optional vector of ADT features; if NULL, use all features present.

- k:

  Number of nearest neighbors.

- distance_metric:

  One of `c("mean_abs", "manhattan", "euclidean", "cosine")`.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell
  values.

- exclude_self:

  Drop the cell itself from its neighbor list if present.

## Value

A single numeric value if `return_mean = TRUE`, or a named numeric
vector of per-cell distances if `return_mean = FALSE`.

## Details

Euclidean distance measures straight-line distance between vectors.
Manhattan distance sums absolute differences across features. Mean
absolute difference is the per-feature mean of absolute differences.
Cosine similarity measures the cosine of the angle between two vectors,
with values closer to 1 indicating greater similarity.
