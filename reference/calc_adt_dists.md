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
  drop_self = TRUE
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

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) so `k` is fixed
  across cells regardless of whether self-matches are present. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

## Value

A data frame with columns: Graph, Assay, Feature, Method, Score (the
canonical schema shared by the `calc_*` evaluation functions). Distances
are computed jointly across all requested features, so `Feature` is the
"+"-joined set of features used (or "all" when `features_adt` is NULL)
rather than one row per feature. `Method` records the distance metric.
If `return_mean = TRUE`, `Score` is the mean across all cells; else it
contains per-cell values with an additional `cell_id` column.

## Details

Euclidean distance measures straight-line distance between vectors.
Manhattan distance sums absolute differences across features. Mean
absolute difference is the per-feature mean of absolute differences.
Cosine similarity measures the cosine of the angle between two vectors,
with values closer to 1 indicating greater similarity.
