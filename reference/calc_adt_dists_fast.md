# Compute mean ADT distance to each cell's k nearest neighbors (faster version)

This function calculates the mean distance of each cell's ADT profile to
its k nearest neighbors in a specified neighbor space (e.g., RNA, BCR,
WNN) using Euclidean distance. It returns a named numeric vector of
per-cell mean distances.

## Usage

``` r
calc_adt_dists_fast(adt_data, features, neighbors, exclude_self = TRUE)
```

## Arguments

- adt_data:

  ADT data matrix (features by cells) to use for distance calculations.

- features:

  Optional vector of ADT features; if missing, use all features present.

- neighbors:

  The kNN slot.

- exclude_self:

  Drop the cell itself from neighbors if present.

## Details

Only returns the mean.
