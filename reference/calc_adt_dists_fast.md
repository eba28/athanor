# Compute per-cell mean ADT distance to its k nearest neighbors from a base assay

Compute per-cell mean ADT distance to its k nearest neighbors from a
base assay

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
