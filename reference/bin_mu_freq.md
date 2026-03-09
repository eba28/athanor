# Bin the mutation frequency

Bin the mutation frequency

## Usage

``` r
bin_mu_freq(seurat_obj, num_bins = c(2, 3, 5))
```

## Arguments

- seurat_obj:

  The Seurat object.

- num_bins:

  The number of bins to split `mu_freq` into. Must be at least one of 2,
  3, or 5.

## Value

The provided Seurat object with a new binned mu_freq column.

## Details

The bins are (most likely) not going to be equal sizes.
