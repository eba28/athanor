# Convert family and gene information to sorted factors

Converts the gene family and gene columns added by
[`add_family_info()`](https://eba28.github.io/athanor/reference/add_family_info.md)
to properly ordered factors using numeric sorting.

## Usage

``` r
factor_family_info(combined_airr)
```

## Arguments

- combined_airr:

  An AIRR-formatted data.frame. This function ensures that gene families
  and genes are ordered correctly (e.g., IGHV1, IGHV2, IGHV10 instead of
  IGHV1, IGHV10, IGHV2).

## Value

A data.frame with up to six columns converted to sorted factors.

## Details

For after
[`add_family_info()`](https://eba28.github.io/athanor/reference/add_family_info.md)
has been run
