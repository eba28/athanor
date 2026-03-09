# Reformat VDJ barcodes to make them unique across samples

This function processes barcode data by separating cell IDs from
suffixes and creating unique cell identifiers by combining sample names
with barcodes.

## Usage

``` r
reformat_vdj_barcode(
  data,
  col_samples = "sample_id",
  col_barcodes = "cell_id",
  col_output = "cell_id"
)
```

## Arguments

- data:

  Input data frame containing barcode information.

- col_samples:

  Column name containing sample names.

- col_barcodes:

  Column name containing cell IDs/barcodes.

- col_output:

  Name for the output column.

## Value

A character vector of unique cell identifiers.
