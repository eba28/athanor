# Reformat VDJ barcodes by adding sample names and removing suffixes

This function adds the sample name to the beginning of barcodes and
removes the "-1" suffixes from the ends to create unique cell
identifiers.

## Usage

``` r
reformat_vdj_barcode_sample(data, sample_name, barcode_col = "barcode")
```

## Arguments

- data:

  The input dataset which contains a barcode column.

- sample_name:

  The sample name to add to make the barcode unique.

- barcode_col:

  The column in the data that contains the barcodes.

## Value

A character vector of reformatted barcodes.
