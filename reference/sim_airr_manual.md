# Simulate an assay made from AIRR embeddings independent of any GEX object

This function simulates an assay of AIRR embeddings by generating a
matrix of random values between -0.6 and 0.6, with a specified number of
cells and embedding dimensions.

## Usage

``` r
sim_airr_manual(num_cells, num_dims, separator = "-")
```

## Arguments

- num_cells:

  Number of cells to simulate.

- num_dims:

  Number of embedding dimensions to simulate.

- separator:

  Separator for cell and dimension names.

## Value

A Seurat Assay

## Details

The choice of value ranges was based off of real data, including having
9 decimal points.
