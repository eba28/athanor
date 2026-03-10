# Manually simulate gene expression data

This function simulates gene expression data by generating a matrix of
counts using a Poisson distribution, where the rate parameter is set to
0.5. The resulting matrix is then converted to a sparse matrix format
and formatted with gene and cell names. Finally, a Seurat object is
created from the counts matrix, and metadata is added with cell
identifiers.

## Usage

``` r
sim_gex_manual(num_genes = 1000, num_cells = 2000, separator = "-")
```

## Arguments

- num_genes:

  Number of genes to simulate.

- num_cells:

  Number of cells to simulate.

- separator:

  Separator for gene and cell names.

## Value

A Seurat object

## Details

Partially based off of [Single Cell workshop: Chapter
6](http://gateway.training.ncgr.org/single-cell-workshop/seurat.md)
