# Simulate gene expression data using Splatter and convert to Seurat object

Simulate gene expression data using Splatter and convert to Seurat
object

## Usage

``` r
sim_gex_splatter(
  num_genes = 1000,
  num_cells = 2000,
  splatter_groups = 1,
  splatter_method = "single",
  seed = 42,
  verbose = FALSE
)
```

## Arguments

- num_genes:

  Number of genes to simulate. Defaults to 1000.

- num_cells:

  Number of cells to simulate. Defaults to 2000.

- splatter_groups:

  Group probabilities for Splatter simulation. Defaults to 1.

- splatter_method:

  Splatter simulation method. Defaults to "single".

- seed:

  andom seed for reproducible simulations. Defaults to 42.

- verbose:

  Whether or not to print verbose output. Defaults to FALSE.

## Value

A Seurat object containing the simulated gene expression data with
metadata and cell identifiers

## Details

Splatter doesn't use a separator for its fake names.
