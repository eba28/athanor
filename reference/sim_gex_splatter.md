# Simulate gene expression data using Splatter

This function simulates gene expression data using the Splatter package,
which generates synthetic single-cell RNA-seq data based on a specified
set of parameters. The function allows for customization of the number
of genes, cells, simulation method, group probabilities, and random seed
for reproducibility. The resulting simulated data is then converted into
a Seurat object with appropriate metadata and cell identifiers.

## Usage

``` r
sim_gex_splatter(
  num_genes = 1000,
  num_cells = 2000,
  splatter_groups = 1,
  splatter_method = "single",
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- num_genes:

  Number of genes to simulate.

- num_cells:

  Number of cells to simulate.

- splatter_groups:

  Group probabilities for Splatter simulation.

- splatter_method:

  Splatter simulation method.

- seed:

  Random seed for reproducible simulations.

- verbose:

  Logical indicating whether or not to print messages.

## Value

A Seurat object containing the simulated gene expression data with
metadata and cell identifiers

## Details

Splatter doesn't use a separator for its fake names.
