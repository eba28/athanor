# Calculate the proportion of a cell's neighbors that are mutated or unmutated

This function calculates the proportion of neighbors that are mutated or
unmutated for each cell in a Seurat object. It uses the neighbor
information stored in the Seurat object to determine the mutation status
of neighboring cells and computes the fraction accordingly. The function
also generates a plot to visualize the distribution of mutated or
unmutated neighbors across different cell types.

## Usage

``` r
calc_nn_frac(seurat_obj, assay = "WNN", plot_title = "", category = "Mutated")
```

## Arguments

- seurat_obj:

  A Seurat object with neighbors calculated.

- assay:

  One of: RNA, BCR, WNN

- plot_title:

  A string to add to the plot title, e.g. "All Cells" or "Memory B Cells
  Only".

- category:

  One of: Mutated, Unmutated

## Value

A data.frame and a plot.

## Details

Assumes that WNN neighbors are stored in "w.nn". Defines mutated as
anything about 0% SHM.
[`setNames()`](https://rdrr.io/r/stats/setNames.html) works better than
a for loop checking all of the neighbors.
