# This function calculation the proportion of the neighbors that are mutated or unmutated.

This function calculation the proportion of the neighbors that are
mutated or unmutated.

## Usage

``` r
calc_nn_frac(seurat_obj, assay = "WNN", plot_title = "", category = "mutated")
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
