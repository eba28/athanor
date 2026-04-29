# Plot a box plot of modality weights per cell type

This function creates a box plot to visualize the distribution of
modality weights (e.g. RNA vs. BCR) across different cell types or
clusters in a post-WNN Seurat object. It should be run after
[`run_wnn()`](https://eba28.github.io/athanor/reference/run_wnn.md).

## Usage

``` r
plot_mws(
  seurat_obj,
  second_assay = "BCR",
  clrs_specific = named_colors$mu_freq_bins,
  split_by = "mu_freq_bins",
  y_axis_label = "SHM Frequency Bins",
  details = ""
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- second_assay:

  List of other assays run through WNN in order.

- clrs_specific:

  The specific color palette (should be named).

- split_by:

  A meta.data column to split the box plots up by.

- y_axis_label:

  Label for the y-axis.

- details:

  An optional custom subtitle.

## Value

A ggplot with the distribution of weights

## Details

Assumes "annotated_clusters" is a column.
