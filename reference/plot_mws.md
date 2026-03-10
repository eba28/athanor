# Plot a box plot of modality weights per cell type

This function creates a box plot to visualize the distribution of
modality weights (e.g., RNA vs. BCR) across different cell types or
clusters in a post-WNN Seurat object.

## Usage

``` r
plot_mws(
  seurat_obj,
  details = "",
  second_assay = "BCR",
  clrs_specific = named_colors$mu_freq_bins,
  split_by = "mu_freq_bins",
  y_axis_label = "SHM Frequency Bins"
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- details:

  Details to add to the plot title.

- second_assay:

  List of other assays run through WNN in order.

- clrs_specific:

  A specific (must have names) color palette.

- split_by:

  A meta.data column to split the box plots up by.

- y_axis_label:

  Label for the y-axis.

## Value

A ggplot with the distribution of weights

## Details

Assumes annotated_clusters is a column
