# Plot a box plot of modality weights per cell type

This function creates a box plot to visualize the distribution of
modality weights (e.g. RNA vs. BCR) across different cell types or
clusters in a post-WNN Seurat object. It should be run after
[`run_wnn()`](https://eba28.github.io/athanor/reference/run_wnn.md).

## Usage

``` r
plot_mws(
  seurat_obj,
  data_source = "",
  second_assay = "BCR",
  clrs_specific = named_colors$mu_freq_bins,
  split_by = "mu_freq_bins",
  facet_by = "annotated_clusters_simpler",
  y_axis_label = "SHM Frequency Bins",
  details = NULL
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- data_source:

  Dataset description.

- second_assay:

  The second assay to compare against RNA (e.g. "BCR"). If multiple
  assays are provided, it will assume the weight column is in the format
  "main_assay.weight" (e.g. "BCR.weight").

- clrs_specific:

  The specific color palette (should be named).

- split_by:

  A meta.data column to split the box plots up by.

- facet_by:

  A meta.data column to facet the box plots by.

- y_axis_label:

  Label for the y-axis.

- details:

  An optional custom subtitle.

## Value

A ggplot with the distribution of weights

## Details

Assumes "annotated_clusters" is a column.
