# Plot an overview of a doublet identification method

This function creates a grid of four plots to visualize the results of a
doublet identification method. The left column contains UMAP plots
colored by doublet/singlet status and by clusters/annotations, while the
right column contains bar plots showing the counts and percentages of
doublets across clusters or annotations.

## Usage

``` r
plot_doublets(
  seurat_obj,
  data_source,
  clrs_specific,
  use_hues = FALSE,
  reduc = "rna.umap",
  meta_col = "seurat_clusters",
  group_label = NULL,
  doublet_col = "scDblFinder.class",
  doublet_package = "scDblFinder",
  details = NULL
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- data_source:

  Dataset description.

- clrs_specific:

  The specific color palette (should be named).

- use_hues:

  Use the `iwanthue` hues instead of the default ggplot colors. Doesn't
  let you set any other settings.

- reduc:

  The reduction to use for plotting e.g. "bpca" or wnn.umap".

- meta_col:

  The column to group by.

- group_label:

  The label for the grouping variable to use in the plot titles and axis
  labels. If NULL, it will be determined based on the meta_col name.

- doublet_col:

  The column containing the doublets information

- doublet_package:

  The doublet method being used.

- details:

  An optional custom subtitle.

## Value

A grid of four plots with UMAPs in the left column and bar plots in the
right column.

## Details

It assumes that `named_colors$doublet` has been defined. Depends on
other plots. The doublets will be plotted "on top" for the first UMAP.
