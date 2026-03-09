# This function plots an overview of a doublet identification method

This function plots an overview of a doublet identification method

## Usage

``` r
layout_doublets(
  seurat_obj,
  tissue_type,
  clrs_specific,
  use_hues = FALSE,
  group_col = "seurat_clusters",
  group_label = NULL,
  doublet_col = "scDblFinder.class",
  doublet_package = "scDblFinder",
  details = NULL
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- tissue_type:

  Blood, Skin.

- clrs_specific:

  The specific color palette (should be named).

- use_hues:

  Use the iwanthue hues instead of the default ggplot colors. Doesn't
  let you set any other settings.

- group_col:

  The column to group by.

- group_label:

  The label for the grouping variable to use in the plot titles and axis
  labels. If NULL, it will be determined based on the group_col name.

- doublet_col:

  The column containing the doublets information

- doublet_package:

  The doublet method being used.

- details:

  The optional subtitle.

## Value

A grid of four plots with UMAPs in the left column and bar plots in the
right column.

## Details

It assumes that named_colors\$doublet has been defined. Depends on other
plots. The doublets will be plotted "on top" for the first UMAP.
