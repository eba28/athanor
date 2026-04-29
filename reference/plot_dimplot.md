# Plot a Seurat UMAP using `DimPlot`

This function generates a UMAP plot from a Seurat object using `DimPlot`
with various customizable options for coloring, labeling, and grouping
the data.

## Usage

``` r
plot_dimplot(
  seurat_obj,
  data_source = "",
  clrs_specific,
  use_hues = FALSE,
  pt_size = 0.2,
  title,
  reduc = "rna.umap",
  meta_col = "annotated_clusters",
  highlight,
  plot_label = TRUE,
  label_size = 3,
  label_box = TRUE,
  include_legend = TRUE,
  legend_label,
  sort_idents = TRUE,
  idents_char = TRUE,
  order = FALSE,
  details,
  ...
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

  Use the iwanthue hues instead of the default ggplot colors. Doesn't
  let you set any other settings.

- pt_size:

  The point size.

- title:

  The plot title.

- reduc:

  The reduction to use for plotting e.g. "bpca" or wnn.umap".

- meta_col:

  Which column in the object metadata to color by. When combined with
  `highlight`, highlights those values as an overlay instead of coloring
  all cells.

- highlight:

  Can overlay clusters of interest e.g. B cell or by \#. Overrides the
  annotated option.

- plot_label:

  Add labels to the plot (or not).

- label_size:

  The size of the plot labels.

- label_box:

  Whether or not to give the labels a background.

- include_legend:

  Include the legend or not.

- legend_label:

  The label for the legend.

- sort_idents:

  Whether or not to sort the idents (for proper ordering of the colors).
  This can mess up the order you want, so be careful.

- idents_char:

  If sorting idents, whether to sort them as characters or numerically
  (e.g. cluster 10 should be after cluster 9, not before).

- order:

  Plot cells on top or not.

- details:

  An optional custom subtitle.

- ...:

  Any other Seurat parameters.

## Value

A Seurat plot of the specified reduction.
