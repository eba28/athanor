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
  assay,
  reduc = "umap",
  plot_label = TRUE,
  label_box = TRUE,
  label_size = 3,
  annotated = FALSE,
  specific_clusters,
  meta_col = "seurat_clusters",
  include_legend = TRUE,
  legend_label,
  factor_idents = TRUE,
  order = FALSE,
  details,
  ...
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- data_source:

  The dataset of origin.

- clrs_specific:

  Specific colors for plotting (make sure it has names).

- use_hues:

  Use the iwanthue hues instead of the default ggplot colors. Doesn't
  let you set any other settings.

- pt_size:

  The point size.

- assay:

  The data type e.g. ADT, GEX, BCR, WNN... Will be used for the plot
  title.

- reduc:

  The reduction to use for plotting e.g. "bpca" or wnn.umap".

- plot_label:

  Add labels to the plot (or not).

- label_box:

  Whether or not to give the labels a background.

- label_size:

  The size of the plot labels.

- annotated:

  If the cell types have been identified.

- specific_clusters:

  Can overlay clusters of interest e.g. B cell or by \#. Overrides the
  annotated option.

- meta_col:

  Which column in the object metadata to color by. When combined with
  `specific_clusters`, highlights those values as an overlay instead of
  coloring all cells.

- include_legend:

  Include the legend or not.

- legend_label:

  The label for the legend.

- factor_idents:

  Whether or not to factorize the idents (for proper ordering of the
  colors). This can mess up the order you want, so be careful.

- order:

  Plot cells on top or not.

- details:

  A custom subtitle.

- ...:

  Any other Seurat parameters.

## Value

A Seurat plot of the specified reduction.
