# This function plots a Seurat UMAP using `DimPlot`

This function plots a Seurat UMAP using `DimPlot`

## Usage

``` r
plot_dimplot(
  seurat_obj,
  data_source = "",
  clrs_specific,
  use_hues = FALSE,
  pt_size = 0.2,
  assay,
  reduc,
  plot_label = TRUE,
  label_box = TRUE,
  label_size = 3,
  annotated = FALSE,
  specific_clusters,
  clusters_col = "seurat_clusters",
  annotations_col = "annotated_clusters",
  include_legend = TRUE,
  legend_label,
  factor_idents = TRUE,
  details
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

  The data type e.g. ADT, GEX, BCR, WNN...

- reduc:

  The reduction to use for plotting e.g. wnn.umap

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

- clusters_col:

  Which column in the object stores the clusters.

- annotations_col:

  Which column in the object stores the cell types.

- include_legend:

  Include the legend or not.

- legend_label:

  The label for the legend.

- factor_idents:

  Whether or not to factorize the idents (for proper ordering of the
  colors). This can mess up the order you want, so be careful.

- details:

  A custom subtitle.

## Value

A Seurat UMAPPlot.
