# This function plots a Seurat UMAP(s) in several useful ways.

This function plots a Seurat UMAP(s) in several useful ways.

## Usage

``` r
plot_umap(
  seurat_obj,
  tissue_type = "",
  clrs_specific,
  use_hues = FALSE,
  plot_by = "all",
  specific_clusters = c(),
  specific_col,
  plot_label = TRUE,
  label_size = 3,
  label_box = TRUE,
  clusters_col = "seurat_clusters",
  annotated = FALSE,
  annotations_col = "annotated_clusters",
  annotations_type,
  order = FALSE,
  include_legend = TRUE,
  ncol = 4,
  details
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

- plot_by:

  What to plot by: by dataset, by sample type (control vs EM/blood), by
  cluster, by sample, with all samples, or with all subjects.

- specific_clusters:

  Can overlay clusters of interest e.g. B cell or by \#.

- specific_col:

  Overlay a specific column in the object e.g. "annotated_clusters" or
  "sample_id". Overrides the other options.

- plot_label:

  Add labels to the plot (or not).

- label_size:

  The size of the plot labels.

- label_box:

  Whether or not to give the labels a background.

- clusters_col:

  Which column in the object stores the clusters.

- annotated:

  If the cell types have been identified.

- annotations_col:

  Which column in the object stores the cell types.

- annotations_type:

  The method for annotation e.g. Manual, singleR, etc..

- order:

  Plot cells on top or not.

- include_legend:

  Include the legend or not.

- ncol:

  The number of columns if outputting multiple plots.

- details:

  A custom subtitle.

## Value

A Seurat UMAPPlot.

## Details

Includes whether or not the object has been annotated with specific cell
types.
