# Plot a specific condition on a Seurat UMAPPlot.

Plot a specific condition on a Seurat UMAPPlot.

## Usage

``` r
plot_umap_condition(
  seurat_obj,
  tissue_type,
  clrs_specific,
  condition_name,
  operator,
  condition_val,
  color_by = "value",
  plot_type = "general",
  label_plot = TRUE,
  include_subtitle = TRUE,
  include_legend = FALSE
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- tissue_type:

  The tissue type e.g. Blood, Skin.

- clrs_specific:

  The specific color palette (should be named).

- condition_name:

  The column in the object that contains the condition of interest e.g.
  "annotated_clusters", "mu_freq", "isotype", etc.

- operator:

  \<, \>, ==

- condition_val:

  The value to compare the condition to. If color_by is "name", this
  should be a name in the condition_name column. If color_by is "value",
  this should be a value in the condition_name column.

- color_by:

  name, value

- plot_type:

  general, overlay (BCR/TCR)

- label_plot:

  Put labels on the plot (or not).

- include_subtitle:

  Include a subtitle (or not).

- include_legend:

  Include a legend (or not).

## Value

A Seurat UMAPPlot.

## Details

Based on: https://github.com/satijalab/seurat/issues/1053 Can be used
for plotting QC metrics, isolating specific cell types, overlaying AIRR
data, overlaying B cell isotypes, etc.
