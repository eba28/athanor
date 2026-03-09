# Visualize how many cells are in each Seurat cluster or cell type.

Visualize how many cells are in each Seurat cluster or cell type.

## Usage

``` r
plot_counts_cluster(
  seurat_obj,
  tissue_type = "",
  clrs_specific,
  clusters_col,
  fill_col,
  fill_col_name,
  x_axis = "Cluster",
  details = NULL
)
```

## Arguments

- seurat_obj:

  The Seurat object containing the clusters and/or cell type annotations
  to plot.

- tissue_type:

  The tissue type of interest e.g. "Blood" or "Skin".

- clrs_specific:

  The specific color palette (should be named).

- clusters_col:

  The column to plot on the x axis (e.g. seurat_clusters).

- fill_col:

  The column to fill by (e.g. annotated_clusters).

- fill_col_name:

  The label for the fill aesthetic.

- x_axis:

  What to plot on the x axis: "Cluster" or "Cell Type".

- details:

  The optional subtitle.

## Value

A ggplot bar plot of cell counts.
