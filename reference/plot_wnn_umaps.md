# This function gives a summary of the post-WNN object.

This function gives a summary of the post-WNN object.

## Usage

``` r
plot_wnn_umaps(
  seurat_obj,
  data_source = "Manual",
  airr_type = "BCR",
  airr_processing = "Embeddings",
  reducs_list = c("rna.umap", "bcr.umap", "wnn.umap"),
  clusters_col = "seurat_clusters",
  plot_label = TRUE,
  pt_size = 0.8,
  clrs_specific
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object.

- data_source:

  The source of the data for the plot title.

- airr_type:

  The type of AIRR data.

- airr_processing:

  The type of AIRR processing; one of `c("Embeddings", "Features)`

- reducs_list:

  Should be in order (GEX, AIRR, WNN).

- clusters_col:

  The metadata column to color the UMAPs by.

- plot_label:

  Whether or not to plot the metadata labels on the UMAPs.

- pt_size:

  The point size for the UMAPs.

- clrs_specific:

  A specific (must have names) color palette for the clusters. If not
  provided, the default Seurat colors will be used.

## Value

A combined plot of the GEX, BCR, and WNN UMAPs colored by the specified
metadata column.

## Details

Should be able to plot ADT instead of BCR too. Only plots one
clustering.
