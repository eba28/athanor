# Map cell types to Seurat clusters

This function maps cell types to Seurat clusters by counting the number
of cells in each cluster that belong to each cell type, and then
assigning the most common cell type to each cluster. It returns a data
frame with the assigned cell types for each cluster.

## Usage

``` r
cell_type_clusters(
  seurat_obj,
  clusters_col = "seurat_clusters",
  annotations_col
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- clusters_col:

  The metadata column with the Seurat clusters.

- annotations_col:

  The metadata column with the cell types.

## Value

A data.frame with a row for each cell type.
