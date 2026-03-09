# This function maps the cell types to the Seurat clusters

This function maps the cell types to the Seurat clusters

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
