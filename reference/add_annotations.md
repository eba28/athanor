# Add a user-specified list of cluster annotations to a Seurat object

This function adds a user-specified list of cluster annotations to the
Seurat object. It can be used for adding both manual and automated
annotations, updating cell identities and adding annotation columns to
the metadata.

## Usage

``` r
add_annotations(
  seurat_obj,
  annotations_df,
  cell_types_col = "CellType",
  relabel = TRUE,
  relocate = TRUE,
  alphabetize = TRUE,
  clusters_col = "seurat_clusters",
  annotations_col = "annotated_clusters"
)
```

## Arguments

- seurat_obj:

  The Seurat object to annotate.

- annotations_df:

  Data frame containing cluster-to-cell-type mappings, typically with
  "Cluster" and "CellType" columns.

- cell_types_col:

  The name of the column containing the cell type annotations.

- relabel:

  Whether to update the Seurat object's active identities. Sometimes you
  just want to add the metadata.

- relocate:

  Whether to relocate the annotation column in metadata.

- alphabetize:

  Whether to alphabetize the cell types.

- clusters_col:

  The name of the column containing cluster IDs.

- annotations_col:

  The name of the new metadata column for the annotations.

## Value

A Seurat object with added annotation information.

## Details

This would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md).
This assumes that cell typing was done on a cluster level, so the
annotations_df should have one row per cluster. If you have cell-level
annotations, you can skip the relabeling and just add the metadata
column to the Seurat object.
