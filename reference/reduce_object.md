# Reduce a Seurat object's size

This function reduces a Seurat object by removing count matrices and
keeping only specified dimensionality reductions. This is especially
useful for creating lightweight objects for Shiny apps or sharing.

## Usage

``` r
reduce_object(
  seurat_obj,
  dim_reducs = "umap",
  meta_cols,
  remove_neighbors = TRUE,
  print_size = TRUE,
  ...
)
```

## Arguments

- seurat_obj:

  Processed Seurat object to reduce.

- dim_reducs:

  Vector of dimensionality reductions to keep.

- meta_cols:

  Vector of metadata column names to keep. If unspecified, keeps all
  metadata columns.

- remove_neighbors:

  Whether or not to remove neighbor graphs from the object.

- print_size:

  Whether to print info about how much the object was reduced.

## Value

A reduced Seurat object with specified reductions kept.

## Details

Modify this as needed if your object is built differently (e.g. tSNE
instead). Uses `DietSeurat` to remove counts while preserving reductions
and metadata. This is especially useful if you are making a Shiny app or
just visualizing the data.
