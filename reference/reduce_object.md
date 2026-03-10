# Reduce a Seurat object's size

This function reduces a Seurat object by removing count matrices and
keeping only specified dimensionality reductions. This is especially
useful for creating lightweight objects for Shiny apps or sharing.

## Usage

``` r
reduce_object(
  seurat_obj,
  dim_reducs = "umap",
  print_size = TRUE,
  load_annotations = FALSE,
  annotations_file
)
```

## Arguments

- seurat_obj:

  Processed Seurat object to reduce.

- dim_reducs:

  Vector of dimensionality reductions to keep.

- print_size:

  Whether to print info about how much the object was reduced.

- load_annotations:

  Whether to load and add cell type annotations.

- annotations_file:

  File path to CSV file containing cluster and cell type annotations.

## Value

A reduced Seurat object with specified reductions kept.

## Details

Modify this as needed if your object is built differently (e.g. tSNE
instead). Change the annotation column name if needed. Uses DietSeurat
to remove counts while preserving reductions and metadata. This is
especially useful if you are making a Shiny app.
