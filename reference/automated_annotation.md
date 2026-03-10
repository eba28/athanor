# Run automated cell type annotation

This function runs automated cell type annotation using either `Azimuth`
or `CellTypist`.

## Usage

``` r
automated_annotation(
  seurat_obj,
  annotation_method,
  reference = "pbmcref",
  azimuth_assay = "RNA",
  azimuth_levels = c("l1", "l2", "l3")
)
```

## Arguments

- seurat_obj:

  The Seurat object. Must be the path to a H5AD object if using
  CellTypist.

- annotation_method:

  Which method to use: c("Azimuth", "CellTypist")

- reference:

  Reference or model to use for prediction. Defaults to "pbmcref" (for
  Azimuth).

- azimuth_assay:

  Assay to use for Azimuth

- azimuth_levels:

  Levels to process for Azimuth e.g. c("l1", "l2", "l3")

## Value

A data.frame with the annotations for each cell

## Details

Supports Azimuth and CellTypist annotation methods. Assumes that the
`Cells()` of `seurat_obj` are properly formatted (i.e. unique). For
`Azimuth` with a Seurat v5 object, all of the layers have to be joined.
For `CellTypist`, assumes the H5AD file and predictions have already
been generated.
