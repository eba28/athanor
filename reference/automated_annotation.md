# Run automated cell type annotation

This function runs automated cell type annotation using `CellTypist`.

## Usage

``` r
automated_annotation(seurat_obj, annotation_method, reference = "pbmcref")
```

## Arguments

- seurat_obj:

  The Seurat object. Must be the path to a H5AD object if using
  CellTypist.

- annotation_method:

  Which method to use: CellTypist"

- reference:

  Reference or model to use for prediction.

## Value

A data.frame with the annotations for each cell

## Details

Supports CellTypist annotation methods. Assumes that the `Cells()` of
`seurat_obj` are properly formatted (i.e. unique). For `CellTypist`,
assumes the H5AD file and predictions have already been generated. This
would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md).
