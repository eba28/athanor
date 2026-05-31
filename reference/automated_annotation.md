# Run automated cell type annotation

This function runs automated cell type annotation using `CellTypist`.

## Usage

``` r
automated_annotation(
  seurat_obj,
  annotation_method,
  reference = "pbmcref",
  majority_voting = FALSE,
  over_clustering = NULL
)
```

## Arguments

- seurat_obj:

  The Seurat object. Must be the path to a H5AD object if using
  CellTypist.

- annotation_method:

  Which method to use: CellTypist

- reference:

  Reference or model to use for prediction.

- majority_voting:

  Whether to enable majority voting for CellTypist predictions, which
  refines predictions based on local cluster information but may
  increase runtime.

## Value

A data.frame with the annotations for each cell

## Details

Supports CellTypist annotation methods. Assumes that the `Cells()` of
`seurat_obj` are properly formatted (i.e. unique). For `CellTypist`,
assumes the H5AD file and predictions have already been generated. This
would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md).

"Majority voting refines the prediction result in a local cell cluster
by choosing the dominant cell type label but may increase the runtime
especially for a large dataset due to the over-clustering step. This
approach usually improves the cell annotation, as voting is conducted in
small subclusters derived from over-clustering (cells belonging to a
given cell type will be assigned the same label regardless of potential
batch effects separating them)."

- https://www.celltypist.org/tutorials/onlineguide
