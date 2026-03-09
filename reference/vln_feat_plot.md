# Plots a Seurat `VlnPlot` and a `FeaturePlot` side by side for the same marker

Plots a Seurat `VlnPlot` and a `FeaturePlot` side by side for the same
marker

## Usage

``` r
vln_feat_plot(
  seurat_obj,
  feature,
  assay = "RNA",
  group_col = NULL,
  rotate = FALSE
)
```

## Arguments

- seurat_obj:

  The Seurat object with GEX data.

- feature:

  The feature of interest.

- assay:

  The assay to search for the feature in.

- group_col:

  What to group by (uses the Idents by default).

- rotate:

  Rotate the labels or not.

## Value

Two patchworked Seurat plots.

## Details

Will put the highest expressing cells on top for the latter.
