# Plot a Seurat `VlnPlot` and a `FeaturePlot` side by side for the same marker

This function generates a side-by-side visualization of a Seurat
`VlnPlot` and a `FeaturePlot` for a specified marker, allowing for a
comprehensive comparison of expression levels across different groups
and spatial distribution on the UMAP. The `VlnPlot` will display the
distribution of expression levels across specified groups, while the
`FeaturePlot` will show the spatial localization of the marker on the
UMAP, with higher expressing cells highlighted on top for better
visibility.

## Usage

``` r
plot_vln_feat(
  seurat_obj,
  clrs_specific,
  feature,
  reduc = "umap",
  meta_col = NULL,
  rotate = FALSE,
  ...
)
```

## Arguments

- seurat_obj:

  The Seurat object with GEX data.

- clrs_specific:

  The specific color palette (should be named).

- feature:

  The feature of interest.

- reduc:

  The reduction to use for plotting e.g. "bpca" or wnn.umap".

- meta_col:

  What to group by (uses the Idents by default).

- rotate:

  Rotate the labels or not.

## Value

Two patchworked Seurat plots.

## Details

Will put the highest expressing cells on top for the latter.
