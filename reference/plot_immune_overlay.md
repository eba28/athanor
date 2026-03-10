# Plot UMAPs with AIRR overlays

This function plots the annotated UMAPs alongside the BCR/TCR overlays.

## Usage

``` r
plot_immune_overlay(
  seurat_obj,
  tissue_type,
  airr_type,
  clrs_specific,
  plot_by = "all",
  barcode_col = "cell_id",
  plot_label = FALSE,
  ncol_sample = 6
)
```

## Arguments

- seurat_obj:

  The Seurat object (must contain Has_BCR/Has_TCR cols)

- tissue_type:

  The tissue type of interest.

- airr_type:

  BCR or TCR

- clrs_specific:

  The specific color palette (should be named).

- plot_by:

  The grouping method: all samples, by dataset, by sample, by isotype.

- barcode_col:

  The barcode column name: barcode, cell_id, or Cell_ID_Unique.

- plot_label:

  Whether or not to include the labels.

- ncol_sample:

  The number of columns for sample-wise plots.

## Value

A Seurat UMAP plot with AIRR cells highlighted by the specified
grouping.

## Details

lightgray is Seurat's default background color and cells_total has to be
a list for the labels to work later. This could probably be replaced
with
[`plot_umap()`](https://eba28.github.io/athanor/reference/plot_umap.md).
