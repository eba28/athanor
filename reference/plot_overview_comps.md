# Plots several UMAPs side by side for a Seurat object.

Plots several UMAPs side by side for a Seurat object.

## Usage

``` r
plot_overview_comps(
  seurat_objs,
  data_source = "",
  pt_size = 0.1,
  second_assay = "BCR",
  assay_name,
  reduction = "wnn.umap",
  comparisons = c("annotated_clusters", "v_call_family", "light_chains", "isotype",
    "mu_freq")
)
```

## Arguments

- seurat_objs:

  List of WNN objects with different embedding types.

- data_source:

  Dataset description.

- pt_size:

  The size of the points in the UMAP.

- second_assay:

  The second assay to use in the title if plotting a WNN reduction.

- assay_name:

  The name of the assay to use in the title. By default, it will be set
  based on the reduction (e.g. "GEX" for "rna.umap", "BCR" for
  "bcr.umap", and "GEX & BCR" for "wnn.umap").

- reduction:

  Which reduction to plot ("rna.umap", "bcr.umap", or "wnn.umap").

- comparisons:

  The labelling of the plots.

## Value

patchwork object with overview plots

## Details

The names of the seurat_objs list correspond with embedding_types.
Assumes that CellTypist is the annotation approach being used.
data_source is set to empty to save on space
