# Plot several UMAPs side by side

This function generates multiple UMAP plots from a Seurat object with
various customizable options for coloring, labeling, and grouping the
data based on different metadata columns, cluster annotations, and
specific clusters of interest. It allows for a comprehensive overview of
the data across different embedding types (e.g. RNA, ADT, BCR) and
comparisons (e.g. annotated clusters, V call families, isotypes), making
it easier to explore and interpret the underlying structure of the data.

## Usage

``` r
plot_overview_comps(
  seurat_objs,
  data_source = "",
  pt_size = 0.1,
  second_assay = "BCR",
  assay_name,
  reduction = "wnn.umap",
  use_adt = FALSE,
  ncol,
  comparisons = c("annotated_clusters_simpler", "v_call_family", "light_chains",
    "isotype", "mu_freq"),
  details_col = "embedding_type",
  custom_colors = list(),
  ...
)
```

## Arguments

- seurat_objs:

  List of Seurat objects.

- data_source:

  Dataset description.

- pt_size:

  The size of the points in the UMAP.

- second_assay:

  The second assay to use in the title if plotting a combined reduction.
  Will usually be "BCR".

- assay_name:

  The name of the assay to use in the title. By default, it will be set
  based on the reduction (e.g. "GEX" for "rna.umap", "BCR" for
  "bcr.umap", and "GEX & BCR" for "wnn.umap").

- reduction:

  Which reduction to plot (e.g. "rpca", "bcr.umap", "wnn.umap").

- use_adt:

  Whether or not the comparisons being plotted represent ADT markers.

- ncol:

  The number of columns to use in the grid.

- comparisons:

  Which metadata columns to plot. By default, it will plot
  "annotated_clusters_simpler", "v_call_family", "light_chains",
  "isotype", and "mu_freq". The first one is the simplified CellTypist
  annotations, and the rest are BCR features.

- details_col:

  Which column in `seurat_obj@misc` to use for the plot subtitles. By
  default, it will use "embedding_type" to show the type of embedding
  being plotted (e.g. "AntiBERTa2").

- custom_colors:

  A named list of color scales for comparisons not covered by
  `named_colors`, keyed by comparison name (e.g.
  `list(cdr3_bins = my_colors)`). Takes precedence over `named_colors`
  when both define the same comparison.

## Value

A patchwork object with overview plots in a grid.

## Details

"L %\|\|% R newly in base is an expressive idiom for the phrases
if(!is.null(L)) L else R or if(is.null(L)) R else L." - R 4.4
