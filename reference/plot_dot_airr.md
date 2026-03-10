# Add AIRR (and other) info along the right side of an existing Seurat `DotPlot`

This function enhances dot plots by adding adaptive immune receptor
repertoire (AIRR) information and other metadata along the right side
for easy comparison. Currently supports cluster size, mean mutation
frequency, and BCR/TCR percentages.

## Usage

``` r
plot_dot_airr(
  plot,
  seurat_obj,
  row_identity = "seurat_clusters",
  facet_col,
  info_to_add = c("cluster_size")
)
```

## Arguments

- plot:

  The generated Seurat `DotPlot`.

- seurat_obj:

  The Seurat object containing the data.

- row_identity:

  The y axis identities.

- facet_col:

  The column to facet by e.g. "Cell_Type_Full".

- info_to_add:

  Vector of information types to add. The options are: `cluster_size`,
  `mean_mu_freq`, `percent_BCR`, `percent_TCR`, and `TRUST4`.

## Value

A Seurat DotPlot with annotations along the right side.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- DotPlot(seurat_obj, features = genes)
plot_dot_airr(p, seurat_obj, info_to_add = c("cluster_size", "mean_mu_freq"))
} # }
```
