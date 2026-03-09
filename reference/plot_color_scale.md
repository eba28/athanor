# Generates a proper color scale for a Seurat DotPlot (with white at zero).

Generates a proper color scale for a Seurat DotPlot (with white at
zero).

## Usage

``` r
plot_color_scale(
  plot,
  palette = rev(pals::brewer.rdbu(n = 7)),
  val_col = "avg.exp.scaled",
  fill_by = "color"
)
```

## Arguments

- plot:

  The generated Seurat DotPlot or ggplot.

- palette:

  A palette of colors to go off of.

- val_col:

  The column in the plot data that contains the values to be plotted
  (e.g. "avg.exp.scaled").

- fill_by:

  One of "color" or "fill".

## Value

A Seurat dot plot or ggplot with an updated color scale.

## Details

Seurat's `col` option frequently is misleading with where the zeroes
fall. I don't want to rescale the expression. You could also just do
something like
`scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B")`
Can also be used to generated a color scale for a general ggplot2. Use
the function via a pipe right after the function call.
