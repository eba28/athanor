# Generate a color scale for a ggplot2 plot with white at zero

It creates a gradient of colors based on the range of values in the plot
or data provided, and applies it to the specified color or fill
aesthetic if applicable.

## Usage

``` r
plot_color_scale(
  plot,
  data,
  val_col = "avg.exp.scaled",
  palette = rev(pals::brewer.rdbu(n = 7)),
  fill_by = "color"
)
```

## Arguments

- plot:

  The generated Seurat DotPlot or ggplot.

- data:

  The data used to generate the plot.

- val_col:

  The column in the plot data that contains the values to be plotted
  (e.g. "avg.exp.scaled").

- palette:

  A palette of colors to go off of.

- fill_by:

  One of "color" or "fill".

## Value

A Seurat dot plot or ggplot with an updated color scale, or just a
vector of colors if a plot is not provided.

## Details

Seurat's `col` option frequently is misleading with where the zeroes
fall. Works well for making a scale for a Seurat `DotPlot` that
accurately reflects the expression value. You could also just do
something like
`scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B")`
Can also be used to generated a color scale for a general ggplot2. Use
the function via a pipe right after the function call.

## Examples

``` r
p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, color = mpg)) +
  ggplot2::geom_point()
plot_color_scale(plot = p, val_col = "mpg")
```
