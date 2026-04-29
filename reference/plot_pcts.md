# Plot percentages in a stacked bar plot

This function plots a stacked bar plot of percentages calculated using
[`calc_pcts()`](https://eba28.github.io/athanor/reference/calc_pcts.md)
with percentages labeled and total counts on top.

## Usage

``` r
plot_pcts(
  pcts,
  data_source,
  clrs_specific,
  plot_type = "All",
  plot_value = "Cell Type",
  x_axis = "sample_id",
  x_axis_label = "Sample",
  fill_type = "annotated_clusters",
  fill_label = fill_type,
  perc_min = 3,
  label_size = 3,
  label_fill = FALSE,
  include_counts = TRUE,
  drop_zeroes = TRUE,
  reverse_order = FALSE,
  total_order = FALSE,
  details = NULL
)
```

## Arguments

- pcts:

  The output of
  [`calc_pcts()`](https://eba28.github.io/athanor/reference/calc_pcts.md).

- data_source:

  Dataset description.

- clrs_specific:

  The specific color palette (should be named).

- plot_type:

  One of `All`, `Binary`.

- plot_value:

  What is being plotted.

- x_axis:

  What to put along the x axis.

- x_axis_label:

  The label for the x axis.

- fill_type:

  What to group the bar plot by.

- fill_label:

  The description of what you're filling by.

- perc_min:

  The minimum percentage to show in the plot.

- label_size:

  The size of the percentage labels.

- label_fill:

  Add a white background for clarity.

- include_counts:

  Plot the counts on top/bottom.

- drop_zeroes:

  Remove percentages of zeroes.

- reverse_order:

  Change the fill order for a stacked plot.

- total_order:

  Rearrange x axis in descending order by totals instead of
  alphabetically.

- details:

  An optional custom subtitle.

## Value

A stacked ggplot bar plot

## Details

Give the percentages already as percents (\* 100 in the calculations).
Make sure `pcts` includes Dataset if you want to split by dataset. Note
that the percentages may seem inaccurate because of the accuracy. There
was a big issue with the combo of `geom_text()` & `facet_grid()`. This
assumes that you want to show the counts for binary plots.
