# This function plots the internal/external clustering metrics.

This function plots the internal/external clustering metrics.

## Usage

``` r
plot_metrics(
  metrics,
  plot_title = "",
  best_score = "higher",
  type = "Internal",
  y_axis = "Labeling",
  round_to = 2,
  details = ""
)
```

## Arguments

- metrics:

  Data frame of metrics to plot, with columns: Embedding, Reduction,
  Labeling, Metric, Score.

- plot_title:

  Title to use for the plot.

- best_score:

  One of "higher" or "lower" to indicate whether higher or lower scores
  are better for the metrics being plotted. This is used to determine
  which scores to outline in the plot.

- type:

  One of "Internal" or "External" to indicate the type of metrics being
  plotted, used for the plot title.

- y_axis:

  The name of the column in the metrics data frame to use for the y-axis
  (e.g., "Labeling").

- round_to:

  Number of decimal places to round the score labels to in the plot.

- details:

  Additional details to include in the plot title (optional).

## Value

A ggplot showing the metrics across embeddings and reductions, with the
best scores outlined.

## Details

Make sure to check that the best score is consistent across your
plotting metrics.
