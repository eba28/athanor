# Create bar plots of read or cell counts for quality control

This function generates bar plots showing read counts, cell counts, or
barcode counts for quality control purposes. Can display data split by
data type, sample, or other grouping variables.

## Usage

``` r
plot_counts(
  summary_df,
  data_types = "All",
  count_type,
  aggregation_state = "Unaggregated",
  fill_type_label = "Data Type",
  x_axis = "sample_id",
  x_axis_label = "Sample",
  clrs_datatype
)
```

## Arguments

- summary_df:

  Output from create_metrics_summary().

- data_types:

  Vector of data types to include (e.g. "BCR", "GEX", "TCR"). Use "All"
  to include all data types.

- count_type:

  Type of count to plot: "Read", "Cell", or "Barcode".

- aggregation_state:

  Whether data is "Aggregated" or "Unaggregated". (from cellranger
  aggr).

- fill_type_label:

  Label for the fill aesthetic.

- x_axis:

  Variable to plot along the x-axis.

- x_axis_label:

  Label for the x-axis.

- clrs_datatype:

  Named vector of colors for data types. Must be given in the correct
  order if plotting barcodes.

## Value

A ggplot2 barplot.
