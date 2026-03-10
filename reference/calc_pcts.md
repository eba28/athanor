# Calculate percentages per metadata group

This function calculates the percentage of occurrences of a specified
focus group within a dataset, grouped by specified metadata columns.

## Usage

``` r
calc_pcts(
  data,
  meta_group_by = c("sample_id", "Dataset"),
  focus_group,
  order_by
)
```

## Arguments

- data:

  A data.frame containing the data to be analyzed e.g. the metadata of a
  Seurat object.

- meta_group_by:

  A character vector specifying the columns to group by.

- focus_group:

  A character string specifying the column to focus on for percentage
  calculation.

- order_by:

  An optional character string specifying the focus group value to order
  the results by. If not provided, no ordering is applied.

## Value

A data.frame with the calculated percentages and counts for each group.

## Details

This is not very elegant, but seems to work (dealing with string args is
weird). Make sure to filter your data as needed (e.g. no NA isotypes)
beforehand. Made for cell types and isotypes. The percentages will be
doubles between 0 and 100. It will fill in missing values (especially
needed for ordering if provided). You could alternatively use the fill
option in geom_bar() instead.
