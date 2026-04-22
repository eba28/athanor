# Calculate percentages per metadata group

This function calculates the percentage of occurrences of a specified
focus group within a dataset, grouped by specified metadata columns. It
can be used in conjunction with
[`plot_pcts()`](https://eba28.github.io/athanor/reference/plot_pcts.md)
to visualize the results.

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

## Examples

``` r
df <- data.frame(
  sample_id = rep(c("S1", "S2"), each = 50),
  cell_type = c(rep(c("B", "T", "NK"), times = c(20, 20, 10)),
                rep(c("B", "T", "NK"), times = c(10, 30, 10)))
)
calc_pcts(df, meta_group_by = "sample_id", focus_group = "cell_type")
#> # A tibble: 6 × 4
#> # Groups:   sample_id [2]
#>   sample_id cell_type Count Percent
#>   <chr>     <chr>     <int>   <dbl>
#> 1 S1        B            20      40
#> 2 S1        NK           10      20
#> 3 S1        T            20      40
#> 4 S2        B            10      20
#> 5 S2        NK           10      20
#> 6 S2        T            30      60
```
