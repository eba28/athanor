# Display a nicely formatted table when the R Markdown file knits

This function creates a formatted table with scrolling capability for
use in R Markdown documents. It applies striped styling and makes the
table scrollable within specified dimensions.

## Usage

``` r
print_kable(table, kable_height = "500px", kable_width = "100%")
```

## Arguments

- table:

  The input data frame to be printed

- kable_height:

  The height of the output table (you can set it to NULL to display the
  full table without scrolling). Default is "500px".

- kable_width:

  The width of the output table. Default is "100%".

## Value

A formatted table.
