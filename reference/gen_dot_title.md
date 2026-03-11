# Generate a title for a `DotPlot`

This function generates a title for a DotPlot by combining a provided
dataset description with a formatted list of marker sources. The marker
sources are sorted, formatted to replace underscores with spaces, and
concatenated into a comma-separated string enclosed in parentheses.

## Usage

``` r
gen_dot_title(plot_title = "", marker_sources)
```

## Arguments

- plot_title:

  Dataset description

- marker_sources:

  The list of marker sources

## Value

A string with sources comma-separated and in parentheses
