# Gets the information needed using `add_info_bar()` on a `DotPlot`

This function takes the filtered markers dataframe and selects the
relevant columns to be used for adding an info bar to a DotPlot. It
renames the Marker column to features.plot and formats the
Cell_Type_Full column for better readability.

## Usage

``` r
get_cell_types(markers_df)
```

## Arguments

- markers_df:

  The markers data.frame filtered to match your input features

## Value

A data.frame with Cell_Type_Full and features.plot columns
