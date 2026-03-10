# Display markers from a filtered marker genes database as a table

This function takes a filtered marker genes dataframe and returns a
formatted table showing markers organized by cell type.

## Usage

``` r
source_markers(filtered_markers_df)
```

## Arguments

- filtered_markers_df:

  A filtered dataframe from the marker genes database containing at
  least Cell_Type and Marker columns.

## Value

A formatted table showing markers grouped by cell type.
