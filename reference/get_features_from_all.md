# Get specific markers from the marker genes database

This function returns specific markers from the "all" marker genes
dataframe based on various filtering criteria including source, cell
types, and tissue types.

## Usage

``` r
get_features_from_all(
  markers_df,
  sources,
  contains,
  tissue_types,
  cell_types,
  alphabetize_types = TRUE,
  alphabetize_all = TRUE
)
```

## Arguments

- markers_df:

  The database of marker genes

- sources:

  Optional vector of sources - who the markers came from if you want
  specific origins.

- contains:

  Optional string to catch multiple cell types (e.g. "mDCs" and "pDCs").

- tissue_types:

  Optional vector of tissue types (e.g. "blood", "skin", etc.).

- cell_types:

  Optional vector of cell types you want markers for.

- alphabetize_types:

  Whether to return the markers alphabetized for each cell type.

- alphabetize_all:

  Whether to return all of the selected markers alphabetized.

## Value

A character vector of unique gene features/markers.
