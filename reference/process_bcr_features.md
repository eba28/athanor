# Process BCR features for integration

Processes BCR features by renaming columns, converting
ordered/categorical variables to numeric representations, and
normalizing all features. This prepares BCR data for integration with
gene expression data in a Seurat object.

## Usage

``` r
process_bcr_features(bcr_features)
```

## Arguments

- bcr_features:

  A data frame containing BCR features (e.g. isotype, mutation
  frequency). May contain ordered factors, numeric, or categorical
  variables.

## Value

A transposed matrix of processed features where:

- Ordered variables are converted to numeric scores and suffixed with
  "-ordered"

- Numeric variables are suffixed with "-scaled"

- Categorical variables are one-hot encoded (all categories kept)

- All numeric predictors are normalized to mean = 0 and sd = 1

- Underscores are removed from feature names (so Seurat doesn't throw a
  warning)

## Details

Uses `recipes` to:

1.  Rename columns to distinguish from existing metadata

2.  Convert ordered factors to ordinal scores OR one-hot encode
    categoricals

3.  Center and normalize all numeric predictors

4.  Transpose the result for compatibility with Seurat assays You could
    do something like janitor::clean_names(bcr_features) to remove the
    underscores from the column names right off the bat, but one-hot
    encoding will add names with underscores automatically unless you
    messing with the `naming` argument in `step_dummy()`.
