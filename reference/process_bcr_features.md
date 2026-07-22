# Process BCR features for integration

Processes BCR features by renaming columns, converting
ordered/categorical variables to numeric representations, and
normalizing all features. This prepares BCR data for integration with
gene expression data in a Seurat object.

## Usage

``` r
process_bcr_features(
  bcr_features,
  scaling = c("z_score", "scale", "center", "range", "none"),
  scale_dummies = TRUE,
  remove_nzv = FALSE,
  verbose = TRUE
)
```

## Arguments

- bcr_features:

  A data frame containing BCR features (e.g. isotype, mutation
  frequency). May contain ordered factors, numeric, or categorical
  variables.

- scaling:

  Character string specifying how to scale numeric predictors. One of:

  - `"z_score"` (default): center to mean = 0 and scale to sd = 1
    (z-score).

  - `"scale"`: scale to sd = 1 without centering.

  - `"center"`: center to mean = 0 without scaling.

  - `"range"`: min-max normalization to \[0, 1\].

  - `"none"`: no scaling applied.

- scale_dummies:

  Logical. If `TRUE` (default), scaling is also applied to one-hot
  encoded dummy variables. If `FALSE`, scaling is applied before
  `step_dummy()` so that dummy columns are left on their original 0/1
  scale.

- remove_nzv:

  Logical. If `TRUE`, near-zero-variance predictors (as defined by
  [`recipes::step_nzv()`](https://recipes.tidymodels.org/reference/step_nzv.html))
  are dropped in addition to exact zero-variance ones. Off by default
  since it can remove rare-but-informative categories (e.g. a rare
  isotype).

- verbose:

  Logical indicating whether or not to print messages.

## Value

A transposed matrix of processed features where:

- Ordered variables are converted to numeric scores and suffixed with
  "-ordered"

- Numeric variables are suffixed with "-scaled"

- Categorical variables are one-hot encoded (all categories kept)

- All numeric predictors are scaled according to `scaling`

- Underscores are removed from feature names (so Seurat doesn't throw a
  warning)

## Details

Uses `recipes` to:

1.  Rename columns to distinguish from existing metadata

2.  Convert ordered factors to ordinal scores OR one-hot encode
    categorical variables

3.  Scale all numeric predictors according to `scaling`

4.  Optionally remove near-zero-variance predictors

5.  Transpose the result for compatibility with Seurat assays You could
    do something like janitor::clean_names(bcr_features) to remove the
    underscores from the column names right off the bat, but one-hot
    encoding will add names with underscores automatically unless you
    messing with the `naming` argument in `step_dummy()`.
