# Extract and validate an ADT matrix from a Seurat object

Extract and validate an ADT matrix from a Seurat object

## Usage

``` r
get_adt_matrix(
  seurat_obj,
  adt_assay = "ADT",
  layer = "data",
  features_adt = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- adt_assay:

  Name of the ADT assay.

- layer:

  Layer to retrieve from the assay.

- features_adt:

  Optional character vector of ADT features to subset. If `NULL`, all
  features are returned.

## Value

A numeric matrix of dimension cells by features.

## Details

The Seurat object must contain an assay with the specified `adt_assay`
name, and that assay must contain a layer with the specified `layer`
name. If `features_adt` is provided, it must be a character vector of
feature names that are present in the assay. The function will return a
matrix of dimension cells by features, where the rows are named with
cell IDs and the columns are named with feature names.
