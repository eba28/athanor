# Add AIRR information to a Seurat object

This function integrates adaptive immune receptor repertoire (AIRR) data
with gene expression (GEX) data in a Seurat object. Currently built to
integrate BCR data, including both heavy and light chain information.

## Usage

``` r
gex_add_airr(
  seurat_obj,
  airr_type = "BCR",
  combined_airr,
  new_cols,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  The Seurat object containing GEX data.

- airr_type:

  Type of immune receptor data. Currently supports "BCR".

- combined_airr:

  BCR AIRR formatted data frame with heavy and light chains.

- new_cols:

  Vector of column names to select from the AIRR data.

- verbose:

  Logical indicating whether or not to print integration summary
  information.

## Value

The Seurat object with AIRR columns added to the metadata, including
Has_BCR, isotype information, mutation frequencies, and pairing status.

## Details

This would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
and before
[`concatenate_gex_bcr()`](https://eba28.github.io/athanor/reference/concatenate_gex_bcr.md).
Right now this is just built to integrate in BCR data and assumes that
the BCR data includes light chains. Assumes that `seurat_obj` contains
`cell_id` and `annotated_clusters_simpler`. Collapsed light chains are
alphabetized. There are columns for IGK and IGL instead of just "light"
since a cell can have both. Adds mutation frequency bins and isotype
staging information for BCR data.
