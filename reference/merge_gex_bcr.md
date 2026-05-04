# Merge a GEX Seurat object with a BCR Seurat object

Combines a gene expression (GEX) Seurat object and a BCR Seurat object
by matching on shared `cell_id` values. Handles the common case where
the two objects have different cell counts (e.g. not every GEX cell has
a BCR sequence).

## Usage

``` r
merge_gex_bcr(gex_obj, bcr_obj, transfer_reductions = TRUE, verbose = TRUE)
```

## Arguments

- gex_obj:

  A Seurat object containing GEX (RNA) data with a `cell_id` metadata
  column.

- bcr_obj:

  A Seurat object containing a BCR assay with a `cell_id` metadata
  column (typically produced by
  [`bcr_embeddings_pipeline()`](https://eba28.github.io/athanor/reference/bcr_embeddings_pipeline.md)).

- transfer_reductions:

  Whether to copy BCR reductions (`bpca`, `bcr.umap`) and graphs
  (`BCR.nn`, `BCR_nn`, `BCR_snn`) from `bcr_obj` into the merged object.

- verbose:

  Whether or not to print a summary of the merge.

## Value

A Seurat object with both RNA and BCR assays and BCR metadata columns.

## Details

Cell matching is done via `cell_id` metadata, not raw Seurat barcodes.
If GEX barcodes differ from `cell_id` (e.g. they carry a sample prefix),
BCR cell names are renamed via
[`Seurat::RenameCells()`](https://satijalab.github.io/seurat-object/reference/RenameCells.html)
to match before the assay is transferred.

Both objects are subset to their shared cells; the result is suitable
for concatenation or WNN or other downstream workflows.

Due to subsetting, graphs, neighbors and reductions maybe have to be
regenerated.
