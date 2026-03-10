# Give a summary of a Seurat object post-WNN

This function generates a summary message about the post-WNN Seurat
object, including the number of cells, details about the assays (e.g.,
number of genes, markers, embedding dimensions), information about the
reductions used for WNN, and the number of clusters identified in each
modality (RNA, BCR, and WNN) based on the largest resolutions.

## Usage

``` r
extract_wnn_vars(
  seurat_obj,
  gex_pca = "rpca",
  other_pca = "bpca",
  other_type = "BCR"
)
```

## Arguments

- seurat_obj:

  The post-WNN Seurat object

- gex_pca:

  The name of the GEX PCA reduction.

- other_pca:

  The name of the BCR/ADT/etc. PCA reduction.

- other_type:

  The second assay. Defaults to "BCR".

## Value

A text message.

## Details

Assumes that embeddings were used (for now) and that the object has RNA.
