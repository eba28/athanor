# This function gives a summary of the post-WNN object.

This function gives a summary of the post-WNN object.

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
