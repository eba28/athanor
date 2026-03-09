# This function runs Weighted Nearest Neighbor (WNN) analysis on combined GEX and BCR data.

This function runs Weighted Nearest Neighbor (WNN) analysis on combined
GEX and BCR data.

## Usage

``` r
run_wnn(
  seurat_obj,
  embeddings,
  embedding_type,
  pc_gex = 20,
  pc_bcr = 20,
  k_param = 20,
  k_main = 20,
  cluster = TRUE,
  cluster_res = list(GEX = 1, BCR = 1, JC = 1),
  modality_weights = NULL,
  show_output = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing GEX data (at the least).

- embeddings:

  Matrix of BCR embeddings (genes x cells).

- embedding_type:

  The embeddings method.

- pc_gex:

  The number of PCs for the GEX assay.

- pc_bcr:

  The number of PCs for the BCR assay.

- k_param:

  The number of neighbors to use for each modality. Can be a single
  value or a vector of values to test.

- k_main:

  The main number of neighbors to use for the final WNN UMAP and
  clustering.

- cluster:

  Whether or not to perform clustering.

- cluster_res:

  Named list of clustering resolutions for GEX, BCR, and joint
  clustering (JC).

- modality_weights:

  Named vector of modality weights. If NULL, Seurat will calculate
  automatically.

- show_output:

  Whether or not to show verbose output from Seurat functions.

## Value

A Seurat object with WNN run.

## Details

- Currently only works for the embeddings approach and BCR data.

- The GEX object must have a `cell_id` metadata column.

- If I end up combining multiple embeddings into one object, then I
  should use something like
  `bcr_assay_name <- "BCR" # paste0("BCR_", embedding_type)`

- The neighbors for the assays can be saved in both the `graphs` slot
  (`compute.SNN` for clustering) and the `neighbors` slot
  (`return.neighbors` for distance calculations later) The `compute.SNN`
  option constructs a shared nearest neighbor graph using Jaccard index.
  Perhaps the GEX and BCR sections should be run if `modality_weights`
  is provided.
