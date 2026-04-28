# Build a Seurat object from BCR embeddings

Creates and processes a BCR-only Seurat object from a matrix of
pre-computed embeddings (e.g. AntiBERTa2, AntiBERTy, BALM-paired, ESM2,
immune2vec, etc.). Produces `bpca`, `BCR.nn`, `BCR_nn`, `BCR_snn`, and
`bcr.umap` reductions/graphs.

## Usage

``` r
bcr_embeddings_pipeline(
  embeddings,
  embedding_type,
  combined_airr = NULL,
  new_cols = NULL,
  num_pcs = 50,
  num_dims = 20,
  k_param = 20,
  verbose = TRUE
)
```

## Arguments

- embeddings:

  Matrix of BCR embeddings (features x cells).

- embedding_type:

  Character label for the embedding method (stored in `Misc`).

- combined_airr:

  Optional data frame passed to
  [`gex_add_airr()`](https://eba28.github.io/athanor/reference/gex_add_airr.md)
  to add AIRR metadata columns. If NULL, the step is skipped.

- new_cols:

  Character vector of columns to add from `combined_airr`. Only used
  when `combined_airr` is provided.

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of PCA dimensions to use for neighbor finding and UMAP.

- k_param:

  Number of nearest neighbors.

- verbose:

  Logical indicating whether to print progress messages.

## Value

A Seurat object with a BCR assay, a new PCA (`bpca`), new neighbor
graphs (`BCR_nn`, `BCR_snn`, `BCR.nn`), and a new UMAP (`bcr.umap`).

## Details

Embeddings are used as-is for `scale.data` (no `ScaleData` call) since
they are already on a comparable scale. The `data` layer is populated
from `counts` so downstream reads do not fail. It is possible that some
embeddings return identical values across all dimensions for different
cells. This can cause `RunUMAP()` to hang on the spectral initialization
step as it struggles to find a good low-dimensional representation of
the data.
