# Build a Seurat object from BCR embeddings

Creates and processes a BCR-only Seurat object from a matrix of
pre-computed embeddings (e.g. ESM2, immune2vec). Produces `bpca`,
`BCR.nn`, `BCR_nn`, `BCR_snn`, and `bcr.umap` reductions/graphs.

## Usage

``` r
bcr_embeddings_pipeline(
  embeddings,
  embedding_type,
  num_pcs = 50,
  num_dims = 20,
  k_param = 20,
  combined_airr = NULL,
  airr_cols = NULL,
  verbose = TRUE
)
```

## Arguments

- embeddings:

  Matrix of BCR embeddings (features x cells).

- embedding_type:

  Character label for the embedding method (stored in `Misc`).

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of PCA dimensions to use for neighbor finding and UMAP.

- k_param:

  Number of nearest neighbors.

- combined_airr:

  Optional data frame passed to
  [`gex_add_airr()`](https://eba28.github.io/athanor/reference/gex_add_airr.md)
  to add AIRR metadata columns. If NULL, the step is skipped.

- airr_cols:

  Character vector of columns to add from `combined_airr`. Only used
  when `combined_airr` is provided.

- verbose:

  Logical indicating whether to print progress messages.

## Value

A Seurat object with BCR assay, PCA (`bpca`), neighbor graphs (`BCR_nn`,
`BCR_snn`, `BCR.nn`), and UMAP (`bcr.umap`).

## Details

Embeddings are used as-is for `scale.data` (no `ScaleData` call) since
they are already on a comparable scale. The `data` layer is populated
from `counts` so downstream reads do not fail.
