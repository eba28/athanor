# Build a Seurat object from BCR embeddings, deduplicating identical embeddings first

Wrapper around
[`bcr_embeddings_pipeline()`](https://eba28.github.io/athanor/reference/bcr_embeddings_pipeline.md)
that detects cells with identical embeddings (e.g. from clonal
expansion), runs PCA, neighbor finding, and UMAP on unique embeddings
only, then copies coordinates back to all cells. This avoids the
spectral initialization hang in `RunUMAP()` caused by degenerate
neighbor graphs from zero-distance duplicate points.

## Usage

``` r
bcr_embeddings_pipeline_dedup(
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

  Logical indicating whether or not to print messages.

## Value

A Seurat object as returned by
[`bcr_embeddings_pipeline()`](https://eba28.github.io/athanor/reference/bcr_embeddings_pipeline.md),
with all original cells present. Cells with identical embeddings receive
the same PCA and UMAP coordinates as their first occurrence.
