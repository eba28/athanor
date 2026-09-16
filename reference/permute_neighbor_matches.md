# Compute a permuted baseline for neighbor matching scores

Runs
[`calc_neighbor_matches()`](https://eba28.github.io/athanor/reference/calc_neighbor_matches.md)
`n_permutations` times with shuffled metadata labels and returns the
mean score across permutations as a random baseline.

## Usage

``` r
permute_neighbor_matches(
  seurat_obj,
  nn_names = names(seurat_obj@neighbors),
  meta_cols = c("annotated_clusters_bcr", "annotated_clusters_binary",
    "annotated_clusters_gex_bcr", "annotated_clusters_simpler", "cdr3_aa_length",
    "clone_id_unique", "isotype_stage", "locus_light", "mu_freq_bins_binary",
    "v_call_family"),
  gex_cols = c("annotated_clusters_binary", "annotated_clusters_gex_bcr",
    "annotated_clusters_simpler", "seurat_clusters"),
  n_permutations = 10,
  return_mean = TRUE,
  drop_self = TRUE,
  verbose = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- nn_names:

  Character vector of neighbor graph names to evaluate. Defaults to all
  graphs in the object.

- meta_cols:

  Character vector of metadata columns to evaluate.

- gex_cols:

  Character vector of `meta_cols` that describe GEX-derived groupings;
  these get `Assay = "GEX"` in the output and every other column gets
  `Assay = "BCR"`. Defaults to the known GEX cluster columns plus
  `"seurat_clusters"`.

- n_permutations:

  Number of permutations to run.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell values
  (averaged across permutations).

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) so `k` is fixed
  across cells regardless of whether self-matches are present. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

- verbose:

  If TRUE, print progress messages.

## Value

Data frame in the same format as
[`calc_neighbor_matches()`](https://eba28.github.io/athanor/reference/calc_neighbor_matches.md).
