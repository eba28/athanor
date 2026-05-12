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
  n_permutations = 10,
  return_mean = TRUE,
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

- n_permutations:

  Number of permutations to run.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell values
  (averaged across permutations).

- verbose:

  If TRUE, print progress messages.

## Value

Data frame in the same format as
[`calc_neighbor_matches()`](https://eba28.github.io/athanor/reference/calc_neighbor_matches.md).
