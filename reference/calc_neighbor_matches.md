# Calculate neighbor matching scores across metadata columns

For each cell, calculates the proportion of its nearest neighbors that
share the same value for each metadata column.

## Usage

``` r
calc_neighbor_matches(
  seurat_obj,
  nn_names = names(seurat_obj@neighbors),
  meta_cols = c("annotated_clusters_bcr", "annotated_clusters_binary",
    "annotated_clusters_gex_bcr", "annotated_clusters_simpler", "cdr3_aa_length",
    "clone_id_unique", "isotype_stage", "locus_light", "mu_freq_bins_binary",
    "v_call_family"),
  cdr3_length_range = 1,
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

- cdr3_length_range:

  Integer range for `cdr3_aa_length` matching. Neighbors within this
  many amino acids of the query cell are counted as matches. Defaults to
  1.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell
  values.

- verbose:

  If TRUE, print progress messages.

## Value

Data frame with columns: Graph, Assay, Feature, Method, Score. If
`return_mean = TRUE`, `Score` is the mean across all cells; else it
contains per-cell values with an additional `cell_id` column.

## Details

`cell_id`s are saved for subsetting later if desired. Use
[`permute_neighbor_matches()`](https://eba28.github.io/athanor/reference/permute_neighbor_matches.md)
to compute a permuted baseline.
