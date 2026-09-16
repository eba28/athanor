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
  cdr3_length_range = 0,
  clonal_cols = character(0),
  min_clone_size = 2,
  gex_cols = c("annotated_clusters_binary", "annotated_clusters_gex_bcr",
    "annotated_clusters_simpler", "seurat_clusters"),
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

- cdr3_length_range:

  Integer range for `cdr3_aa_length` matching. Neighbors within this
  many amino acids of the query cell are counted as matches.

- clonal_cols:

  Character vector of `meta_cols` (typically `"clone_id_unique"`) that
  should be scored as clone-size-weighted matches instead of the usual
  proportion-of-k. For these columns, the score is the number of
  matching neighbors divided by `min(clone_size - 1, k)`, i.e.
  normalized by the number of same-clone cells actually available to be
  found (capped at k) rather than by k itself. This avoids penalizing
  small clones that can never fill all k neighbor slots. Cells whose
  clone is smaller than `min_clone_size` get a score of NA.

- min_clone_size:

  Minimum clone size (inclusive) for a cell to receive a score on a
  `clonal_cols` column; smaller clones get NA. Ignored for columns not
  listed in `clonal_cols`.

- gex_cols:

  Character vector of `meta_cols` that describe GEX-derived groupings
  (e.g. transcriptomic clusters). These get `Assay = "GEX"` in the
  output; every other column gets `Assay = "BCR"`. Defaults to the known
  GEX cluster columns plus `"seurat_clusters"`.

- return_mean:

  If TRUE, return the mean across all cells; else return per-cell
  values.

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) so `k` is fixed
  across cells regardless of whether self-matches are present. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

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
