# Calculate the proportion of neighbors that match metadata categories or ADT expression thresholds

Calculate the proportion of neighbors that match metadata categories or
ADT expression thresholds

## Usage

``` r
calc_neighbor_matches(
  seurat_obj,
  nn_name,
  meta_cols = c("annotated_clusters_bcr", "annotated_clusters_binary",
    "annotated_clusters_gex_bcr", "annotated_clusters_simpler", "cdr3_aa_length",
    "clone_id_unique", "isotype_stage", "locus_light", "mu_freq_bins_binary",
    "v_call_family"),
  adt_features = NULL,
  adt_range = 0.1,
  adt_methods = c("mean_abs", "range"),
  permute = FALSE,
  n_permutations = 10,
  previous_matches,
  path_save
)
```

## Arguments

- seurat_obj:

  The Seurat object, with details added to the Misc() slot.

- nn_name:

  Name of the nearest neighbor graph slot.

- meta_cols:

  Character vector of metadata columns to evaluate.

- adt_features:

  Character vector of ADT feature names to evaluate.

- adt_range:

  Numeric vector of ADT expression range threshold(s).

- adt_methods:

  How to calculate the ADT metric(s).

- permute:

  Shuffle labels for each meta column and ADT expression per feature
  before computing matches to get a random baseline.

- n_permutations:

  The number of times to permute labels.

- previous_matches:

  Data frame of previous matches to combine with the new results
  (optional).

- path_save:

  Where to save the results.

## Value

Data frame with these columns: Full_Name, Category, Category_Details,
Assay, Meta_Col, Method, Matches

## Details

`cell_id`s are saved for subsetting later if desired.
