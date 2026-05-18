# Concatenate GEX and BCR data in a Seurat object

Creates a combined representation of gene expression (GEX) and B-cell
receptor (BCR) data, then runs the standard Seurat pipeline on the
result. Supports two input types and two PCA stages (see Details).

## Usage

``` r
concatenate_gex_bcr(
  seurat_obj,
  pca_stage = c("raw", "reduced", "embed"),
  input_type = c("features", "embeddings"),
  cols_to_include,
  embeddings = NULL,
  gex_reduction = "rpca",
  filter_genes,
  ensembl_version = NULL,
  cache_file = NULL,
  normalize = TRUE,
  num_features = 2000,
  num_pcs = 50,
  num_dims = 20,
  k_param = 20,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object with RNA assay and BCR metadata.

- pca_stage:

  One of `"raw"`, `"reduced"`, or `"embed"`. `"raw"` concatenates at the
  count/feature level and runs a joint PCA from scratch. `"reduced"`
  appends BCR features onto transposed GEX PCA embeddings and runs a
  joint PCA. `"embed"` column-binds existing GEX and BCR PCA spaces
  directly. See Details for all combinations.

- input_type:

  `"features"` to use processed BCR metadata columns; `"embeddings"` to
  use a pre-computed embedding matrix or the BCR assay from a merged
  object.

- cols_to_include:

  Character vector of BCR metadata column names to use as features (e.g.
  `c("mu_freq", "isotype")`). Required for `input_type = "features"`.

- embeddings:

  A features-by-cells matrix of BCR embeddings. Required for
  `input_type = "embeddings"` when no merged object is provided.

- gex_reduction:

  Name of the GEX PCA reduction to use as the GEX component for
  `pca_stage = "reduced"` and `pca_stage = "embed"`. Defaults to
  `"rpca"`. Use `"integrated"` to use a batch-corrected reduction (e.g.
  from Harmony).

- filter_genes:

  If specified, filter out genes from this category (e.g. `"IG"` and/or
  `"TR"`). Only applies for `pca_stage = "raw"`.

- ensembl_version:

  Ensembl version for gene annotations (e.g. `"GRCh38.104"`). If `NULL`,
  auto-detected from `Misc(seurat_obj, "ensembl_version")` when
  available.

- cache_file:

  Passed to
  [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md).
  Path to a cached RDS result to use instead of querying Ensembl live.

- normalize:

  If `TRUE`, normalize the `RNA_BCR` assay before scaling. Set to
  `FALSE` if the data has already been normalized. Ignored for
  `pca_stage = "reduced"` (GEX PCs do not need normalization).

- num_features:

  Number of variable features for the `"Before"` paths.

- num_pcs:

  Number of principal components to compute.

- num_dims:

  Number of PCA dimensions to use for neighbor finding. For `"After"`
  and `"Before_PCs"`, also controls how many GEX PCs are taken from
  `rpca` before concatenation.

- k_param:

  Number of nearest neighbors.

- verbose:

  Whether to show output from Seurat functions.

## Value

A Seurat object with:

- New `RNA_BCR` assay (for `pca_stage = "Before"`)

- PCA reduction (`rna_bcr.pca`)

- UMAP reduction (`rna_bcr.umap`)

- Neighbor graphs computed on the combined data

## Details

This would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
and
[`gex_add_airr()`](https://eba28.github.io/athanor/reference/gex_add_airr.md).

The six combinations of `input_type` and `pca_stage` offer different
tradeoffs:

**`input_type = "embeddings"`, `pca_stage = "raw"`**  
Same as above but BCR data comes from a pre-computed embedding matrix or
the `BCR` assay of a merged object (from
[`merge_gex_bcr()`](https://eba28.github.io/athanor/reference/merge_gex_bcr.md)).
Supply `embeddings` directly or pass a merged object and the BCR assay
is detected automatically. Embedding dimensions tend to have more
comparable scales to RNA features than raw metadata columns, but scale
mismatch still applies.

**`input_type = "features"`, `pca_stage = "raw"`** (default)  
BCR metadata columns (via `cols_to_include`) are processed by
[`process_bcr_features()`](https://eba28.github.io/athanor/reference/process_bcr_features.md)
into a numeric features-by-cells matrix and row-bound onto the RNA
count/data matrix. A new `RNA_BCR` assay is created and the full Seurat
pipeline (normalize, scale, PCA, neighbors, UMAP) is run from scratch on
the combined data. The main limitation is scale mismatch: log-normalized
RNA values and BCR metadata live in different ranges, so PCA may be
dominated by whichever modality has higher total variance even after
scaling.

**`input_type = "embeddings"`, `pca_stage = "reduced"`**  
Same as the features variant above but uses BCR embeddings instead of
metadata columns. BCR embedding dimensions are appended to transposed
GEX PCA embeddings and a joint PCA is run on the combined matrix.

**`input_type = "features"`, `pca_stage = "reduced"`**  
A middle ground between `"raw"` and `"embed"`. Instead of row-binding
BCR features onto the raw RNA matrix, they are appended to the
transposed GEX PCA embeddings (`rpca`, subset to `num_dims` PCs). The
combined matrix (n_gex_pcs + n_bcr_features rows) is stored as a new
`RNA_BCR` assay, scaled, and a joint PCA is run. This avoids the scale
mismatch of `"raw"` (GEX PCs and BCR metadata are in more comparable
ranges) while still performing a new projection that can mix the two
modalities. Normalization is skipped since GEX PCs are already
processed. `filter_genes` does not apply.

**`input_type = "embeddings"`, `pca_stage = "embed"`**  
Column-binds the existing `rpca` and `bpca` reductions directly.
Requires a merged object from
[`merge_gex_bcr()`](https://eba28.github.io/athanor/reference/merge_gex_bcr.md)
with both reductions already computed. This is the most efficient path
when BCR embeddings are already available. As with the features "embed"
path, `num_dims` controls how many PCs are taken from each reduction
before joining.

**`input_type = "features"`, `pca_stage = "embed"`**  
BCR metadata features are first embedded into their own PCA space: a
`BCR` assay is created, scaled, and PCA is run to produce `bpca` (capped
at `nrow(bcr_features) - 1` PCs). The resulting BCR PCA embeddings are
then column-bound with the existing `rpca` embeddings to form a joint
PCA space. This is more principled than `"raw"` for metadata features
because both modalities are in comparable PCA spaces before being
joined, and each modality's internal variance structure is preserved.
Use `num_dims` to control how many PCs are taken from each side.

For the `"raw"` and `"reduced"` paths, variable features are always set
by appending BCR features onto the existing RNA variable features rather
than re-running
[`Seurat::FindVariableFeatures()`](https://satijalab.org/seurat/reference/FindVariableFeatures.html).
Re-running would drop BCR features from the selection and trigger Seurat
warnings about underscores in feature names and missing count layers.

The `normalize`, `num_features`, `num_pcs`, `num_dims`, and `k_param`
arguments are passed to
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
for the `"raw"` and `"reduced"` paths. For `"reduced"` and `"embed"`,
`num_dims` additionally controls how many GEX PCs are taken from `rpca`
before concatenation.
