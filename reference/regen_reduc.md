# Regenerate neighbor graphs and UMAPs.

Regenerates PCA, neighbor graphs, and UMAP for a Seurat object. This is
useful if you subset a Seurat object and need to recompute these
reductions for the new object. It can also be used to regenerate these
reductions if you have modified the object in a way that requires them
to be redone (e.g. filtering out cells or genes).

## Usage

``` r
regen_reduc(
  seurat_obj,
  pca_name = "rpca",
  assay = "RNA",
  num_dims = 20,
  k_param = 20,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- pca_name:

  Name of the PCA reduction to use for neighbor finding and UMAP. This
  should be the name of an existing PCA reduction in the Seurat object
  (e.g. "pca" or "rpca").

- assay:

  Name of the assay to use for neighbor finding and UMAP.

- num_dims:

  Number of PCA dimensions to use for neighbor finding.

- k_param:

  Number of nearest neighbors.

- verbose:

  Print out Seurat's progress messages.

## Value

A processed Seurat object with the `graphs`, `neighbors`, and
`reductions` slots filled in or updated.

## Details

If you are providing a reduction based on batch effect integration (e.g.
Harmony, RPCA), then you should use the name of that reduction for
`pca_name`. Although this function is called "regen_reduc", it can also
be run to generate reductions for the first time (e.g. if you have
already run PCA and just want to calculate neighbors and run UMAP).
