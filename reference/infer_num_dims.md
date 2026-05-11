# Infer num_dims (PCA dims used for neighbor finding) from existing RNA neighbors.

Infers the number of dimensions used for neighbor finding from the
existing RNA neighbors slot, or returns a default if not found.

## Usage

``` r
infer_num_dims(seurat_obj, default = 20, verbose = TRUE)
```

## Arguments

- seurat_obj:

  A Seurat object containing a neighbors slot.

- default:

  The default number of dimensions to return if not found in the
  neighbors slot.

- verbose:

  Logical indicating whether or not to print a message about the source
  of the inferred num_dims.

## Value

An integer representing the number of dimensions used for neighbor
finding.
