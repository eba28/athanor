# Infer k_param from existing RNA neighbors, or return a default.

Infers the k_param used for neighbor finding from the existing RNA
neighbors slot, or returns a default if not found.

## Usage

``` r
infer_k_param(seurat_obj, default = 20, verbose = TRUE)
```

## Arguments

- seurat_obj:

  A Seurat object containing a neighbors slot.

- default:

  The default k_param to return if not found in the neighbors slot.

- verbose:

  Logical indicating whether or not to print a message about the source
  of the inferred k_param.

## Value

An integer representing the k_param used for neighbor finding.
