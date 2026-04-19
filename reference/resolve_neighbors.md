# Resolve the nearest-neighbor index matrix from a Seurat object

Resolve the nearest-neighbor index matrix from a Seurat object

## Usage

``` r
resolve_neighbors(seurat_obj, base_assay)
```

## Arguments

- seurat_obj:

  A Seurat object containing a neighbors slot.

- base_assay:

  Character string passed to
  [`map_assay_name()`](https://eba28.github.io/athanor/reference/map_assay_name.md)
  to determine which neighbor slot to retrieve.

## Value

An integer matrix of nearest-neighbor indices (cells by k).
