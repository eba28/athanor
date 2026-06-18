# Calculate shared neighbors between two neighbor graphs

For each cell in a Seurat object, counts how many of its nearest
neighbors are shared between two specified neighbor graphs.

## Usage

``` r
calc_neighbor_overlap(
  seurat_obj,
  nn1,
  nn2,
  exclude_self = FALSE,
  return_mean = TRUE,
  verbose = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- nn1:

  Name of the first neighbor graph (e.g. "RNA.nn").

- nn2:

  Name of the second neighbor graph (e.g. "BCR.nn").

- exclude_self:

  If TRUE, remove the cell's own index from each neighbor list before
  counting shared neighbors.

- return_mean:

  If TRUE, return the mean shared neighbor count across all cells; else
  return a named numeric vector of per-cell counts.

- verbose:

  If TRUE, print a summary message.

## Value

A single numeric value if `return_mean = TRUE`, or a named numeric
vector of per-cell shared neighbor counts if `return_mean = FALSE`.
