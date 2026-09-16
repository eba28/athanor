# Calculate shared neighbors between two neighbor graphs

For each cell in a Seurat object, counts how many of its nearest
neighbors are shared between two specified neighbor graphs.

## Usage

``` r
calc_neighbor_overlap(
  seurat_obj,
  nn1,
  nn2,
  drop_self = TRUE,
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

- drop_self:

  Logical indicating whether to drop one neighbor per cell (the
  self-match if present, else the farthest neighbor) from each graph so
  `k` is fixed across cells regardless of whether self-matches are
  present. See
  [`drop_self_neighbors()`](https://eba28.github.io/athanor/reference/drop_self_neighbors.md)
  for details.

- return_mean:

  If TRUE, return the mean shared neighbor count across all cells; else
  return a named numeric vector of per-cell counts.

- verbose:

  If TRUE, print a summary message.

## Value

Data frame with columns: Graph, Assay, Feature, Method, Score (the
canonical schema shared by the `calc_*` evaluation functions). `Graph`
is `"{nn1} vs {nn2}"`, `Assay` and `Feature` are NA (no ADT feature is
involved), and `Method` is "Overlap". If `return_mean = TRUE`, `Score`
is the mean shared-neighbor count across all cells; else it contains
per-cell counts with an additional `cell_id` column.
