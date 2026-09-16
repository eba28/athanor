# Calculate Moran's i

Calculates the global Moran's i index for an ADT feature given a
neighbor graph.

## Usage

``` r
calc_adt_moran(
  seurat_obj,
  features_adt,
  adt_assay = "ADT",
  graph_name,
  row_standardize = TRUE,
  verbose = FALSE
)
```

## Arguments

- seurat_obj:

  The Seurat object. Must have `FindNeighbors()` already run.

- features_adt:

  Name of the ADT feature to evaluate (e.g. "CD27.1").

- adt_assay:

  Name of the assay containing ADT data.

- graph_name:

  Name of the neighbor graph slot to use for the weights matrix (e.g.
  "RNA.nn", "BCR.nn", "w.nn").

- row_standardize:

  Whether to row-standardize the weights matrix (i.e. make each row sum
  to 1).

- verbose:

  If TRUE, print progress messages.

## Value

A one-row data frame with columns: Graph, Assay, Feature, Method, Score
(the canonical schema shared by the `calc_*` evaluation functions),
where `Score` is the observed Moran's i and `Method` is "Moran_I".

## Details

Row standardization makes sure the resulting score will always be
between -1 and 1.
