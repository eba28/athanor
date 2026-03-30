# Calculate Moran's i for a Seurat object

This function calculates the global Moran's i index.

## Usage

``` r
calc_moran(seurat_obj, feature, graph_name, row_standardize = TRUE)
```

## Arguments

- seurat_obj:

  The Seurat object. Must have `FindNeighbors()` already run and an
  assay named "ADT".

- feature:

  Name of the ADT feature to evaluate (e.g. "CD27.1", "CD38").

- graph_name:

  Name of the neighbor graph slot to use for the weights matrix (e.g.
  "RNA.nn", "BCR.nn", "w.nn").

- row_standardize:

  Whether or not to row-standardize the weights matrix (i.e. make each
  row sum to 1).

## Value

A single numeric value representing the observed Moran's i index for the
specified feature and neighbor graph.

## Details

We are using `MERINGUE`'s implementation instead of `ape`'s because it
runs faster. However, `MERINGUE` is not on CRAN, which means this
package could not be published on CRAN. Row standardization makes sure
that the resulting score will always be between -1 and 1.
