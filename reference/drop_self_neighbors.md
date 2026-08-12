# Drop self-neighbors from a k-NN index matrix, keeping k fixed for every cell

A cell can appear as its own nearest neighbor (e.g. a distance-0 self
match), and whether this happens can vary by cell or even be absent
entirely, depending on how the neighbor graph was built. To make the
output width predictable regardless, this *always* drops exactly one
neighbor per row: the self-entry if present, otherwise the single
farthest neighbor (the last column). Every row therefore always ends up
with `k - 1` neighbors, whether or not that cell (or any cell) had a
self-match.

Because of this, always build the neighbor graph with one extra neighbor
than you actually want (e.g. `k.param = 21` to end up with 20 real
neighbors per cell) — you don't need to know in advance whether
self-matches are present.

## Usage

``` r
drop_self_neighbors(nn_idx)
```

## Arguments

- nn_idx:

  A cell-by-k matrix of neighbor indices (`@nn.idx` from a Seurat
  neighbor graph), where row `i` gives the neighbors of cell `i`.

## Value

A cell-by-(k - 1) matrix of neighbor indices with no self-matches.
