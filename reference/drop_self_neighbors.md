# Drop self-neighbors from a k-NN index matrix, keeping k fixed for every cell

A cell can appear as its own nearest neighbor (e.g. a distance-0 self
match). This drops that self-entry from each row so it doesn't trivially
"match" itself, while keeping the neighbor count the same for every
cell: rows without a self-match instead drop their single farthest
neighbor (the last column), so all rows end up with `k - 1` neighbors.

## Usage

``` r
drop_self_neighbors(nn_idx)
```

## Arguments

- nn_idx:

  A cell-by-k matrix of neighbor indices (`@nn.idx` from a Seurat
  neighbor graph), where row `i` gives the neighbors of cell `i`.

## Value

A cell-by-(k - 1) matrix of neighbor indices with no self-matches. If no
cell is its own neighbor, `nn_idx` is returned unchanged.
