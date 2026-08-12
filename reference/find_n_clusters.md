# Find the right clustering resolution to obtain the desired number of clusters

This function iteratively tests clustering resolutions in Seurat to find
the resolution that yields the desired number of clusters. It uses the
specified graph and returns the Seurat object with clusters if
successful, or stops if the desired number is exceeded or not found.

## Usage

``` r
find_n_clusters(
  seurat_obj,
  graph_name = "RNA_snn",
  desired_n,
  res_range = seq(0.1, 2, by = 0.1)
)
```

## Arguments

- seurat_obj:

  The Seurat object.

- graph_name:

  The name of the graph to use for clustering.

- desired_n:

  The desired number of clusters.

- res_range:

  The range of resolutions to test.

## Value

The Seurat object with clusters at the resolution that matches
`desired_n`.

## Details

This would typically be used after
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md).
