# Find the right clustering resolution to obtain the desired number of clusters

This function iteratively tests clustering resolutions in Seurat to find
the resolution that yields the desired number of clusters. It uses the
specified graph and returns the Seurat object with clusters if
successful, or stops if the desired number is exceeded or not found.

## Usage

``` r
find_k_clusters(seurat_obj, graph_name = "RNA_snn", desired_k)
```

## Arguments

- seurat_obj:

  The Seurat object.

- graph_name:

  The name of the graph to use for clustering.

- desired_k:

  The desired number of clusters.

## Value

The Seurat object with clusters at the resolution that matches
desired_k.
