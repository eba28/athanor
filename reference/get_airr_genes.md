# Get IG and TR genes from Ensembl using biomaRt

This function retrieves IG and TR genes from Ensembl using the biomaRt
package. This allows for accurate IG and TR gene names instead of using
a search for genes that begin with "IG" or "TR".

## Usage

``` r
get_airr_genes(
  genome = "hsapiens",
  ensembl_version = NULL,
  category = c("IG", "TR"),
  cache_file = NULL
)
```

## Arguments

- genome:

  The genome to use for gene annotation (e.g. "hsapiens" or
  "mmusculus").

- ensembl_version:

  The Ensembl version to use for gene annotation (e.g. "114").

- category:

  The category of genes to retrieve: "IG" for immunoglobulin genes, "TR"
  for T cell receptor genes, or both.

- cache_file:

  Optional path to an RDS file. If the file exists, its contents are
  returned directly without querying Ensembl. After a successful Ensembl
  query the result is saved to this path for future offline use.

## Value

A character vector of IG and/or TR gene names to be filtered out from
the most variable features.

## Details

This will usually just be used as part of
[`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md).
