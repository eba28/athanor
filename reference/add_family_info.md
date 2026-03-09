# Add in family and gene information from `alakazam.`

Extracts and adds V, D, and J gene family and gene information to an
AIRR-formatted dataframe using the alakazam package functions.

## Usage

``` r
add_family_info(combined_airr)
```

## Arguments

- combined_airr:

  An AIRR-formatted data.frame.

## Value

A data.frame with six new columns containing gene family and gene
information.

## Details

biomaRt also has a getGene function, so we have to be specific.
