# Convert from the output of an embeddings method to a `Matrix` matrix.

Convert from the output of an embeddings method to a `Matrix` matrix.

## Usage

``` r
convert_embeddings(embeddings, combined_airr, combined_airr_input)
```

## Arguments

- embeddings:

  The data.frame embeddings output.

- combined_airr:

  The combined output from airrflow/Immcantation. Must contain columns
  called "cell_id" and "cell_id_original".

- combined_airr_input:

  The data.frame provided to immune2vec; contains translated sequences.

## Value

A `Matrix` of embeddings

## Details

Assume that all inputs can be provided as tsv files. The AMULETy outputs
always have a column named "cell_id". immune2vec does not include the
cell ids in the output. Uses `Matrix()` instead of
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) since the latter
only returns dense matrices. AntiBERTy runs on 512 dimensions,
AntiBERTa2 and BALM-paired run on 1024 dimensions, and ESM2 runs on 1280
dimensions (through AMULETy). immune2vec is usually run with 100
dimensions.
