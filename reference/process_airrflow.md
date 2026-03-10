# Read in and process output files from nf-core/airrflow

This function reads airrflow's repertoire files, processes them to add
subject and sample information, computes CDR3 amino acid properties, and
adds gene family information.

## Usage

``` r
process_airrflow(dataset_path, version_airrflow)
```

## Arguments

- dataset_path:

  The path to the dataset directory.

- version_airrflow:

  The airrflow version (as a string).

## Value

A processed AIRR-formatted data.frame with several columns added.

## Details

Only written for BCR data right now. v4.0 has `clone_size_count` and
`clone_size_freq` columns. v4.3.1 has `duplicate_count` and
`light_only_cell` columns. There can still be some `NA` `c_call`s.
