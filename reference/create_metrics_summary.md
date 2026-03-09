# Create a table with info from 10x's metric summary file(s)

This function reads and combines metrics summary files from 10x Genomics
Cell Ranger outputs, extracting key quality control metrics like
estimated number of cells, number of reads, and reads per cell for
multiple samples and data types.

## Usage

``` r
create_metrics_summary(meta, path_data_specific, data_types, samples_list)
```

## Arguments

- meta:

  The metadata file containing sample and dataset information.

- path_data_specific:

  Where the data is located (path to the data directory).

- data_types:

  Vector of data types to process (e.g. "BCR", "GEX", "TCR"). If
  missing, uses all data types found in meta.

- samples_list:

  Vector of sample names to examine. If missing, uses all samples found
  in meta.

## Value

A tibble with columns: sample_id, subject_id, SampleType, DataType,
Dataset, EstimatedNumberofCells, NumberofReads, ReadsPerCell.

## Details

You can provide a specific list of samples or data types to create
summaries for if you don't want to use everything in the given metadata
file.
