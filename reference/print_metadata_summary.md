# Give an overview of the data in terms of files, samples and datasets.

This function provides a summary of metadata files, printing information
about the number of subjects, samples, and total files for each disease
in the dataset. Can handle metadata files representing multiple diseases
and datasets.

## Usage

``` r
print_metadata_summary(meta_file)
```

## Arguments

- meta_file:

  The analyst-created metadata csv file containing disease, sample, and
  subject information.

## Value

A text description of the provided metadata file.
