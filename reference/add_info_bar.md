# Add an information bar on top of a given plot.

This function adds informational strips (facets) to the top of an
existing `ggplot2`.

## Usage

``` r
add_info_bar(
  plot,
  method = "contains",
  info_type,
  info,
  sort = TRUE,
  text_size = 8,
  label = TRUE,
  angle = 0,
  top_side = TRUE
)
```

## Arguments

- plot:

  A ggplot2 object to which information bars will be added.

- method:

  The method for getting the info bar information; one of `add`,
  `contains`, `join`.

- info_type:

  What to add e.g. Dataset, Cell_Type, etc.

- info:

  Data frame containing the information to join/add, must contain the
  same column as plot\$data to join by.

- sort:

  Sort the info alphanumerically.

- text_size:

  The size of the text in the bars.

- label:

  Whether or not to include what you are adding info for (good for
  datasets).

- angle:

  Rotation degree.

- top_side:

  Plot on top (by default) or on the right side.

## Value

A ggplot plot.

## Details

This is very useful for splitting by dataset, showing cell types, etc.
`contains` is not really needed. Note that providing a named list of
features to a DotPlot will automatically group them without the need for
this function.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- ggplot(data, aes(x, y)) + geom_point()
add_info_bar(p, info_type = "Dataset", info = dataset_info)
} # }
```
