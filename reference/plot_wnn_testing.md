# Plot the results of the manual WNN testing.

Plot the results of the manual WNN testing.

## Usage

``` r
plot_wnn_testing(manual_wnn_test, sim_var, other_vars, count_range)
```

## Arguments

- manual_wnn_test:

  A data frame containing the results of the manual WNN testing, with
  columns "Count" and "Passed".

- sim_var:

  The variable that was varied in the manual WNN testing (e.g., "Genes",
  "Cells", "Dimensions", "GEX PCs", or "BCR PCs").

- other_vars:

  A named list of the other variables that were held constant in the
  manual WNN testing, with names corresponding to the variable names
  (e.g., "Genes", "Cells", "Dimensions", "GEX PCs", or "BCR PCs") and
  values corresponding to the constant values used in the testing.

- count_range:

  A numeric vector specifying the range of counts that were tested in
  the manual WNN testing, used for setting the x-axis breaks in the
  plot.

## Value

A ggplot object showing the step plot of the WNN testing results, with
the x-axis representing the count of the varied variable and the y-axis
representing whether Seurat's WNN computation passed (1), failed (0), or
gave a warning (0.5). The plot includes a title indicating the variable
that was varied and a subtitle listing the other variables that were
held constant.
