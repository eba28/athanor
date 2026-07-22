# Run manual WNN simulations by varying a specific variable while keeping others constant

This function runs a series of WNN simulations by varying a specified
variable (e.g. number of genes, cells, dimensions, GEX PCs, BCR PCs)
while keeping other parameters constant. It captures whether each
simulation run passes without errors or warnings and returns a data
frame summarizing the results.

## Usage

``` r
run_wnn_sims(count_range, sim_var, other_vars, show_progress = FALSE)
```

## Arguments

- count_range:

  A vector of values for the variable to be varied in the simulations.

- sim_var:

  The name of the variable being varied (e.g. "Genes", "Cells",
  "Dimensions", "GEX PCs", "BCR PCs").

- other_vars:

  A list of other variables and their constant values to be used in the
  simulations.

- show_progress:

  Whether to display a progress bar during the simulations.

## Value

A data frame with columns "Count" and "Passed", where "Count"
corresponds to the values in `count_range` and "Passed" indicates
whether the simulation run passed without errors (1), had warnings
(0.5), or had errors (0).

## Details

It uses both
[`sim_airr_manual()`](https://eba28.github.io/athanor/reference/sim_airr_manual.md)
and
[`sim_gex_manual()`](https://eba28.github.io/athanor/reference/sim_gex_manual.md)
to generate the necessary data for the WNN simulations. The function
also includes error handling to capture any issues that arise during the
simulation runs and provides an option to display a progress bar.
