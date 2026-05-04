# Named color palettes for athanor plots

A list of named character vectors mapping categorical values to
hexadecimal colors, typically used as the `clrs_specific` argument in
plotting functions. These color mappings are used to ensure consistent
and meaningful coloring across different plots in the athanor package,
especially for categories that are commonly visualized (e.g. cell types,
isotypes, SHM frequencies).

## Usage

``` r
named_colors
```

## Details

The `named_colors` list contains color mappings for various categories
relevant to BCR and single-cell analysis, including:

- `c_call`: Colors for IG constant region calls (e.g. IGHA, IGHD, etc.)

- `cdr3`: Colors for CDR3 lengths from 4 to 41

- `cell_types_celltypist`, `cell_types_simpler`: Colors for cell type
  annotations from CellTypist and a simpler scheme

- `datatype`: Colors for different data types (e.g. ADT, BCR, GEX, etc.)

- `doublet`: Colors for singlets vs doublets in doublet detection

- `embeddings`: Colors for different embedding methods (e.g. antiberta2,
  balm-paired, etc.)

- `isotype_stage`: Colors for isotype switching stages (i.e. Unswitched,
  Switched)

- `isotype`: Colors for immunoglobulin isotypes (e.g. IgA, IgD, etc.)

- `light`: Colors for light chain types (i.e. IGK, IGL, IGK/IGL)

- `mu_freq_bins`, `mu_freq_bins_fewer`, `mu_freq_bins_binary`: Colors
  for binned mutation frequencies

- `mu_freq_iso`: Colors for combinations of mutation frequency bins and
  isotype switching stages

- `v_call_family`, `d_call_family`, `j_call_family`: Colors for V, D, J
  gene families respectively

- `weights`: Colors for WNN weights from 0 to 1

When possible, colors are chosen to be colorblind-friendly.
