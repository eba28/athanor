# Package index

## Metadata

Functions for handling metadata.

- [`print_metadata_summary()`](https://eba28.github.io/athanor/reference/print_metadata_summary.md)
  : Give an overview of the data in terms of files, samples and datasets

## Quality Control

Functions for QC and filtering.

- [`create_metrics_summary()`](https://eba28.github.io/athanor/reference/create_metrics_summary.md)
  : Create a table with info from 10x Genomics' metric summary file(s)
- [`layout_doublets()`](https://eba28.github.io/athanor/reference/layout_doublets.md)
  : Plot an overview of a doublet identification method
- [`plot_counts()`](https://eba28.github.io/athanor/reference/plot_counts.md)
  : Create bar plots of read or cell counts for quality control
- [`reformat_vdj_barcode()`](https://eba28.github.io/athanor/reference/reformat_vdj_barcode.md)
  : Reformat VDJ barcodes to make them unique across samples
- [`reformat_vdj_barcode_sample()`](https://eba28.github.io/athanor/reference/reformat_vdj_barcode_sample.md)
  : Reformat VDJ barcodes by adding sample names and removing suffixes

## Gene Expression

Functions for processing and analyzing gene expression data.

- [`add_annotations()`](https://eba28.github.io/athanor/reference/add_annotations.md)
  : Add a user-specified list of cluster annotations to a Seurat object

- [`automated_annotation()`](https://eba28.github.io/athanor/reference/automated_annotation.md)
  : Run automated cell type annotation

- [`cell_type_clusters()`](https://eba28.github.io/athanor/reference/cell_type_clusters.md)
  : Map cell types to Seurat clusters

- [`find_k_clusters()`](https://eba28.github.io/athanor/reference/find_k_clusters.md)
  : Find the right clustering resolution to obtain the desired number of
  clusters

- [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md)
  : Get IG and TR genes from Ensembl using biomaRt.

- [`gen_dot_title()`](https://eba28.github.io/athanor/reference/gen_dot_title.md)
  :

  Generate a title for a `DotPlot`

- [`get_cell_types()`](https://eba28.github.io/athanor/reference/get_cell_types.md)
  :

  Get the information needed for using
  [`add_info_bar()`](https://eba28.github.io/athanor/reference/add_info_bar.md)
  on a `DotPlot`

- [`get_features_from_all()`](https://eba28.github.io/athanor/reference/get_features_from_all.md)
  : Get specific markers from a marker genes database

- [`plot_counts_cluster()`](https://eba28.github.io/athanor/reference/plot_counts_cluster.md)
  : Visualize how many cells are in each Seurat cluster or cell type

- [`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
  : Run Seurat's standard pipeline

- [`source_markers()`](https://eba28.github.io/athanor/reference/source_markers.md)
  : Display markers from a filtered marker genes database as a table

## Adaptive Immune Receptors

Functions for processing and analyzing adaptive immune receptor data.

- [`add_family_info()`](https://eba28.github.io/athanor/reference/add_family_info.md)
  :

  Add in family and gene information from `alakazam`

- [`bin_mu_freq()`](https://eba28.github.io/athanor/reference/bin_mu_freq.md)
  : Bin the mutation frequency

- [`convert_embeddings()`](https://eba28.github.io/athanor/reference/convert_embeddings.md)
  : Convert the output of an embedding method to a matrix

- [`factor_family_info()`](https://eba28.github.io/athanor/reference/factor_family_info.md)
  : Convert family and gene information to sorted factors

- [`plot_immune_overlay()`](https://eba28.github.io/athanor/reference/plot_immune_overlay.md)
  : Plot UMAPs with AIRR overlays

- [`process_airrflow()`](https://eba28.github.io/athanor/reference/process_airrflow.md)
  : Read in and process output files from nf-core/airrflow

- [`process_bcr_features()`](https://eba28.github.io/athanor/reference/process_bcr_features.md)
  : Process BCR features for integration

## GEX and AIRR Integration

Functions for integrating gene expression and adaptive immune receptor
data.

### Metadata

- [`gex_add_airr()`](https://eba28.github.io/athanor/reference/gex_add_airr.md)
  : Add AIRR information to a Seurat object

### Concatenation

- [`concatenate_gex_bcr()`](https://eba28.github.io/athanor/reference/concatenate_gex_bcr.md)
  : Concatenate GEX and BCR data in a Seurat object

### Weighted Nearest Neighbors

- [`extract_wnn_vars()`](https://eba28.github.io/athanor/reference/extract_wnn_vars.md)
  : Give a summary of a Seurat object post-WNN
- [`plot_mws()`](https://eba28.github.io/athanor/reference/plot_mws.md)
  : Plot a box plot of modality weights per cell type
- [`plot_wnn_testing()`](https://eba28.github.io/athanor/reference/plot_wnn_testing.md)
  : Plot the results of manual WNN testing
- [`plot_wnn_umaps()`](https://eba28.github.io/athanor/reference/plot_wnn_umaps.md)
  : Plot UMAPs of a Seurat object post-WNN
- [`run_wnn()`](https://eba28.github.io/athanor/reference/run_wnn.md) :
  Run Weighted Nearest Neighbors (WNN) analysis on combined GEX and BCR
  data

## Evaluation

Functions for evaluating the performance of models and analyses.

- [`calc_adt_dists()`](https://eba28.github.io/athanor/reference/calc_adt_dists.md)
  : Compute mean ADT distance to each cell's k nearest neighbors
- [`calc_adt_dists_fast()`](https://eba28.github.io/athanor/reference/calc_adt_dists_fast.md)
  : Compute mean ADT distance to each cell's k nearest neighbors (faster
  version)
- [`calc_adt_nn_within_range()`](https://eba28.github.io/athanor/reference/calc_adt_nn_within_range.md)
  : Calculate the proportion of neighbors within an ADT marker's
  expression range
- [`calc_adt_quantile()`](https://eba28.github.io/athanor/reference/calc_adt_quantile.md)
  : Calculate the proportion of neighbors within an ADT marker's
  quantile by expression
- [`calc_adt_scores()`](https://eba28.github.io/athanor/reference/calc_adt_scores.md)
  : Calculate homogeneity scores for binary ADT features across
  embeddings
- [`calc_distances()`](https://eba28.github.io/athanor/reference/calc_distances.md)
  : Calculate cluster distances for a Seurat object
- [`calc_ext_metrics()`](https://eba28.github.io/athanor/reference/calc_ext_metrics.md)
  : Calculate external clustering metrics for a Seurat object
- [`calc_int_metrics()`](https://eba28.github.io/athanor/reference/calc_int_metrics.md)
  : Calculate internal clustering metrics for a Seurat object
- [`calc_moran()`](https://eba28.github.io/athanor/reference/calc_moran.md)
  : Calculate Moran's i for a Seurat object
- [`calc_neighbor_matches()`](https://eba28.github.io/athanor/reference/calc_neighbor_matches.md)
  : Calculate neighbor matching scores across metadata columns
- [`plot_metrics()`](https://eba28.github.io/athanor/reference/plot_metrics.md)
  : Plot the internal or external clustering metrics

## Applications

Functions for applying the package to specific use cases.

### Reclassification

- [`calc_nn_frac()`](https://eba28.github.io/athanor/reference/calc_nn_frac.md)
  : Calculate the proportion of a cell's neighbors that are mutated or
  unmutated

## Plotting

Functions for visualizing data and results.

- [`add_info_bar()`](https://eba28.github.io/athanor/reference/add_info_bar.md)
  : Add an information bar on top of a given ggplot

- [`calc_pcts()`](https://eba28.github.io/athanor/reference/calc_pcts.md)
  : Calculate percentages per metadata group

- [`plot_color_scale()`](https://eba28.github.io/athanor/reference/plot_color_scale.md)
  :

  Generate a color scale for a Seurat `DotPlot` with white at zero

- [`plot_counts()`](https://eba28.github.io/athanor/reference/plot_counts.md)
  : Create bar plots of read or cell counts for quality control

- [`plot_counts_cluster()`](https://eba28.github.io/athanor/reference/plot_counts_cluster.md)
  : Visualize how many cells are in each Seurat cluster or cell type

- [`plot_dimplot()`](https://eba28.github.io/athanor/reference/plot_dimplot.md)
  :

  Plot a Seurat UMAP using `DimPlot`

- [`plot_dot_airr()`](https://eba28.github.io/athanor/reference/plot_dot_airr.md)
  :

  Add AIRR (and other) info along the right side of an existing Seurat
  `DotPlot`

- [`plot_immune_overlay()`](https://eba28.github.io/athanor/reference/plot_immune_overlay.md)
  : Plot UMAPs with AIRR overlays

- [`plot_metrics()`](https://eba28.github.io/athanor/reference/plot_metrics.md)
  : Plot the internal or external clustering metrics

- [`plot_mws()`](https://eba28.github.io/athanor/reference/plot_mws.md)
  : Plot a box plot of modality weights per cell type

- [`plot_overview_comps()`](https://eba28.github.io/athanor/reference/plot_overview_comps.md)
  : Plot several UMAPs side by side for a Seurat object

- [`plot_pcts()`](https://eba28.github.io/athanor/reference/plot_pcts.md)
  : Plot percentages in a stacked bar plot

- [`plot_umap()`](https://eba28.github.io/athanor/reference/plot_umap.md)
  : Plot Seurat UMAP(s) in several useful ways

- [`plot_umap_condition()`](https://eba28.github.io/athanor/reference/plot_umap_condition.md)
  : Plot a specific condition on a Seurat UMAP

- [`plot_wnn_testing()`](https://eba28.github.io/athanor/reference/plot_wnn_testing.md)
  : Plot the results of manual WNN testing

- [`plot_wnn_umaps()`](https://eba28.github.io/athanor/reference/plot_wnn_umaps.md)
  : Plot UMAPs of a Seurat object post-WNN

- [`vln_feat_plot()`](https://eba28.github.io/athanor/reference/vln_feat_plot.md)
  :

  Plot a Seurat `VlnPlot` and a `FeaturePlot` side by side for the same
  marker

## Utilities

Utility functions for various tasks.

- [`print_dt()`](https://eba28.github.io/athanor/reference/print_dt.md)
  :

  Display a sortable, scrollable `DataTable` table

- [`print_kable()`](https://eba28.github.io/athanor/reference/print_kable.md)
  : Display a nicely formatted table

- [`reduce_object()`](https://eba28.github.io/athanor/reference/reduce_object.md)
  : Reduce a Seurat object's size

## Data

Simulate data.

- [`run_wnn_sims()`](https://eba28.github.io/athanor/reference/run_wnn_sims.md)
  : Run manual WNN simulations by varying a specific variable while
  keeping others constant
- [`sim_bcr_manual()`](https://eba28.github.io/athanor/reference/sim_bcr_manual.md)
  : Simulate an assay made from BCR embeddings independent of any GEX
  object
- [`sim_gex_manual()`](https://eba28.github.io/athanor/reference/sim_gex_manual.md)
  : Manually simulate gene expression data
- [`sim_gex_splatter()`](https://eba28.github.io/athanor/reference/sim_gex_splatter.md)
  : Simulate gene expression data using Splatter
