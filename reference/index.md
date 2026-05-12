# Package index

## General

General functions for processing and analyzing single-cell data.

- [`object_overview()`](https://eba28.github.io/athanor/reference/object_overview.md)
  : Print an overview of a Seurat object.
- [`regen_reduc()`](https://eba28.github.io/athanor/reference/regen_reduc.md)
  : Regenerate neighbor graphs and UMAPs.
- [`seurat_pipeline()`](https://eba28.github.io/athanor/reference/seurat_pipeline.md)
  : Run Seurat's standard pipeline

## Quality Control

Functions for QC and filtering.

- [`reformat_vdj_barcode()`](https://eba28.github.io/athanor/reference/reformat_vdj_barcode.md)
  : Reformat VDJ barcodes to make them unique across samples
- [`reformat_vdj_barcode_sample()`](https://eba28.github.io/athanor/reference/reformat_vdj_barcode_sample.md)
  : Reformat VDJ barcodes by adding sample names and removing suffixes

## Gene Expression

Functions for processing and analyzing gene expression data.

### Pipeline

- [`filter_variable_features()`](https://eba28.github.io/athanor/reference/filter_variable_features.md)
  : Filter AIRR genes from variable features
- [`get_airr_genes()`](https://eba28.github.io/athanor/reference/get_airr_genes.md)
  : Get IG and TR genes from Ensembl using biomaRt

### Cell type annotation

- [`add_annotations()`](https://eba28.github.io/athanor/reference/add_annotations.md)
  : Add a user-specified list of cluster annotations to a Seurat object
- [`automated_annotation()`](https://eba28.github.io/athanor/reference/automated_annotation.md)
  : Run automated cell type annotation

### Clustering-related

- [`cell_type_clusters()`](https://eba28.github.io/athanor/reference/cell_type_clusters.md)
  : Map cell types to Seurat clusters
- [`find_k_clusters()`](https://eba28.github.io/athanor/reference/find_k_clusters.md)
  : Find the right clustering resolution to obtain the desired number of
  clusters

## Adaptive Immune Receptors

Functions for processing and analyzing adaptive immune receptor data.

### Read in data

- [`process_airrflow()`](https://eba28.github.io/athanor/reference/process_airrflow.md)
  : Read in and process output files from nf-core/airrflow

### BCR embeddings

- [`bcr_embeddings_pipeline()`](https://eba28.github.io/athanor/reference/bcr_embeddings_pipeline.md)
  : Build a Seurat object from BCR embeddings
- [`bcr_embeddings_pipeline_dedup()`](https://eba28.github.io/athanor/reference/bcr_embeddings_pipeline_dedup.md)
  : Build a Seurat object from BCR embeddings, deduplicating identical
  embeddings first
- [`convert_embeddings()`](https://eba28.github.io/athanor/reference/convert_embeddings.md)
  : Convert the output of an embedding method to a matrix

### BCR features

- [`add_family_info()`](https://eba28.github.io/athanor/reference/add_family_info.md)
  :

  Add in family and gene information from `alakazam`

- [`bin_mu_freq()`](https://eba28.github.io/athanor/reference/bin_mu_freq.md)
  : Bin the mutation frequency

- [`factor_family_info()`](https://eba28.github.io/athanor/reference/factor_family_info.md)
  : Convert family and gene information to sorted factors

- [`process_bcr_features()`](https://eba28.github.io/athanor/reference/process_bcr_features.md)
  : Process BCR features for integration

## GEX and AIRR Integration

Functions for integrating gene expression and adaptive immune receptor
data.

### Metadata

- [`gex_add_airr()`](https://eba28.github.io/athanor/reference/gex_add_airr.md)
  : Add AIRR information to a Seurat object

### Seurat objects

- [`infer_k_param()`](https://eba28.github.io/athanor/reference/infer_k_param.md)
  : Infer k_param from existing RNA neighbors, or return a default.
- [`infer_num_dims()`](https://eba28.github.io/athanor/reference/infer_num_dims.md)
  : Infer num_dims (PCA dims used for neighbor finding) from existing
  RNA neighbors.
- [`merge_gex_bcr()`](https://eba28.github.io/athanor/reference/merge_gex_bcr.md)
  : Merge a GEX Seurat object with a BCR Seurat object

### Concatenation

- [`concatenate_gex_bcr()`](https://eba28.github.io/athanor/reference/concatenate_gex_bcr.md)
  : Concatenate GEX and BCR data in a Seurat object

### Weighted Nearest Neighbors

- [`run_wnn()`](https://eba28.github.io/athanor/reference/run_wnn.md) :
  Run Weighted Nearest Neighbors (WNN) analysis on combined GEX and BCR
  data

## Evaluation

Functions for evaluating the performance of models and analyses.

### ADT-based on k-nearest neighbors

- [`calc_adt_correlation()`](https://eba28.github.io/athanor/reference/calc_adt_correlation.md)
  : Calculate the correlation between each cell's expression and the
  mean of its neighbors' expression
- [`calc_adt_dists()`](https://eba28.github.io/athanor/reference/calc_adt_dists.md)
  : Compute mean ADT distance to each cell's k nearest neighbors
- [`calc_adt_moran()`](https://eba28.github.io/athanor/reference/calc_adt_moran.md)
  : Calculate Moran's i
- [`calc_adt_nn_within_range()`](https://eba28.github.io/athanor/reference/calc_adt_nn_within_range.md)
  : Calculate the proportion of neighbors within an ADT marker's
  expression range
- [`calc_adt_quantile()`](https://eba28.github.io/athanor/reference/calc_adt_quantile.md)
  : Calculate the proportion of neighbors within an ADT marker's
  quantile by expression
- [`calc_neighbor_matches()`](https://eba28.github.io/athanor/reference/calc_neighbor_matches.md)
  : Calculate neighbor matching scores across metadata columns

### Clustering-based

- [`calc_distances()`](https://eba28.github.io/athanor/reference/calc_distances.md)
  : Calculate cluster distances
- [`calc_ext_metrics()`](https://eba28.github.io/athanor/reference/calc_ext_metrics.md)
  : Calculate external clustering metrics
- [`calc_int_metrics()`](https://eba28.github.io/athanor/reference/calc_int_metrics.md)
  : Calculate internal clustering metrics

## Applications

Functions for applying the package to specific use cases.

### Reclassification

## Plotting

Functions for visualizing data and results.

### Accessory

- [`add_info_bar()`](https://eba28.github.io/athanor/reference/add_info_bar.md)
  : Add an information bar on top of a given ggplot
- [`plot_color_scale()`](https://eba28.github.io/athanor/reference/plot_color_scale.md)
  : Generate a color scale for a ggplot2 plot with white at zero

### UMAPs

- [`plot_dimplot()`](https://eba28.github.io/athanor/reference/plot_dimplot.md)
  :

  Plot a Seurat UMAP using `DimPlot`

- [`plot_overview_comps()`](https://eba28.github.io/athanor/reference/plot_overview_comps.md)
  : Plot several UMAPs side by side

### Other plots

- [`calc_pcts()`](https://eba28.github.io/athanor/reference/calc_pcts.md)
  : Calculate percentages per metadata group

- [`plot_doublets()`](https://eba28.github.io/athanor/reference/plot_doublets.md)
  : Plot an overview of a doublet identification method

- [`plot_mws()`](https://eba28.github.io/athanor/reference/plot_mws.md)
  : Plot a box plot of modality weights per cell type

- [`plot_pcts()`](https://eba28.github.io/athanor/reference/plot_pcts.md)
  : Plot percentages in a stacked bar plot

- [`plot_vln_feat()`](https://eba28.github.io/athanor/reference/plot_vln_feat.md)
  :

  Plot a Seurat `VlnPlot` and a `FeaturePlot` side by side for the same
  marker

## Simulations

Simulate GEX and AIRR data.

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

## Utilities

Utilities for various tasks.

### Colors

- [`named_colors`](https://eba28.github.io/athanor/reference/named_colors.md)
  : Named color palettes for athanor plots

### Seurat object

- [`reduce_object()`](https://eba28.github.io/athanor/reference/reduce_object.md)
  : Reduce a Seurat object's size

### Themes

- [`clean_dimplot`](https://eba28.github.io/athanor/reference/clean_dimplot.md)
  : Clean theme for DimPlots
- [`clean_dimplot2`](https://eba28.github.io/athanor/reference/clean_dimplot2.md)
  : Clean theme for DimPlots
- [`labels_rotate_x`](https://eba28.github.io/athanor/reference/labels_rotate_x.md)
  : Rotate x-axis labels 45 degrees
- [`labels_standard`](https://eba28.github.io/athanor/reference/labels_standard.md)
  : Standard label sizes for athanor plots
- [`labels_standard_vln`](https://eba28.github.io/athanor/reference/labels_standard_vln.md)
  : Standard label sizes for violin plots
- [`labels_standard_vln_rotate`](https://eba28.github.io/athanor/reference/labels_standard_vln_rotate.md)
  : Violin plot labels with horizontal x-axis text
- [`plot_anno`](https://eba28.github.io/athanor/reference/plot_anno.md)
  : Patchwork annotation with Roman numeral panel labels
- [`theme_bw_custom`](https://eba28.github.io/athanor/reference/theme_bw_custom.md)
  : ggplot2 theme with horizontal grid lines only
