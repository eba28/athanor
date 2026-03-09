#' Create a table with info from 10x's metric summary file(s)
#'
#' @description
#' This function reads and combines metrics summary files from 10x Genomics Cell Ranger
#' outputs, extracting key quality control metrics like estimated number of cells,
#' number of reads, and reads per cell for multiple samples and data types.
#'
#' @details
#' You can provide a specific list of samples or data types to create summaries for if you don't want to use everything in the given metadata file.
#'
#'
#' @param meta The metadata file containing sample and dataset information.
#' @param path_data_specific Where the data is located (path to the data directory).
#' @param data_types Vector of data types to process (e.g. "BCR", "GEX", "TCR").
#'   If missing, uses all data types found in meta.
#' @param samples_list Vector of sample names to examine. If missing, uses all
#'   samples found in meta.
#'
#' @returns A tibble with columns: sample_id, subject_id, SampleType, DataType, Dataset,
#'   EstimatedNumberofCells, NumberofReads, ReadsPerCell.
#' @export
create_metrics_summary <- function(meta, path_data_specific,
                                   data_types, samples_list) {
  metrics_summary <- tibble()
  if (rlang::is_missing(data_types)) data_types <- unique(meta$DataType) # depends on dataset
  if (rlang::is_missing(samples_list)) samples_list <- sort(unique(meta$sample_id))

  for (sample in samples_list) {
    for (data_type in data_types) {
      # metadata info
      meta_sample <- filter(meta, sample_id == sample, DataType == data_type)
      sample_type <- unique(meta_sample$SampleType)
      subject <- unique(meta_sample$subject_id)
      dataset <- unique(meta_sample$Dataset)

      # read in the summary table
      path_data <- file.path(path_data_specific, paste0("dataset", dataset),
                             "default")
      data <- read_csv(file.path(path_data, meta_sample$SampleDir,
                                 "outs", "metrics_summary.csv"),
                       show_col_types = FALSE)

      # set the columns of interest
      estimated_number_of_cells <- data[["Estimated Number of Cells"]]
      num_reads <- data[grepl("Number of Read", colnames(data))][[1]] # Number of Read Pairs (VDJ), Number of Reads (GEX)
      reads_per_cell <- num_reads %/% estimated_number_of_cells # mean read pairs per cell is also a col

      # update the combined summary table
      metrics_summary <- bind_rows(metrics_summary,
                                   tibble(sample, subject, sample_type, data_type,
                                          dataset, estimated_number_of_cells,
                                          num_reads, reads_per_cell))
    }
  }

  colnames(metrics_summary) <- c("SampleName", "Subject", "SampleType",
                                 "DataType", "Dataset",
                                 "EstimatedNumberofCells", "NumberofReads",
                                 "ReadsPerCell")
  return(metrics_summary)
}


#' Reformat VDJ barcodes by adding sample names and removing suffixes
#'
#' @description
#' This function adds the sample name to the beginning of barcodes and removes
#' the "-1" suffixes from the ends to create unique cell identifiers.
#'
#' @param data The input dataset which contains a barcode column.
#' @param sample_name The sample name to add to make the barcode unique.
#' @param barcode_col The column in the data that contains the barcodes.
#'
#' @returns A character vector of reformatted barcodes.
#' @export
reformat_vdj_barcode_sample <- function(data, sample_name,
                                        barcode_col = "barcode") {
  unlist(lapply(strsplit(as.character(data[[barcode_col]]), "-"),
                function(x) paste0(sample_name, "_", x[1])))
}

#' Reformat VDJ barcodes to make them unique across samples
#'
#' @description
#' This function processes barcode data by separating cell IDs from suffixes and
#' creating unique cell identifiers by combining sample names with barcodes.
#'
#' @param data Input data frame containing barcode information.
#' @param col_samples Column name containing sample names.
#' @param col_barcodes Column name containing cell IDs/barcodes.
#' @param col_output Name for the output column.
#'
#' @returns A character vector of unique cell identifiers.
#' @export
reformat_vdj_barcode <- function(data, col_samples = "sample_id",
                                 col_barcodes = "cell_id",
                                 col_output = "cell_id") {
  data <- data %>%
    select(all_of(col_samples), all_of(col_barcodes)) %>%
    separate(!!sym(col_barcodes), sep = "-",
             into = "barcodes", remove = FALSE, extra = "drop") %>%
    unite(col_output, c(!!sym(col_samples), "barcodes"), sep = "_",
          remove = FALSE)

  return(data$col_output)
}


#' Create bar plots of read or cell counts for quality control
#'
#' @description
#' This function generates bar plots showing read counts, cell counts, or barcode
#' counts for quality control purposes. Can display data split by data type,
#' sample, or other grouping variables.
#'
#' @param summary_df Output from create_metrics_summary().
#' @param data_types Vector of data types to include (e.g. "BCR", "GEX", "TCR").
#'   Use "All" to include all data types.
#' @param count_type Type of count to plot: "Read", "Cell", or "Barcode".
#' @param aggregation_state Whether data is "Aggregated" or "Unaggregated".
#'   (from cellranger aggr).
#' @param fill_type_label Label for the fill aesthetic.
#' @param x_axis Variable to plot along the x-axis.
#' @param x_axis_label Label for the x-axis.
#' @param clrs_datatype Named vector of colors for data types. Must be given in
#'   the correct order if plotting barcodes.
#'
#' @returns A ggplot2 barplot.
#' @export
plot_counts <- function(summary_df, data_types = "All", count_type,
                        aggregation_state = "Unaggregated",
                        fill_type_label = "Data Type",
                        x_axis = "sample_id", x_axis_label = "Sample",
                        clrs_datatype) {
  # data details
  if ("All" %in% data_types) {
    data_types <- unique(summary_df$DataType)
    data_type_desc <- "All"
  } else {
    data_type_desc <- paste(data_types, collapse = "/")
  }

  # setup the aesthetics
  if (count_type != "Barcode") {
    if (count_type == "Read") y_axis_label <- "Number of Reads"
    else if (count_type == "reads_per_cell") y_axis_label <- "Reads Per Cell"
    else y_axis_label <- "Estimated Number of Cells"
    y_axis <- str_replace_all(y_axis_label, " ", "") # remove spaces

    summary_df <- filter(summary_df, DataType %in% data_types)
    data_type_desc <- paste(data_type_desc, aggregation_state)
  } else {
    y_axis <- "Count"
    y_axis_label <- "Barcode Count"
  }

  fill_type <- str_replace_all(fill_type_label, " ", "") # remove spaces

  # set up the plot
  p <- ggplot(summary_df, aes(x = !!sym(x_axis), y = !!sym(y_axis),
                              fill = !!sym(fill_type))) +
    labs(title = paste(count_type, "Counts for All Samples and",
                       data_type_desc, "Data"),
         x = x_axis_label, y = y_axis_label)

  # split bars for reads/counts, stack for barcodes
  if (count_type == "Barcode") {
    p <- p + geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
      scale_fill_manual(name = fill_type_label,
                        values = unname(clrs_datatype))
  } else {
    p <- p + geom_bar(stat = "identity", position = "dodge",
                      color = "black", linewidth = 0.2) + #
      scale_fill_manual(name = fill_type_label,
                        values = unname(clrs_datatype[data_types]))
  }

  # remove scientific notation
  p <- p + scale_y_continuous(labels = label_comma(), breaks = breaks_pretty())

  # add remaining styling
  p <- p + theme_bw + labels_rotate_x + labels_standard

  return(p)
}


#' This function plots an overview of a doublet identification method
#'
#' @details
#' It assumes that named_colors$doublet has been defined.
#' Depends on other plots.
#' The doublets will be plotted "on top" for the first UMAP.
#'
#' @param seurat_obj The Seurat object.
#' @param tissue_type Blood, Skin.
#' @param clrs_specific The specific color palette (should be named).
#' @param use_hues Use the iwanthue hues instead of the default ggplot colors. Doesn't let you set any other settings.
#' @param group_col The column to group by.
#' @param group_label The label for the grouping variable to use in the plot titles and axis labels. If NULL, it will be determined based on the group_col name.
#' @param doublet_col The column containing the doublets information
#' @param doublet_package The doublet method being used.
#' @param details The optional subtitle.
#'
#' @returns A grid of four plots with UMAPs in the left column and bar plots in the right column.
#' @export
layout_doublets <- function(seurat_obj, tissue_type, clrs_specific,
                            use_hues = FALSE, group_col = "seurat_clusters",
                            group_label = NULL,
                            doublet_col = "scDblFinder.class",
                            doublet_package = "scDblFinder", details = NULL) {
  # if you want to use default ggplot2 or generated iwanthue colors
  if (rlang::is_missing(clrs_specific)) {
    n_colors <- n_distinct(seurat_obj[[group_col]])

    if (use_hues) clrs_specific <- hues::iwanthue(n_colors)
    else clrs_specific <- hue_pal()(n_colors)
  }

  # determine the legend label
  if (is.null(group_label)) {
    if (grepl("annotated", group_col)) {
      cluster_legend <- "Cell Type"
      # subclustered only
      if (group_col == "annotated_subclusters") {
        cluster_legend <- paste("Subclustered", cluster_legend)
      }
    } else if (grepl("cluster", group_col)) {
      cluster_legend <- "Cluster"
      # subclustered only
      if (group_col == "seurat_subclusters") {
        cluster_legend <- "Subcluster"
      }
    } else {
      cluster_legend <- group_col
    }
  } else {
    cluster_legend <- group_label
  }

  # set identities
  Idents(seurat_obj) <- group_col

  # UMAP colored by doublets/singlets
  p1 <- plot_umap(seurat_obj = seurat_obj, tissue_type = tissue_type,
                  clrs_specific = named_colors$doublet, specific_clusters = "doublet",
                  specific_col = doublet_col,
                  plot_by = "cluster_all", plot_label = FALSE, order = TRUE,
                  include_legend = FALSE) +
    labs(subtitle = details)

  # bar plot of doublets with raw counts
  p2 <- data.frame(table(seurat_obj[[]] %>%
                           select(all_of(group_col), all_of(doublet_col)))) %>%
    ggplot(aes(x = !!rlang::sym(group_col), y = Freq,
               fill = !!rlang::sym(doublet_col))) +
    geom_bar(stat = "identity", position = "dodge",
             color = "black", linewidth = 0.2) +
    geom_text(aes(label = Freq), position = position_dodge(width = 0.9),
              vjust = -1, size = 3) +
    labs(title = paste(tissue_type, "Doublets by", cluster_legend),
         subtitle = details, x = cluster_legend,
         y = "Count", fill = doublet_package) +
    scale_fill_manual(values = named_colors$doublet) +
    theme_bw + labels_standard

  # UMAP colored by clusters/annotations
  if (grepl("annotated", group_col)) {
    # rotate x-axis labels for annotation plots
    p2 <- p2 + labels_rotate_x

    p3 <- plot_umap(seurat_obj = seurat_obj, tissue_type = tissue_type,
                    clrs_specific = clrs_specific,
                    plot_by = "cluster_all", annotated = TRUE,
                    annotations_col = group_col,
                    include_legend = FALSE)
  } else {
    p3 <- plot_umap(seurat_obj = seurat_obj, tissue_type = tissue_type,
                    clrs_specific = clrs_specific,
                    plot_by = "cluster_all", clusters_col = group_col,
                    include_legend = FALSE)
  }

  # bar plot of doublets with percentage counts
  p4 <- plot_pcts(pcts = calc_pcts(data = seurat_obj[[]],
                                   meta_group_by = group_col,
                                   focus_group = doublet_col),
                  tissue_type = tissue_type, plot_type = "All",
                  plot_value = "Doublets",
                  x_axis = group_col, x_axis_label = cluster_legend,
                  fill_type = doublet_col, fill_label = doublet_package,
                  clrs_specific = named_colors$doublet, details = details)

  # combine plots
  p1 + p2 + p3 + p4 +
    plot_layout(guides = "collect", widths = c(1, 3)) + plot_anno &
    theme(plot.tag = element_text(face = "plain", size = 12))
}
