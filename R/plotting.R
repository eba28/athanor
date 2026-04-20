#' Plot UMAPs with AIRR overlays
#'
#' @description
#' This function plots the annotated UMAPs alongside the BCR/TCR overlays.
#'
#' @details
#' lightgray is Seurat's default background color and cells_total has to be a list for the labels to work later.
#' This could probably be replaced with `plot_umap()`.
#'
#' @param seurat_obj The Seurat object (must contain Has_BCR/Has_TCR cols)
#' @param tissue_type The tissue type of interest.
#' @param airr_type BCR or TCR
#' @param clrs_specific The specific color palette (should be named).
#' @param barcode_col The barcode column name: barcode, cell_id, or Cell_ID_Unique.
#' @param plot_label Whether or not to include the labels.
#' @param plot_by The grouping method: all samples, by dataset, by sample, by isotype.
#' @param ncol_sample The number of columns for sample-wise plots.
#'
#' @returns A Seurat UMAP plot with AIRR cells highlighted by the specified grouping.
#' @export
plot_immune_overlay <- function(seurat_obj, tissue_type, airr_type,
                                clrs_specific, plot_by = "all",
                                barcode_col = "cell_id",
                                plot_label = FALSE, ncol_sample = 6) {
  # plot options
  pt_size <- 0.2
  label_size <- 4
  sizes_highlight <- 0.2 # should probably be the same size as the points
  plot_title <- paste(tissue_type, airr_type, "Overlay")

  # if you want to use the default Seurat colors
  # if (rlang::is_missing(clrs_specific)) clrs_specific <- hue_pal()(num_clusts)

  # AIRR info
  airr_col <- paste0("Has_", airr_type)
  combined_vdj_gex <- filter(seurat_obj[[]], !!rlang::sym(airr_col))

  if (plot_by == "dataset" || plot_by == "sample") {
    # set up the cells to be highlighted
    cells_total <- c()
    datasets <- as.character(unique(combined_vdj_gex$Dataset))

    for (dataset in datasets) {
      cells_dataset <- filter(combined_vdj_gex,
                              Dataset == dataset)[[barcode_col]]
      cells_total[dataset] <- list(cells_dataset)
    }

    # main plot (coloring the samples by dataset is helpful)
    if (plot_by == "sample") {
      p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                    split.by = "sample_id",
                    repel = TRUE,
                    cells.highlight = cells_total,
                    sizes.highlight = sizes_highlight,
                    ncol = ncol_sample, raster = FALSE) +
        labs(title = paste(plot_title, "by Sample"), color = "Dataset")
    } else { # dataset
      p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                    label = plot_label, label.size = label_size,
                    repel = TRUE,
                    cells.highlight = cells_total,
                    sizes.highlight = sizes_highlight,
                    raster = FALSE) +
        labs(title = paste(plot_title, "by Dataset"), color = "Dataset")
    }

    # the legend
    p <- p + scale_color_manual(labels = c(paste0("non-", airr_type), rev(datasets)),
                                values = c("lightgray", unname(clrs_specific[rev(datasets)])))
  } else if (plot_by == "isotype") {
    # select the cells to be highlighted
    # Idents(seurat_obj) <- "isotype"
    # cells_total <- CellsByIdentities(seurat_obj)
    # cells_total <- lapply(cells_total, function(x) x[!is.na(x)])
    # cells_total["NA"] <- c()
    # isotypes <- rev(names(cells_total)) # need to reverse for plotting
    cells_total <- c()
    isotypes <- sort(unique(combined_vdj_gex$isotype))

    for (isotype in isotypes) {
      cells_isotype <- filter(combined_vdj_gex, isotype == isotype)[[barcode_col]]
      cells_total[isotype] <- list(cells_isotype)
    }

    # main plot
    p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                  label = plot_label, label.size = label_size,
                  repel = TRUE,
                  cells.highlight = cells_total,
                  sizes.highlight = sizes_highlight,
                  raster = FALSE) +
      labs(title = paste(plot_title, "by Isotype"), color = "Isotype")

    # the legend
    p <- p + scale_color_manual(labels = c(paste0("non-", airr_type), rev(isotypes)),
                                values = c("lightgray", unname(clrs_specific[rev(isotypes)])))
  } else { # all
    # set up the cells to be highlighted
    cells_total <- combined_vdj_gex[[barcode_col]]

    # plot all samples
    p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                  label = plot_label, label.size = label_size,
                  repel = TRUE,
                  cols.highlight = clrs_specific[airr_type],
                  cells.highlight = cells_total,
                  sizes.highlight = sizes_highlight,
                  raster = FALSE) +
      labs(title = plot_title, color = "Data Type")

    # the legend
    p <- p + scale_color_manual(name = "Data Type",
                                labels = c(paste0("non-", airr_type), airr_type),
                                values = c("lightgray", unname(clrs_specific[airr_type])))
  }

  # standardize the labels
  p <- p & labels_standard & clean_umap

  return(p)
}

#' Add an information bar on top of a given ggplot
#'
#' @description
#' This function adds informational strips (facets) to the top of an existing `ggplot2`.
#'
#' @details
#' This is very useful for splitting by dataset, showing cell types, etc.
#' `contains` is not really needed.
#' Note that providing a named list of features to a DotPlot will automatically group them without the need for this function.
#'
#' @param plot A ggplot2 object to which information bars will be added.
#' @param method The method for getting the info bar information; one of `add`, `contains`, `join`.
#' @param info_type What to add e.g. Dataset, Cell_Type, etc.
#' @param info Data frame containing the information to join/add, must contain
#'   the same column as plot$data to join by.
#' @param sort Sort the info alphanumerically.
#' @param text_size The size of the text in the bars.
#' @param label Whether or not to include what you are adding info for (good for datasets).
#' @param angle Rotation degree.
#' @param top_side Plot on top (by default) or on the right side.
#'
#' @returns A ggplot plot.
#'
#' @examples
#' \dontrun{
#' p <- ggplot(data, aes(x, y)) + geom_point()
#' add_info_bar(p, info_type = "Dataset", info = dataset_info)
#' }
#' @export
add_info_bar <- function(plot, method = "contains", info_type, info,
                         sort = TRUE, text_size = 8, label = TRUE, angle = 0,
                         top_side = TRUE) {
  # add the information to the plot data if needed
  if (method == "add") {
    plot$data[info_type] <- info[[info_type]]
  } else if (method == "join") {
    plot_data <- suppressMessages(left_join(plot$data, info,
                                            multiple = "all",
                                            relationship = "many-to-many"))

    # sort each group of features alphanumerically within each group
    # (especially useful for markers on a DotPlot)
    # might not be accurate for features that show up in multiple groups
    # if it wasn't for numbers, we could get rid of sort and just do
    if (sort) {
      plot_data <- plot_data %>%
        group_by(!!sym(info_type)) %>%
        mutate(features.plot =
                 factor(features.plot,
                        levels = str_sort(unique(features.plot),
                                          numeric = TRUE))) %>%
        ungroup()
    }

    plot$data <- plot_data
  }

  # explain the label or not
  if (label) {
    labels <- partial(label_both, sep = " ")
  } else {
    if (angle == 90) {
      # adjust the width as desired
      labels <- label_wrap_gen(width = 10, multi_line = TRUE)
    } else {
      labels <- labeller()
    }
  }

  # split up the plot
  if (top_side) {
    plot + facet_grid(cols = vars(!!rlang::sym(info_type)),
                      scales = "free_x", space = "free", labeller = labels) +
      theme(strip.text.x = element_text(size = text_size, angle = angle))
  } else {
    plot + facet_grid(rows = vars(!!rlang::sym(info_type)),
                      scales = "free_y", space = "free", labeller = labels) +
      theme(strip.text.y = element_text(size = text_size, angle = angle))
  }
}


#' Calculate percentages per metadata group
#'
#' @description
#' This function calculates the percentage of occurrences of a specified focus group within a dataset, grouped by specified metadata columns.
#'
#' @details
#' This is not very elegant, but seems to work (dealing with string args is weird).
#' Make sure to filter your data as needed (e.g. no NA isotypes) beforehand.
#' Made for cell types and isotypes.
#' The percentages will be doubles between 0 and 100.
#' It will fill in missing values (especially needed for ordering if provided).
#' You could alternatively use the fill option in geom_bar() instead.
#'
#' @param data A data.frame containing the data to be analyzed e.g. the metadata of a Seurat object.
#' @param meta_group_by A character vector specifying the columns to group by.
#' @param focus_group A character string specifying the column to focus on for percentage calculation.
#' @param order_by An optional character string specifying the focus group value to order the results by. If not provided, no ordering is applied.
#'
#' @returns A data.frame with the calculated percentages and counts for each group.
#' @export
calc_pcts <- function(data, meta_group_by = c("sample_id", "Dataset"),
                      focus_group, order_by) {
  data <- data %>%
    group_by(across(all_of(meta_group_by)), !!sym(focus_group)) %>%
    summarize(Count = n(), .groups = "drop_last") %>%
    dplyr::rename(focus = !!sym(focus_group)) %>%
    complete(focus = unique(data[[focus_group]]),
             fill = list(Count = 0)) %>%
    mutate(Percent = 100 * Count / sum(Count)) %>%
    dplyr::rename(!!sym(focus_group) := focus)

  # useful for showing a trend
  if (!rlang::is_missing(order_by)) {
    new_order <- data %>%
      filter(!!sym(focus_group) == order_by) %>%
      arrange(desc(Percent)) %>%
      pull(meta_group_by) %>%
      as.character()

    data[[meta_group_by]] <- fct_relevel(data[[meta_group_by]], new_order)
  }

  return(data)
}


#' Plot percentages in a stacked bar plot
#'
#' @description
#' This function plots a stacked bar plot of percentages calculated using `calc_pcts()` with percentages
#' labeled and total counts on top.
#'
#' @details
#' Give the percentages already as percents (* 100 in the calculations).
#' Make sure `pcts` includes Dataset if you want to split by dataset.
#' Note that the percentages may seem inaccurate because of the accuracy.
#' There was a big issue with the combo of `geom_text()` & `facet_grid()`.
#' This assumes that you want to show the counts for binary plots.
#'
#' @param pcts The output of `calc_pcts()`.
#' @param tissue_type The type of tissue being plotted e.g. Blood or Skin.
#' @param clrs_specific A specific (must have names) color palette.
#' @param plot_type One of `All`, `Binary`.
#' @param plot_value What is being plotted.
#' @param x_axis What to put along the x axis.
#' @param x_axis_label The label for the x axis.
#' @param fill_type What to group the bar plot by.
#' @param fill_label The description of what you're filling by.
#' @param perc_min The minimum percentage to show in the plot.
#' @param label_size The size of the percentage labels.
#' @param label_fill Add a white background for clarity.
#' @param include_counts Plot the counts on top/bottom.
#' @param drop_zeroes Remove percentages of zeroes.
#' @param reverse_order Change the fill order for a stacked plot.
#' @param total_order Rearrange x axis in descending order by totals instead of alphabetically.
#' @param details The optional subtitle.
#'
#' @returns A stacked ggplot bar plot
#' @export
plot_pcts <- function(pcts, tissue_type, clrs_specific,
                      plot_type = "All", plot_value = "Cell Type",
                      x_axis = "sample_id", x_axis_label = "Sample",
                      fill_type = "annotated_clusters", fill_label = fill_type,
                      perc_min = 3, label_size = 3, label_fill = FALSE,
                      include_counts = TRUE, drop_zeroes = TRUE,
                      reverse_order = FALSE, total_order = FALSE,
                      details = NULL) {
  # if you want to use the default ggplot2 colors
  if (rlang::is_missing(clrs_specific) || is.null(clrs_specific)) {
    clrs_specific <- hue_pal()(n_distinct(pcts %>%
                                            ungroup() %>%
                                            dplyr::select(all_of(fill_type))))
    clrs_specific <- setNames(clrs_specific, unique(pcts[[fill_type]]))
  }

  # mainly affects the legend
  if (drop_zeroes) pcts <- filter(pcts, Percent != 0)

  # set up the plotting layers properly
  if (plot_type == "All") {
    # add the total counts to the data
    pcts <- pcts %>%
      group_by(across(all_of(x_axis))) %>%
      mutate(Total = sum(Count), Total = replace(Total, -n(), ""))

    if (reverse_order) {
      pcts[[fill_type]] <- factor(pcts[[fill_type]],
                                  levels = rev(levels(factor(pcts[[fill_type]]))))
      # clrs_specific <- rev(clrs_specific)
    }

    if (total_order) {
      new_order <- pcts %>%
        mutate(Total = as.numeric(Total)) %>%
        arrange(desc(Total)) %>%
        filter(!is.na(Total)) %>%
        pull(x_axis) %>%
        as.character()

      pcts[[x_axis]] <- fct_relevel(pcts[[x_axis]], new_order)
    }
  } else if (plot_type == "Binary") {
    # so that TRUE plots on top of FALSE
    pcts[[fill_type]] <- factor(pcts[[fill_type]], levels = c(TRUE, FALSE))
  } else {
    stop("Please enter a valid plot type.")
  }

  # set up the plot
  p <- ggplot(data = pcts,
              aes(x = !!sym(x_axis), y = Percent, fill = !!sym(fill_type))) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    labs(title = paste(tissue_type, "Percent of",
                       plot_value, "per", x_axis_label),
         subtitle = details,
         x = x_axis_label, y = paste(plot_value, "Percentage"),
         fill = fill_label)

  # add aesthetics
  p <- p + scale_fill_manual(values = clrs_specific, limits = force) +
    scale_y_continuous(labels = scales::label_percent(scale = 1),
                       expand = expansion(mult = 0.05))

  # add labels
  if (!label_fill) {
    p <- p + geom_text(aes(label = scales::label_percent(accuracy = 1, scale = 1)(ifelse(Percent > perc_min, Percent, NA)),
                           group = !!sym(fill_type)),
                       size = label_size,
                       position = position_stack(vjust = 0.5))
  } else {
    p <- p + geom_label(aes(label = scales::label_percent(accuracy = 1, scale = 1)(ifelse(Percent > perc_min, Percent, NA)),
                            group = !!sym(fill_type)),
                        color = "black", fill = "white", size = label_size,
                        position = position_stack(vjust = 0.5))
  }

  if (include_counts) {
    # actually add in the total counts
    if (plot_type == "All") {
      p <- p + geom_text(mapping = aes(x = !!sym(x_axis), y = 100,
                                       label = Total, fill = NULL),
                         vjust = -0.5, color = "black", size = label_size)
    } else if (plot_type == "Binary") { # just for TRUE v FALSE plots
      p <- p +
        # counts on top
        geom_text(mapping = aes(label = ifelse(!!sym(fill_type) == TRUE,
                                               Count, "")),
                  position = position_stack(), vjust = -0.5,
                  color = "black", size = label_size) +
        # counts on bottom
        geom_text(mapping = aes(label = ifelse(!!sym(fill_type) == FALSE,
                                               Count, ""), y = -1),
                  color = "black", size = label_size)
    }
  }

  # add remaining styling
  p <- p + theme_bw + labels_standard

  # remove horizontal gridlines
  p <- p + theme(panel.grid.major.y = element_blank())

  # fix the legend if needed
  if (reverse_order) {
    p <- p + guides(fill = guide_legend(reverse = TRUE))
  }

  # rotate (typically) long names
  if (x_axis %in% c("sample_id", "annotated_clusters",
                    "annotated_subclusters")) {
    p <- p + labels_rotate_x
  }

  return(p)
}


#' Generate a color scale for a ggplot2 plot with white at zero
#'
#' @description
#' It creates a gradient of colors based on the range of values in the plot or data provided, and applies it to the specified color or fill aesthetic if applicable.
#'
#' @details
#' Seurat's `col` option frequently is misleading with where the zeroes fall.
#' Works well for making a scale for a Seurat `DotPlot` that accurately reflects the expression value.
#' You could also just do something like `scale_color_gradient2(low = "#2166AC",
#' mid = "white", high = "#B2182B")`
#' Can also be used to generated a color scale for a general ggplot2.
#' Use the function via a pipe right after the function call.
#'
#' @param plot The generated Seurat DotPlot or ggplot.
#' @param data The data used to generate the plot.
#' @param val_col The column in the plot data that contains the values to be plotted (e.g. "avg.exp.scaled").
#' @param palette A palette of colors to go off of.
#' @param fill_by One of "color" or "fill".
#'
#' @returns A Seurat dot plot or ggplot with an updated color scale, or just a vector of colors if a plot is not provided.
#' @export
plot_color_scale <- function(plot, data, val_col = "avg.exp.scaled",
                             palette = rev(pals::brewer.rdbu(n = 7)),
                             fill_by = "color") {
  # TODO: just return the scale if plot is not provided and some more values are given

  if (!rlang::is_missing(plot) & !rlang::is_missing(data)) {
    stop("Please only provide either a plot or the data to be plotted, not both.")
  }

  # use data from the plot if one is provided
  if (!rlang::is_missing(plot)) {
    plot_data <- plot$data[[val_col]]
  } else if (!rlang::is_missing(data)) {
    plot_data <- data

    # using base R instead of rlang so that it doesn't consider the default val
    if (!missing(val_col)) plot_data <- plot_data[[val_col]]
  } else {
    stop("Please provide either a plot or the data to be plotted.")
  }

  # expression ranges
  max_val <- max(plot_data, na.rm = TRUE)
  min_val <- min(plot_data, na.rm = TRUE)

  # make a full range of colors
  nclrs <- 100
  abs_val <- max(abs(max_val), abs(min_val))
  breaks <- seq(-abs_val, abs_val, length.out = nclrs + 1) # should be odd
  colors <- colorRampPalette(palette)(length(breaks))
  names(colors) <- breaks
  colors["0"] <- "white" # manually set zero to be white

  # now cut the range to only cover the actual expression
  # might be slightly off from the actual min/max value, but whatever
  if (abs(max_val) >= abs(min_val)) {
    min_val_loc <- which.min(abs(breaks - min_val))
    colors <- colors[min_val_loc:(nclrs + 1)]
  } else {
    max_val_loc <- which.min(abs(breaks - max_val))
    colors <- colors[0:max_val_loc]
  }

  if (!rlang::is_missing(plot)) {
    # return the plot with the new color scale
    if (fill_by == "color") plot + scale_color_gradientn(colors = colors)
    else plot + scale_fill_gradientn(colors = colors)
  } else {
    # just return the scale
    return(unname(colors))
  }
}


#' Plot a Seurat `VlnPlot` and a `FeaturePlot` side by side for the same marker
#'
#' @description
#' This function generates a side-by-side visualization of a Seurat `VlnPlot` and a `FeaturePlot` for a specified marker, allowing for a comprehensive comparison of expression levels across different groups and spatial distribution on the UMAP.
#' The `VlnPlot` will display the distribution of expression levels across specified groups, while the `FeaturePlot` will show the spatial localization of the marker on the UMAP, with higher expressing cells highlighted on top for better visibility.
#'
#' @details
#' Will put the highest expressing cells on top for the latter.
#'
#' @param seurat_obj The Seurat object with GEX data.
#' @param feature The feature of interest.
#' @param assay The assay to search for the feature in.
#' @param group_col What to group by (uses the Idents by default).
#' @param rotate Rotate the labels or not.
#'
#' @returns Two patchworked Seurat plots.
#' @export
plot_vln_feat <- function(seurat_obj, feature, assay = "RNA", group_col = NULL,
                          rotate = FALSE) {
  # set the assay and idents
  if (!is.null(group_col)) Idents(seurat_obj) <- group_col
  DefaultAssay(seurat_obj) <- assay

  p1 <- VlnPlot(object = seurat_obj, features = feature, pt.size = 0.1,
                raster = FALSE) +
    NoLegend()

  if (rotate) p1 <- p1 + labels_standard_vln
  else p1 <- p1 + labels_standard_vln_rotate

  p2 <- FeaturePlot(object = seurat_obj, features = feature, pt.size = 0.1,
                    order = TRUE, min.cutoff = 0, label = TRUE, label.size = 3,
                    raster = FALSE) +
    scale_color_viridis_c(option = "G") +
    labels_standard + clean_umap

  (p1 | p2) & plot_layout(nrow = 1, widths = c(2, 1))
}


#' Plot Seurat UMAP(s) in several useful ways
#'
#' @description
#' This function generates UMAP plots from a Seurat object with various customizable options for coloring, labeling, and grouping the data.
#' It allows for flexible visualization based on different metadata columns, cluster annotations, and specific clusters of interest, making it easier to explore and interpret the underlying structure of the data.
#'
#' @details
#' Includes whether or not the object has been annotated with specific cell types.
#'
#' @param seurat_obj The Seurat object.
#' @param tissue_type Blood, Skin.
#' @param clrs_specific The specific color palette (should be named).
#' @param use_hues Use the `iwanthue` hues instead of the default ggplot colors. Doesn't let you set any other settings.
#' @param plot_by What to plot by: by dataset, by sample type (control vs EM/blood), by cluster, by sample, with all samples, or with all subjects.
#' @param specific_clusters Can overlay clusters of interest e.g. B cell or by #.
#' @param specific_col Overlay a specific column in the object e.g. "annotated_clusters" or "sample_id". Overrides the other options.
#' @param plot_label Add labels to the plot (or not).
#' @param label_box Whether or not to give the labels a background.
#' @param label_size The size of the plot labels.
#' @param clusters_col Which column in the object stores the clusters.
#' @param annotated If the cell types have been identified.
#' @param annotations_col Which column in the object stores the cell types.
#' @param annotations_type The method for annotation e.g. Manual, singleR, etc..
#' @param order Plot cells on top or not.
#' @param include_legend Include the legend or not.
#' @param ncol The number of columns if outputting multiple plots.
#' @param details A custom subtitle.
#'
#' @returns A Seurat UMAPPlot.
#' @export
plot_umap <- function(seurat_obj, tissue_type = "", clrs_specific,
                      use_hues = FALSE, plot_by = "all", specific_clusters = c(),
                      specific_col, plot_label = TRUE, label_size = 3,
                      label_box = TRUE, clusters_col = "seurat_clusters",
                      annotated = FALSE, annotations_col = "annotated_clusters",
                      annotations_type, order = FALSE,
                      include_legend = TRUE, ncol = 4, details) {
  # make sure that the object is properly formatted
  # all(colnames(seurat_obj) == rownames(seurat_obj[[]]))

  # set identities in case they aren't already set
  if (annotated) {
    # don't include "Manual" in the title since we can consider it implied
    if (!rlang::is_missing(annotations_type)) {
      plot_title <- paste0("Annotated (", annotations_type, ") ", tissue_type)
    } else {
      plot_title <- paste("Annotated", tissue_type)
    }

    if (grepl("annotated_", annotations_col)) {
      cluster_legend <- "Cell Type"

      # special case
      if (annotations_col == "annotated_subclusters") {
        cluster_legend <- paste("Subclustered", cluster_legend)
      }
    } else {
      cluster_legend <- annotations_col
    }

    Idents(seurat_obj) <- annotations_col
    # alphabetize just in case
    Idents(seurat_obj) <- factor(Idents(seurat_obj),
                                 levels = sort(levels(seurat_obj)))
  } else {
    plot_title <- tissue_type
    cluster_legend <- "Cluster"

    # subclustered only
    if (clusters_col == "seurat_subclusters") {
      cluster_legend <- "Subcluster"
    }

    Idents(seurat_obj) <- clusters_col
  }

  # like an overlay
  if (!rlang::is_missing(specific_col)) {
    Idents(seurat_obj) <- specific_col
    cluster_legend <- specific_col
  }

  # plot options
  pt_size <- 0.2
  sizes_highlight <- 0.2

  # if you want to use default ggplot2 or generated iwanthue colors
  # this needs fixing
  if (rlang::is_missing(clrs_specific)) {
    if (plot_by == "all") {
      if (use_hues) clrs_specific <- hues::iwanthue(length(unique(seurat_obj$sample_id)))
      else clrs_specific <- hue_pal()(length(unique(seurat_obj$sample_id)))
    } else {
      if (use_hues) clrs_specific <- hues::iwanthue(nlevels(seurat_obj))
      else clrs_specific <- hue_pal()(nlevels(seurat_obj))

      clrs_specific <- setNames(clrs_specific, levels(Idents(seurat_obj)))
    }
  }

  # this option needs a little more setup
  if (plot_by == "cluster_specific") {
    # has to be a list for the labels to work later
    cells_total <- CellsByIdentities(seurat_obj, idents = specific_clusters)

    clrs_highlight <- clrs_specific[specific_clusters]
    clrs_highlight <- rev(clrs_highlight) # for some reason
  }

  # specific plots
  p <- switch(plot_by,
              "all" = UMAPPlot(seurat_obj, cols = clrs_specific,
                               pt.size = pt_size, order = order,
                               group.by = "sample_id",
                               label = plot_label,
                               label.size = label_size, label.box = label_box,
                               raster = FALSE) +
                labs(subtitle = "All Samples",
                     color = "Sample"),

              "subject" = UMAPPlot(seurat_obj, cols = clrs_specific,
                                   pt.size = pt_size, order = order,
                                   group.by = "subject_id",
                                   label = plot_label,
                                   label.size = label_size, label.box = label_box,
                                   raster = FALSE) +
                labs(subtitle = "All Subjects",
                     color = "Subject"),

              "dataset" = UMAPPlot(seurat_obj, cols = clrs_specific,
                                   pt.size = pt_size, order = order,
                                   group.by = "Dataset",
                                   label.box = label_box, raster = FALSE) +
                labs(subtitle = "by Dataset",
                     color = "Dataset"),

              "sample_type" = UMAPPlot(seurat_obj, cols = clrs_specific,
                                       pt.size = pt_size, order = order,
                                       group.by = "SampleType",
                                       label.box = label_box, raster = FALSE) +
                labs(subtitle = "by Sample Type",
                     color = "Sample Type"),

              "cluster_all" = UMAPPlot(seurat_obj, cols = clrs_specific,
                                       pt.size = pt_size, order = order,
                                       label = plot_label,
                                       label.size = label_size,
                                       label.box = label_box,
                                       repel = TRUE, raster = FALSE) +
                labs(subtitle = paste("by", cluster_legend),
                     color = cluster_legend),

              "cluster_sample" = UMAPPlot(seurat_obj, cols = clrs_specific,
                                          pt.size = pt_size, order = order,
                                          split.by = "sample_id",
                                          label.box = label_box, ncol = ncol,
                                          raster = FALSE) +
                labs(subtitle = paste(cluster_legend, "by Sample"),
                     color = cluster_legend),

              # only if annotated, plot all cells on top
              "cluster_specific" = UMAPPlot(seurat_obj, pt.size = pt_size,
                                            order = TRUE, label = plot_label,
                                            label.size = label_size,
                                            label.box = label_box,
                                            cells.highlight = cells_total,
                                            cols.highlight = clrs_highlight,
                                            sizes.highlight = sizes_highlight,
                                            repel = TRUE, raster = FALSE) +
                labs(subtitle = toString(specific_clusters),
                     color = "Highlighted"),

              "sample_type_split" = UMAPPlot(seurat_obj, cols = clrs_specific,
                                             pt.size = pt_size, order = order,
                                             split.by = "SampleType",
                                             label = plot_label, label.size = label_size,
                                             label.box = label_box, repel = TRUE,
                                             raster = FALSE) +
                labs(subtitle = paste(cluster_legend,
                                      "Split by Sample Type"),
                     color = cluster_legend)
  )

  # add the plot title
  p <- p + labs(title = plot_title)

  # don't plot super long subtitles (change nchar as needed)
  if (nchar(toString(specific_clusters)) > 80) p <- p + labs(subtitle = NULL)

  # custom subtitle
  if (!rlang::is_missing(details)) p <- p + labs(subtitle = details)

  # remove the legend if desired
  if (!include_legend) p <- p + NoLegend()

  # give white background to the boxes
  if (label_box) {
    p <- p + scale_fill_manual(values = rep("white", length(clrs_specific)))
    # okay this seems to be broken right now, so do this too
    if ("geom.use" %in% names(p@layers)) {
      p@layers$geom.use$aes_params$fill <- rep("white", length(clrs_specific))
    }
  }

  # standardize the labels
  p <- p & labels_standard & clean_umap

  return(p)
}


#' Plot a Seurat UMAP using `DimPlot`
#'
#' @description
#' This function generates a UMAP plot from a Seurat object using `DimPlot` with various customizable options for coloring, labeling, and grouping the data.
#'
#' @param seurat_obj The Seurat object.
#' @param data_source The dataset of origin.
#' @param clrs_specific Specific colors for plotting (make sure it has names).
#' @param use_hues Use the iwanthue hues instead of the default ggplot colors. Doesn't let you set any other settings.
#' @param pt_size The point size.
#' @param assay The data type e.g. ADT, GEX, BCR, WNN...
#' @param reduc The reduction to use for plotting e.g. wnn.umap
#' @param plot_label Add labels to the plot (or not).
#' @param label_box Whether or not to give the labels a background.
#' @param label_size The size of the plot labels.
#' @param annotated If the cell types have been identified.
#' @param specific_clusters Can overlay clusters of interest e.g. B cell or by #. Overrides the annotated option.
#' @param clusters_col Which column in the object stores the clusters.
#' @param annotations_col Which column in the object stores the cell types.
#' @param include_legend Include the legend or not.
#' @param legend_label The label for the legend.
#' @param factor_idents Whether or not to factorize the idents (for proper ordering of the colors). This can mess up the order you want, so be careful.
#' @param details A custom subtitle.
#' @param ... Any other Seurat parameters.
#'
#' @returns A Seurat UMAPPlot.
#' @export
plot_dimplot <- function(seurat_obj, data_source = "", clrs_specific,
                         use_hues = FALSE, pt_size = 0.2, assay, reduc,
                         plot_label = TRUE, label_box = TRUE, label_size = 3,
                         annotated = FALSE, specific_clusters, # order = FALSE
                         clusters_col = "seurat_clusters",
                         annotations_col = "annotated_clusters",
                         include_legend = TRUE, legend_label,
                         factor_idents = TRUE, details, ...) {
  # set identities in case they aren't already set (needed for group.by)
  if (annotated) {
    cluster_legend <- "Cell Type"
    Idents(seurat_obj) <- annotations_col

    # alphabetize (and factorize) just in case
    # be careful, sometimes this can mess up the order you want
    if (factor_idents) {
      Idents(seurat_obj) <-
        factor(Idents(seurat_obj), levels = sort(levels(seurat_obj)))
    }
  } else {
    cluster_legend <- "Cluster"
    Idents(seurat_obj) <- clusters_col

    # alphabetize (and factorize) just in case
    # be careful, sometimes this can mess up the order you want
    if (factor_idents) {
      Idents(seurat_obj) <-
        factor(Idents(seurat_obj),
               levels = str_sort(levels(seurat_obj), numeric = TRUE))
    }
  }

  # specific legend
  # TODO: use the column name instead
  if (!rlang::is_missing(legend_label)) cluster_legend <- legend_label

  # if you want to use default ggplot2 or generated iwanthue colors
  if (rlang::is_missing(clrs_specific)) {
    if (use_hues) clrs_specific <- hues::iwanthue(nlevels(seurat_obj))
    else clrs_specific <- hue_pal()(nlevels(seurat_obj))

    # the idents will have to be a factor
    clrs_specific <- setNames(clrs_specific, levels(Idents(seurat_obj)))
  }

  # to overlay clusters of interest
  if (!rlang::is_missing(specific_clusters)) {
    cells_total <- c()
    for (cluster in specific_clusters) {
      cells_cluster <- WhichCells(seurat_obj, idents = cluster)
      cells_total[cluster] <- list(cells_cluster) # has to be a list for the labels to work later
    }
    clrs_highlight <- clrs_specific[specific_clusters]
    clrs_highlight <- rev(clrs_highlight) # for some reason
  }

  # base plot
  if (rlang::is_missing(specific_clusters)) {
    p <- DimPlot(object = seurat_obj, cols = clrs_specific, pt.size = pt_size,
                 reduction = reduc, label = plot_label,
                 label.size = label_size, label.box = label_box,
                 repel = TRUE, na.value = "lightgray", raster = FALSE, ...) +
      labs(title = assay, subtitle = data_source, color = cluster_legend)
  } else {
    p <- DimPlot(object = seurat_obj, pt.size = pt_size,
                 reduction = reduc, order = TRUE, label = plot_label,
                 label.size = label_size, label.box = label_box,
                 cells.highlight = cells_total, cols.highlight = clrs_highlight,
                 sizes.highlight = pt_size,
                 repel = TRUE, na.value = "lightgray", raster = FALSE, ...) +
      labs(title = assay, subtitle = data_source, color = "Highlighted") # subtitle = toString(specific_clusters)
  }

  # custom subtitle
  if (!rlang::is_missing(details)) {
    p <- p + labs(title = paste(data_source, assay), subtitle = details)
  }

  # remove the legend if desired
  if (!include_legend) p <- p + NoLegend()

  # give white background to the boxes
  if (label_box) {
    p <- p + scale_fill_manual(values = rep("white", length(clrs_specific)))
    # okay this seems to be broken right now, so do this too
    if ("geom.use" %in% names(p@layers)) {
      p@layers$geom.use$aes_params$fill <- rep("white", length(clrs_specific))
    }
  }

  # standardize the labels
  p <- p + labels_standard & clean_umap

  return(p)
}


#' Plot a specific condition on a Seurat UMAP
#'
#' @description
#' This function generates a UMAP plot from a Seurat object using `DimPlot` with various customizable options for coloring, labeling, and grouping the data based on a specific condition of interest.
#'
#' @details
#' Based on: https://github.com/satijalab/seurat/issues/1053
#' Can be used for plotting QC metrics, isolating specific cell types, overlaying AIRR data, overlaying B cell isotypes, etc.
#'
#' @param seurat_obj The Seurat object.
#' @param tissue_type The tissue type e.g. Blood, Skin.
#' @param clrs_specific The specific color palette (should be named).
#' @param condition_name The column in the object that contains the condition of interest e.g. "annotated_clusters", "mu_freq", "isotype", etc.
#' @param operator <, >, ==
#' @param condition_val The value to compare the condition to. If color_by is "name", this should be a name in the condition_name column. If color_by is "value", this should be a value in the condition_name column.
#' @param color_by name, value
#' @param plot_type general, overlay (BCR/TCR)
#' @param label_plot Put labels on the plot (or not).
#' @param include_subtitle Include a subtitle (or not).
#' @param include_legend Include a legend (or not).
#'
#' @returns A Seurat UMAPPlot.
#' @export
plot_umap_condition <- function(seurat_obj, tissue_type, clrs_specific,
                                condition_name, operator, condition_val,
                                color_by = "value", plot_type = "general",
                                label_plot = TRUE, include_subtitle = TRUE,
                                include_legend = FALSE) {
  # plot options
  pt_size <- 0.2
  label_size <- 4
  sizes_highlight <- 0.2

  # cells of interest
  # condition <- FetchData(seurat_obj, vars = condition_name)
  seurat_obj <-
    switch(operator,
           "<" = subset(seurat_obj,
                        subset = !!sym(condition_name) < condition_val),
           ">" = subset(seurat_obj,
                        subset = !!sym(condition_name) > condition_val),
           "==" = subset(seurat_obj,
                         subset = !!sym(condition_name) == condition_val),
           stop("Invalid operator, please try again."))
  cells_highlight <- Cells(seurat_obj)

  # color by
  if (color_by == "name") {
    clrs_highlight <- clrs_specific[[condition_name]]
    plot_title <- paste(condition_name, operator, condition_val)
  } else if (color_by == "value") {
    clrs_highlight <- clrs_specific[[condition_val]]
    plot_title <- condition_val
  } else {
    cli::cli_abort("Invalid {.arg color_by} value: {.val {color_by}}. Please try again.")
  }

  # set up the plot
  p <- UMAPPlot(seurat_obj, pt.size = pt_size,
                label = label_plot, label.size = label_size,
                cells.highlight = cells_highlight,
                cols.highlight = clrs_highlight,
                sizes.highlight = sizes_highlight,
                raster = FALSE) +
    labs(title = paste(tissue_type, paste0("(", plot_title, ")"))) +
    # cols.highlight isn't working
    scale_color_manual(values = clrs_highlight)

  # include a subtitle counting the condition
  if (include_subtitle) {
    p <- p + labs(subtitle = paste("Count:", length(cells_highlight)))
  }

  # remove the legend if desired
  if (!include_legend) {
    p <- p + NoLegend()
  } else {
    if (plot_type == "overlay") {
      data_type <- substr(condition_name, nchar(condition_name) - 2,
                          nchar(condition_name)) # BCR or TCR

      p <- p + scale_color_manual(name = "Data Type",
                                  labels = c(paste0("non-", data_type), data_type),
                                  values = c("lightgray",
                                             unname(clrs_specific[data_type]))
      )
    } else {
      p <- p + scale_color_manual(name = paste(condition_name, operator, condition_val),
                                  labels = c("False", "True"),
                                  values = c("lightgray",
                                             unname(clrs_specific[condition_name]))
      )
    }
  }

  # standardize the labels
  p <- p + labels_standard + clean_umap

  return(p)
}


#' Plot several UMAPs side by side for a Seurat object
#'
#' @description
#' This function generates multiple UMAP plots from a Seurat object with various customizable options for coloring, labeling, and grouping the data based on different metadata columns, cluster annotations, and specific clusters of interest.
#' It allows for a comprehensive overview of the data across different embedding types (e.g. RNA, ADT, BCR) and comparisons (e.g. annotated clusters, V call families, isotypes), making it easier to explore and interpret the underlying structure of the data.
#'
#' @details
#' The names of the seurat_objs list correspond with embedding_types.
#' Assumes that CellTypist is the annotation approach being used.
#' `data_source` is set to empty to save on space.
#'
#' @param seurat_objs List of WNN objects with different embedding types.
#' @param data_source Dataset description.
#' @param pt_size The size of the points in the UMAP.
#' @param second_assay The second assay to use in the title if plotting a WNN reduction.
#' @param assay_name The name of the assay to use in the title. By default, it will be set based on the reduction (e.g. "GEX" for "rna.umap", "BCR" for "bcr.umap", and "GEX & BCR" for "wnn.umap").
#' @param reduction Which reduction to plot ("rna.umap", "bcr.umap", or "wnn.umap").
#' @param use_adt Whether or not the comparisons being plotted represent ADT markers.
#' @param comparisons The labeling of the plots.
#'
#' @return A patchwork object with overview plots in a grid.
#' @export
plot_overview_comps <- function(seurat_objs, data_source = "", pt_size = 0.1,
                                second_assay = "BCR", assay_name,
                                reduction = "wnn.umap", use_adt = FALSE,
                                comparisons = c("annotated_clusters_simpler",
                                                "v_call_family", "light_chains",
                                                "isotype", "mu_freq")) {
  # validate inputs
  # if (!reduction %in% c("rna.umap", "adt.umap", "bcr.umap", "wnn.umap")) {
  #   stop("reduction must be one of: 'rna.umap', 'adt.umap', 'bcr.umap', 'wnn.umap'")
  # }
  if (typeof(seurat_objs) == "S4") {
    # make a temp list so the rest of the code works
    seurat_objs <- list("obj" = seurat_objs)
  }

  # set assay name by the reduction
  # technically you could just look at the axis titles but people don't often
  # think to do that
  if (rlang::is_missing(assay_name)) {
    assay_name <- switch(reduction, "rna.umap" = "GEX", "adt.umap" = "ADT",
                         "bcr.umap" = "BCR",
                         "wnn.umap" = paste("GEX &", second_assay))
  }

  plots_overview <- list()
  for (type in names(seurat_objs)) {
    seurat_obj <- seurat_objs[[type]]

    # don't require using the embeddings approach
    if (!"embedding_type" %in% names(seurat_obj@misc) |
        reduction == "rna.umap") {
      details <- NULL
    } else {
      details <- embedding_types[[seurat_obj@misc$embedding_type]]
    }

    if (use_adt) {
      # plots_overview[[type]] <-
      #   # plot them together so that the scales are consistent
      #   suppressMessages(FeaturePlot(seurat_obj, features = str_c("adt_", comparisons),
      #                                order = TRUE, reduction = reduction, # keep.scale = "all",
      #                                ncol = 1, raster = FALSE) &
      #                      labs(subtitle = details) &
      #                      scale_color_viridis_c(option = "G", direction = -1) &
      #                      labels_standard & clean_umap)

      for (comparison in comparisons) {
        plots_overview[[paste0(comparison, "_", type)]] <-
          suppressMessages(
            FeaturePlot(seurat_obj, features = str_c("adt_", comparison),
                        order = TRUE, reduction = reduction, raster = FALSE) +
              labs(title = paste(assay_name, comparison), subtitle = details) +
              scale_color_viridis_c(name = comparison, option = "G",
                                    direction = -1) +
              labels_standard + clean_umap
            )
      }
    } else {
      # CellTypist
      if ("annotated_clusters" %in% comparisons) {
        plots_overview[[paste0("CellType_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = data_source, pt_size = pt_size,
                       clrs_specific = seurat_obj@misc$colors_annotated,
                       assay = assay_name,
                       reduc = reduction,
                       plot_label = FALSE, annotated = TRUE,
                       annotations_col = "annotated_clusters",
                       legend_label = "Cell Type", details = details)
      }

      if ("annotated_clusters_simpler" %in% comparisons) {
        plots_overview[[paste0("CellTypeSimpler_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = data_source, pt_size = pt_size,
                       clrs_specific = named_colors$cell_types_celltypist,
                       assay = assay_name,
                       reduc = reduction,
                       plot_label = FALSE, annotated = TRUE,
                       annotations_col = "annotated_clusters_simpler",
                       legend_label = "Cell Type", details = details)
      }

      # V call families
      if ("v_call_family" %in% comparisons) {
        plots_overview[[paste0("v_call_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$v_call_family,
                       assay = paste(assay_name, "V Call Families"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "v_call_family",
                       legend_label = "V Call Family", details = details)
      }

      # light chain types
      if ("light_chains" %in% comparisons) {
        plots_overview[[paste0("light_chain_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$light,
                       assay = paste(assay_name, "Light Chain Types"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "locus_light",
                       legend_label = "Light Chain Type", details = details)
      }

      # isotypes
      if ("isotype" %in% comparisons) {
        plots_overview[[paste0("isotype_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$isotype,
                       assay = paste(assay_name, "Isotypes"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "isotype",
                       legend_label = "Isotype", details = details)
      }

      # subisotypes
      if ("c_call" %in% comparisons) {
        plots_overview[[paste0("c_call_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$c_call,
                       assay = paste(assay_name, "Subisotypes"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "c_call",
                       legend_label = "Subisotype", details = details)
      }

      # SHM frequencies
      if ("mu_freq" %in% comparisons) {
        plots_overview[[paste0("mu_freq_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$mu_freq_bins,
                       assay = paste(assay_name, "SHM Frequencies (Binned)"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "mu_freq_bins",
                       legend_label = "Bins", details = details)
      }

      # CDR3 amino acid length
      if ("cdr3_aa_length" %in% comparisons) {
        plots_overview[[paste0("cdr3_aa_length_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$cdr3,
                       assay = paste(assay_name, "CDR3 Length"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "cdr3_aa_length",
                       legend_label = "CDR3 Length", details = details)
      }

      if ("weight_assay" %in% comparisons) {
        plots_overview[[paste0("weight_assay_", type)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = "", pt_size = pt_size,
                       clrs_specific = named_colors$weight_assay,
                       assay = paste(assay_name, "Modality Weights"),
                       reduc = reduction,
                       plot_label = FALSE, clusters_col = "weight_assay",
                       legend_label = "Chosen Assay", details = details)
      }
    }
  }

  # combine all of the plots
  if (length(seurat_objs) == 1) {
    wrap_plots(plots_overview, ncol = min(length(comparisons), 5)) + plot_anno
  } else {
    # if (use_adt) {
    #   wrap_plots(plots_overview, ncol = length(seurat_objs)) +
    #     plot_anno # + plot_layout(guides = "collect")
    # } else {
    #   wrap_plots(plots_overview, nrow = length(comparisons), byrow = FALSE) +
    #     plot_anno + plot_layout(guides = "collect")
    # }
    wrap_plots(plots_overview, nrow = length(comparisons), byrow = FALSE) +
      plot_anno + plot_layout(guides = "collect")
  }
}


#' Plot a box plot of modality weights per cell type
#'
#' @description
#' This function creates a box plot to visualize the distribution of modality weights (e.g. RNA vs. BCR) across different cell types or clusters in a post-WNN Seurat object.
#'
#' @details
#' Assumes annotated_clusters is a column
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param details Details to add to the plot title.
#' @param second_assay List of other assays run through WNN in order.
#' @param clrs_specific A specific (must have names) color palette.
#' @param split_by A meta.data column to split the box plots up by.
#' @param y_axis_label Label for the y-axis.
#'
#' @returns A ggplot with the distribution of weights
#' @export
plot_mws <- function(seurat_obj, details = "", second_assay = "BCR",
                     clrs_specific = named_colors$mu_freq_bins,
                     split_by = "mu_freq_bins",
                     y_axis_label = "SHM Frequency Bins") {
  main_assay <- ifelse(length(second_assay) > 1, second_assay[-1], second_assay)

  weight <- grep(paste0("^", main_assay, ".*\\.weight.*$"),
                 colnames(seurat_obj[[]]), value = TRUE)

  n_assay <- length(second_assay) + 1

  if (!is.factor(seurat_obj[[]][split_by])) {
    seurat_obj[[]][split_by] <- factor(seurat_obj[[]][[split_by]])
  }

  # don't require using the embeddings approach
  if ("embedding_type" %in% names(seurat_obj@misc)) {
    subtitle <- embedding_types[[seurat_obj@misc$embedding_type]]
  } else {
    subtitle <- NULL
  }

  p <- ggplot(seurat_obj[[]],
              aes(x = !!sym(weight), y = !!sym(split_by),
                  fill = !!sym(split_by))) +
         geom_boxplot(outlier.size = 0.5) +
         geom_jitter(size = 0.2) +
         labs(title = paste(details, "Weights by Cell Type"),
              subtitle = subtitle, x = "Weights", y = y_axis_label) +
         scale_fill_manual(values = clrs_specific) +
         facet_wrap(vars(annotated_clusters), scales = "fixed") +
         theme_bw + labels_standard + theme(legend.position = "none")

  if (n_assay == 2) {
      p <- p +
             geom_vline(xintercept = 0.50, linetype = "dashed") +
             scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                                labels = c("0 [GEX]", "0.25", "0.50", "0.75",
                                           paste0("1 [", second_assay, "]")),
                                limits = c(0, 1))
  } else if (n_assay == 3) {
      p <- p +
             geom_vline(xintercept = 1/3, linetype = "dashed") +
             geom_vline(xintercept = 2/3, linetype = "dashed") +
             scale_x_continuous(breaks = seq(0, 1, by = 0.1),
                                labels = c("0.0 [GEX]", "0.1", "0.2",
                                           paste0("0.3 [GEX:", second_assay[1], "]"),
                                           "0.4", "0.5",
                                           paste0("0.6 [", second_assay[1], ":",
                                                  second_assay[2], "]"),
                                           "0.7", "0.8", "0.9",
                                           paste0("1 [", second_assay[2], "]")),
                                limits = c(0, 1))
  }

  p
}


#' Plot an overview of a doublet identification method
#'
#' @description
#' This function creates a grid of four plots to visualize the results of a doublet identification method.
#' The left column contains UMAP plots colored by doublet/singlet status and by clusters/annotations, while the right column contains bar plots showing the counts and percentages of doublets across clusters or annotations.
#'
#' @details
#' It assumes that `named_colors$doublet` has been defined.
#' Depends on other plots.
#' The doublets will be plotted "on top" for the first UMAP.
#'
#' @param seurat_obj The Seurat object.
#' @param tissue_type Blood, Skin.
#' @param clrs_specific The specific color palette (should be named).
#' @param use_hues Use the `iwanthue` hues instead of the default ggplot colors. Doesn't let you set any other settings.
#' @param group_col The column to group by.
#' @param group_label The label for the grouping variable to use in the plot titles and axis labels. If NULL, it will be determined based on the group_col name.
#' @param doublet_col The column containing the doublets information
#' @param doublet_package The doublet method being used.
#' @param details The optional subtitle.
#'
#' @returns A grid of four plots with UMAPs in the left column and bar plots in the right column.
#' @export
plot_doublets <- function(seurat_obj, tissue_type, clrs_specific,
                          use_hues = FALSE, group_col = "seurat_clusters",
                          group_label = NULL, doublet_col = "scDblFinder.class",
                          doublet_package = "scDblFinder", details = NULL) {
  # TODO: get rid of this function?

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
  # TODO: add an option to not plot the labels
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
