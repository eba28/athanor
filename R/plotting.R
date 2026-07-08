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
#' p <- Seurat::DotPlot(seurat_obj, features = c("CD19", "MS4A1"))
#' add_info_bar(p, method = "join", info_type = "Dataset", info = dataset_df)
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
#' It can be used in conjunction with [plot_pcts()] to visualize the results.
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
#'
#' @examples
#' df <- data.frame(
#'   sample_id = rep(c("S1", "S2"), each = 50),
#'   cell_type = c(rep(c("B", "T", "NK"), times = c(20, 20, 10)),
#'                 rep(c("B", "T", "NK"), times = c(10, 30, 10)))
#' )
#' calc_pcts(df, meta_group_by = "sample_id", focus_group = "cell_type")
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
#'
#' @examples
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, color = mpg)) +
#'   ggplot2::geom_point()
#' plot_color_scale(plot = p, val_col = "mpg")
#' @export
plot_color_scale <- function(plot, data, val_col = "avg.exp.scaled",
                             palette = rev(pals::brewer.rdbu(n = 7)),
                             fill_by = "color") {
  # TODO: just return the scale if plot is not provided and some more values are given

  if (!rlang::is_missing(plot) & !rlang::is_missing(data)) {
    cli::cli_abort("Please only provide either a plot or the data to be plotted, not both.")
  }

  # use data from the plot if one is provided
  if (!rlang::is_missing(plot)) {
    plot_data <- plot$data[[val_col]]
  } else if (!rlang::is_missing(data)) {
    plot_data <- data

    # using base R instead of rlang so that it doesn't consider the default val
    if (!missing(val_col)) plot_data <- plot_data[[val_col]]
  } else {
    cli::cli_abort("Please provide either a plot or the data to be plotted.")
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


#' Plot a Seurat UMAP using `DimPlot`
#'
#' @description
#' This function generates a UMAP plot from a Seurat object using `DimPlot` with various customizable options for coloring, labeling, and grouping the data.
#'
#' @param seurat_obj The Seurat object.
#' @param data_source Dataset description.
#' @param clrs_specific The specific color palette (should be named).
#' @param use_hues Use the iwanthue hues instead of the default ggplot colors. Doesn't let you set any other settings.
#' @param pt_size The point size.
#' @param plot_title The plot title.
#' @param reduc The reduction to use for plotting e.g. "bpca" or wnn.umap".
#' @param meta_col Which column in the object metadata to color by. When combined with `highlight`, highlights those values as an overlay instead of coloring all cells.
#' @param highlight Can overlay clusters of interest e.g. B cell or by #. Overrides the annotated option.
#' @param plot_label Add labels to the plot (or not).
#' @param label_size The size of the plot labels.
#' @param label_box Whether or not to give the labels a background.
#' @param include_legend Include the legend or not.
#' @param legend_label The label for the legend.
#' @param sort_idents Whether or not to sort the idents (for proper ordering of the colors). This can mess up the order you want, so be careful.
#' @param idents_char If sorting idents, whether to sort them as characters or numerically (e.g. cluster 10 should be after cluster 9, not before).
#' @param order Plot cells on top or not.
#' @param details An optional custom subtitle.
#' @param fix_aspect Fix the aspect ratio to 1:1 via `clean_dimplot2`.
#' @param simplify_titles Whether or not to convert x and y axis titles to simplified versions e.g. "UMAP_1" instead of "bcrUMAP_1".
#' @param ... Any other Seurat parameters.
#'
#' @returns A Seurat plot of the specified reduction.
#' @export
plot_dimplot <- function(seurat_obj, data_source = "", clrs_specific,
                         use_hues = FALSE, pt_size = 0.2, plot_title,
                         reduc = "rna.umap", meta_col = "annotated_clusters",
                         highlight, plot_label = FALSE, label_size = 3,
                         label_box = TRUE, include_legend = TRUE, legend_label,
                         sort_idents = TRUE, idents_char = TRUE, order = FALSE,
                         details, fix_aspect = TRUE, simplify_titles = FALSE,
                         ...) {
  # check parameters
  if (!meta_col %in% names(seurat_obj[[]])) {
    cli::cli_abort("{meta_col} is not a valid metadata column name. Please select one of: {names(seurat_obj[[]])}.")
  }
  if (!reduc %in% names(seurat_obj@reductions)) {
    cli::cli_abort(c("Reduction '{reduc}' not found in the provided Seurat object. Please make sure it is present and try again. Available reductions: {names(seurat_obj@reductions)}"))
  }


  # set identities in case they aren't already set
  Idents(seurat_obj) <- meta_col # instead of using group.by in DimPlot
  # factors help with plotting
  # be careful, sometimes this can mess up the order you want
  if (sort_idents) {
    sort_fn <- if (idents_char) sort else function(x) str_sort(x, numeric = TRUE)

    unique_idents <- as.character(unique(Idents(seurat_obj)))
    Idents(seurat_obj) <- factor(Idents(seurat_obj),
                                 levels = sort_fn(unique_idents))
  } else {
    # if it's already a factor, leave the levels as is
    if (!is.factor(Idents(seurat_obj))) {
      Idents(seurat_obj) <- factor(Idents(seurat_obj))
    }
  }

  # make the title and legend labels nicer if not provided
  # e.g. "annotated_clusters" -> "Annotated Clusters"
  meta_col <- stringr::str_to_title(stringr::str_replace_all(meta_col, "_", " "))

  # fill in missing arguments if needed
  if (is.null(data_source) || data_source == "") data_source <- NULL # don't show the subtitle
  if (rlang::is_missing(title)) plot_title <- meta_col
  # if you want to use default ggplot2 or generated iwanthue colors
  if (rlang::is_missing(clrs_specific)) {
    if (use_hues) clrs_specific <- hues::iwanthue(nlevels(seurat_obj))
    else clrs_specific <- scales::hue_pal()(nlevels(seurat_obj))

    # the idents will have to be a factor
    clrs_specific <- setNames(clrs_specific, levels(Idents(seurat_obj)))
  }

  # base plot
  base_args <- list(object = seurat_obj, pt.size = pt_size, reduction = reduc,
                    label = plot_label, label.size = label_size,
                    label.box = label_box, repel = TRUE, na.value = "lightgray",
                    raster = FALSE, ...)
  if (rlang::is_missing(highlight)) {
    extra_args <- list(cols = clrs_specific, order = order)
    legend_label <-
      if (!rlang::is_missing(legend_label)) legend_label else meta_col
  } else {
    # TODO: don't show labels for unselected groups
    cells_total <- CellsByIdentities(seurat_obj, idents = highlight)
    cells_total <- cells_total[highlight] # removes `NA`s
    clrs_highlight <- clrs_specific[highlight]
    clrs_highlight <- rev(clrs_highlight) # for some reason

    # deal with unnamed colors if needed
    if (length(clrs_specific) == 1) {
      if (is.na(clrs_highlight)) {
        clrs_highlight <- setNames(clrs_specific, highlight)
      }
    }

    # if specific clusters are being plotted, order them to be on top
    extra_args <- list(cells.highlight = cells_total,
                       cols.highlight = clrs_highlight,
                       sizes.highlight = pt_size, order = TRUE)
    legend_label <-
      if (!rlang::is_missing(legend_label)) legend_label else "Highlighted"
  }

  # start assembling the plot
  p <- do.call(DimPlot, c(base_args, extra_args)) +
    labs(title = plot_title, subtitle = data_source, color = legend_label)

  # custom subtitle if provided (otherwise the data source is the subtitle)
  if (!rlang::is_missing(details)) {
    p <- p + labs(title = paste(data_source, plot_title), subtitle = details)
  }

  # remove the legend if desired
  if (!include_legend) p <- p + NoLegend()

  # give white background to the boxes (skip in highlight modes)
  if (label_box && rlang::is_missing(highlight)) {
    # Idents will always be a factor due to previous code in this function,
    # so using nlevels is safe
    p <- p + scale_fill_manual(values = rep("white", nlevels(Idents(seurat_obj))))
    # okay this seems to be broken right now, so do this too
    if ("geom.use" %in% names(p@layers)) {
      p@layers$geom.use$aes_params$fill <- rep("white", nlevels(Idents(seurat_obj)))
    }
  }

  # standardize the labels
  p <- p + labels_standard
  if (fix_aspect) {
    p <- p + clean_dimplot2
  } else {
    p <- p + clean_dimplot
  }

  # TODO: add tSNE support
  if (simplify_titles) {
    if (grepl("umap", tolower(reduc))) {
      p <- p + labs(x = "UMAP1", y = "UMAP2")
    } else if (grepl("pca", tolower(reduc))) {
      p <- p + labs(x = "PCA1", y = "PCA2")
    } else {
      p <- p + labs(x = paste0(reduc, "1"), y = paste0(reduc, "2"))
    }
  }

  return(p)
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
#' @param data_source Dataset description.
#' @param clrs_specific The specific color palette (should be named).
#' @param use_hues Use the `iwanthue` hues instead of the default ggplot colors. Doesn't let you set any other settings.
#' @param reduc The reduction to use for plotting e.g. "bpca" or wnn.umap".
#' @param meta_col The column to group by.
#' @param group_label The label for the grouping variable to use in the plot titles and axis labels. If NULL, it will be determined based on the meta_col name.
#' @param doublet_col The column containing the doublets information
#' @param doublet_package The doublet method being used.
#' @param details An optional custom subtitle.
#'
#' @returns A grid of four plots with UMAPs in the left column and bar plots in the right column.
#' @export
plot_doublets <- function(seurat_obj, data_source = "", clrs_specific,
                          use_hues = FALSE, reduc = "rna.umap",
                          meta_col = "seurat_clusters",
                          group_label = NULL, doublet_col = "scDblFinder.class",
                          doublet_package = "scDblFinder", details = NULL) {
  # TODO: get rid of this function?

  # if you want to use default ggplot2 or generated iwanthue colors
  if (rlang::is_missing(clrs_specific)) {
    n_colors <- n_distinct(seurat_obj[[meta_col]])

    if (use_hues) clrs_specific <- hues::iwanthue(n_colors)
    else clrs_specific <- hue_pal()(n_colors)
  }

  # determine the legend label
  if (is.null(group_label)) {
    if (grepl("annotated", meta_col)) {
      cluster_legend <- "Cell Type"
      # subclustered only
      if (meta_col == "annotated_subclusters") {
        cluster_legend <- paste("Subclustered", cluster_legend)
      }
    } else if (grepl("cluster", meta_col)) {
      cluster_legend <- "Cluster"
      # subclustered only
      if (meta_col == "seurat_subclusters") {
        cluster_legend <- "Subcluster"
      }
    } else {
      cluster_legend <- meta_col
    }
  } else {
    cluster_legend <- group_label
  }

  # set identities
  Idents(seurat_obj) <- meta_col

  # UMAP colored by doublets/singlets
  p1 <- plot_dimplot(seurat_obj = seurat_obj, plot_title = data_source,
                     clrs_specific = named_colors$doublet,
                     reduc = reduc, highlight = "doublet",
                     meta_col = doublet_col,
                     plot_label = FALSE, order = TRUE,
                     include_legend = FALSE) +
    labs(subtitle = details)

  # bar plot of doublets with raw counts
  p2 <- data.frame(table(seurat_obj[[]] %>%
                           select(all_of(meta_col), all_of(doublet_col)))) %>%
    ggplot(aes(x = !!rlang::sym(meta_col), y = Freq,
               fill = !!rlang::sym(doublet_col))) +
    geom_bar(stat = "identity", position = "dodge",
             color = "black", linewidth = 0.2) +
    geom_text(aes(label = Freq), position = position_dodge(width = 0.9),
              vjust = -1, size = 3) +
    labs(title = paste(data_source, "Doublets by", cluster_legend),
         subtitle = details, x = cluster_legend,
         y = "Count", fill = doublet_package) +
    scale_fill_manual(values = named_colors$doublet) +
    theme_bw_custom + labels_standard

  # UMAP colored by clusters/annotations
  # TODO: add an option to not plot the labels
  if (grepl("annotated", meta_col)) {
    # rotate x-axis labels for annotation plots
    p2 <- p2 + labels_rotate_x

    p3 <- plot_dimplot(seurat_obj = seurat_obj, plot_title = data_source,
                       clrs_specific = clrs_specific, reduc = reduc,
                       meta_col = meta_col, include_legend = FALSE)
  } else {
    p3 <- plot_dimplot(seurat_obj = seurat_obj, plot_title = data_source,
                       clrs_specific = clrs_specific, reduc = reduc,
                       meta_col = meta_col, include_legend = FALSE)
  }

  # bar plot of doublets with percentage counts
  p4 <- plot_pcts(pcts = calc_pcts(data = seurat_obj[[]],
                                   meta_group_by = meta_col,
                                   focus_group = doublet_col),
                  data_source = data_source, plot_type = "All",
                  plot_value = "Doublets",
                  x_axis = meta_col, x_axis_label = cluster_legend,
                  fill_type = doublet_col, fill_label = doublet_package,
                  clrs_specific = named_colors$doublet, details = details)

  # combine plots
  p1 + p2 + p3 + p4 +
    plot_layout(guides = "collect", widths = c(1, 3)) + plot_anno &
    theme(plot.tag = element_text(face = "plain", size = 12))
}


#' Plot a box plot of modality weights per cell type
#'
#' @description
#' This function creates a box plot to visualize the distribution of modality weights (e.g. RNA vs. BCR) across different cell types or clusters in a post-WNN Seurat object.
#' It should be run after [run_wnn()].
#'
#' @details
#' Assumes "annotated_clusters" is a column.
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param data_source Dataset description.
#' @param second_assay The second assay to compare against RNA (e.g. "BCR"). If multiple assays are provided, it will assume the weight column is in the format "main_assay.weight" (e.g. "BCR.weight").
#' @param clrs_specific The specific color palette (should be named).
#' @param split_by A meta.data column to split the box plots up by.
#' @param facet_by A meta.data column to facet the box plots by.
#' @param y_axis_label Label for the y-axis.
#' @param details An optional custom subtitle.
#'
#' @returns A ggplot with the distribution of weights
#' @export
plot_mws <- function(seurat_obj, data_source = "", second_assay = "BCR",
                     clrs_specific = named_colors$mu_freq_bins,
                     split_by = "mu_freq_bins",
                     facet_by = "annotated_clusters_simpler",
                     y_axis_label = "SHM Frequency Bins", details = NULL) {
  # parameter check
  if (!split_by %in% names(seurat_obj[[]])) {
    cli::cli_abort("{split_by} is not a valid metadata column name. Please select one of: {names(seurat_obj[[]])}.")
  }

  main_assay <- ifelse(length(second_assay) > 1, second_assay[-1], second_assay)

  # could return BCR.weight and BCR_equal
  # weight <- grep(paste0("^", main_assay, ".*\\.weight.*$"),
  #                colnames(seurat_obj[[]]), value = TRUE)

  # assuming the weight colname is in the format "main_assay.weight"
  weight <- paste0(main_assay, ".weight")
  if (!weight %in% names(seurat_obj[[]])) {
    cli::cli_abort("{weight} is not a valid metadata column name. \\
                   Please make sure you have run WNN and that the weight \\
                   column is named in the format 'main_assay.weight' \\
                   (e.g. 'RNA.weight' or 'BCR.weight'). \\
                   Available columns: {names(seurat_obj[[]])}.")
  }

  n_assay <- length(second_assay) + 1

  if (!is.factor(seurat_obj[[]][split_by])) {
    seurat_obj[[]][split_by] <- factor(seurat_obj[[]][[split_by]])
  }

  p <- ggplot(seurat_obj[[]],
              aes(x = !!sym(weight), y = !!sym(split_by),
                  fill = !!sym(split_by))) +
         geom_boxplot(outlier.size = 0.5) +
         geom_jitter(size = 0.2) +
         labs(title = paste(data_source, "Weights by Cell Type"),
              subtitle = details, x = "Weights", y = y_axis_label) +
         scale_fill_manual(values = clrs_specific) +
         facet_wrap(vars(!!sym(facet_by)), scales = "fixed") +
         theme_bw_custom + labels_standard + theme(legend.position = "none")

  if (n_assay == 2) {
    p <-
      p +
      geom_vline(xintercept = 0.50, linetype = "dashed") +
      scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                        labels = c("0 [GEX]", "0.25", "0.50", "0.75",
                                   paste0("1 [", second_assay, "]")),
                        limits = c(0, 1))
  } else if (n_assay == 3) {
    p <-
      p +
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


#' Plot several UMAPs side by side
#'
#' @description
#' This function generates multiple UMAP plots from a Seurat object with various customizable options for coloring, labeling, and grouping the data based on different metadata columns, cluster annotations, and specific clusters of interest.
#' It allows for a comprehensive overview of the data across different embedding types (e.g. RNA, ADT, BCR) and comparisons (e.g. annotated clusters, V call families, isotypes), making it easier to explore and interpret the underlying structure of the data.
#'
#' @param seurat_objs List of Seurat objects.
#' @param data_source Dataset description.
#' @param pt_size The size of the points in the UMAP.
#' @param second_assay The second assay to use in the title if plotting a combined reduction. Will usually be "BCR".
#' @param assay_name The name of the assay to use in the title. By default, it will be set based on the reduction (e.g. "GEX" for "rna.umap", "BCR" for "bcr.umap", and "GEX & BCR" for "wnn.umap").
#' @param reduction Which reduction to plot (e.g. "rpca", "bcr.umap", "wnn.umap").
#' @param use_adt Whether or not the comparisons being plotted represent ADT markers.
#' @param ncol The number of columns to use in the grid.
#' @param comparisons Which metadata columns to plot. By default, it will plot "annotated_clusters_simpler", "v_call_family", "light_chains", "isotype", and "mu_freq". The first one is the simplified CellTypist annotations, and the rest are BCR features.
#' @param details_col Which column in `seurat_obj@misc` to use for the plot subtitles. By default, it will use "embedding_type" to show the type of embedding being plotted (e.g. "AntiBERTa2").
#'
#' @return A patchwork object with overview plots in a grid.
#' @export
plot_overview_comps <- function(seurat_objs, data_source = "", pt_size = 0.1,
                                second_assay = "BCR", assay_name,
                                reduction = "wnn.umap", use_adt = FALSE, ncol,
                                comparisons = c("annotated_clusters_simpler",
                                                "v_call_family", "light_chains",
                                                "isotype", "mu_freq"),
                                details_col = "embedding_type", ...) {
  possible_comps <- c("annotated_clusters_simpler", "annotated_clusters",
                      "c_call", "cdr3_aa_length", "isotype", "light_chains",
                      "mu_freq", "v_call_family", "weight_assay")

  # validate inputs
  # if (!reduction %in% c("rna.umap", "adt.umap", "bcr.umap", "wnn.umap")) {
  #   cli::cli_abort("reduction must be one of: 'rna.umap', 'adt.umap', 'bcr.umap', 'wnn.umap'")
  # }
  if (is.null(seurat_objs)) {
    cli::cli_abort("seurat_objs is NULL. Please pass a Seurat object or a list of Seurat objects.")
  }
  if (typeof(seurat_objs) == "S4") {
    # make a temp list so the rest of the code works
    seurat_objs <- list("obj" = seurat_objs)
  }
  if (is.null(names(seurat_objs))) {
    names(seurat_objs) <- paste("Object", seq_along(seurat_objs))
  }

  # resolve per-object titles (used as the data_source/subtitle for each object's plots)
  obj_titles <- setNames(
    lapply(names(seurat_objs), function(nm) if (nm == "obj") data_source else nm),
    names(seurat_objs)
  )

  # set assay name by the reduction
  # technically you could just look at the axis titles but people don't often
  # think to do that
  if (rlang::is_missing(assay_name)) {
    assay_name <- switch(reduction, "rpca" = "GEX", "rna.umap" = "GEX",
                         "bpca" = "BCR", "bcr.umap" = "BCR",
                         "adt.umap" = "ADT",
                         "rna_bcr.umap" = paste("GEX &", second_assay),
                         "wnn.umap" = paste("GEX &", second_assay))
  }

  plots_overview <- list()
  for (obj_name in names(seurat_objs)) {
    seurat_obj <- seurat_objs[[obj_name]]
    obj_data_source <- obj_titles[[obj_name]]

    if (!reduction %in% names(seurat_obj@reductions)) {
      cli::cli_abort("Reduction {reduction} not found in Seurat object {obj_name}. Please make sure it is present and try again.")
    }

    # plot an informative subtitle if possible
    if (!details_col %in% names(seurat_obj@misc) | reduction == "rna.umap") {
      details <- NULL
    } else {
      details <- seurat_obj@misc[[details_col]]
      details <- str_replace_all(details, "_", " ") # make it look nicer
    }

    if (use_adt) {
      for (comparison in comparisons) {
        plots_overview[[paste0(comparison, "_", obj_name)]] <-
          suppressMessages(
            FeaturePlot(seurat_obj, features = str_c("adt_", comparison),
                        order = TRUE, reduction = reduction, raster = FALSE) +
              labs(title = paste(assay_name, comparison), subtitle = details) +
              scale_color_viridis_c(name = comparison, option = "G",
                                    direction = -1) +
              labels_standard + clean_dimplot
          )
      }
    } else {
      # use colors from the seurat obj instead of named_colors because the
      # schemes (and therefore the cell type names) could vary by annotation method
      if (is.null(seurat_obj@misc$colors_annotated)) {
        cli::cli_inform("colors_annotated not found in seurat_obj@misc for {obj_name}. Using generated colors instead.")

        cell_types <- unique(seurat_obj[[]]$annotated_clusters)

        # if using CellTypist, use our color definitions
        if (all(cell_types %in% names(named_colors$cell_types_celltypist))) {
          seurat_obj@misc$colors_annotated <-
            named_colors$cell_types_celltypist
        } else {
          seurat_obj@misc$colors_annotated <-
            setNames(hues::iwanthue(length(cell_types)), cell_types)
        }
      }

      if ("annotated_clusters" %in% comparisons) {
        plots_overview[[paste0("annotated_clusters_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = seurat_obj@misc$colors_annotated,
                       plot_title = paste(assay_name, "Cell Types"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "annotated_clusters",
                       legend_label = "Cell Type", details = details, ...)
      }

      if ("annotated_clusters_simpler" %in% comparisons) {
        plots_overview[[paste0("annotated_clusters_simpler_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = seurat_obj@misc$colors_annotated,
                       # clrs_specific = named_colors$cell_types_celltypist,
                       plot_title = paste(assay_name, "Cell Types"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "annotated_clusters_simpler",
                       legend_label = "Cell Type", details = details, ...)
      }

      # V call families
      if ("v_call_family" %in% comparisons) {
        plots_overview[[paste0("v_call_family_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$v_call_family,
                       plot_title = paste(assay_name, "V Call Families"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "v_call_family",
                       legend_label = "V Call Family", details = details, ...)
      }

      # light chain types
      # TODO: switch to locus_light for consistency
      if ("light_chains" %in% comparisons) {
        plots_overview[[paste0("light_chains_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$locus_light,
                       plot_title = paste(assay_name, "Light Chain Types"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "locus_light",
                       legend_label = "Light Chain Type", details = details, ...)
      }

      # isotypes
      if ("isotype" %in% comparisons) {
        plots_overview[[paste0("isotype_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$isotype,
                       plot_title = paste(assay_name, "Isotypes"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "isotype",
                       legend_label = "Isotype", details = details, ...)
      }

      # subisotypes
      if ("c_call" %in% comparisons) {
        plots_overview[[paste0("c_call_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$c_call,
                       plot_title = paste(assay_name, "Subisotypes"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "c_call",
                       legend_label = "Subisotype", details = details, ...)
      }

      # SHM frequencies
      # TODO: switch this to actual mutation frequency
      if ("mu_freq" %in% comparisons) {
        plots_overview[[paste0("mu_freq_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$mu_freq_bins,
                       plot_title = paste(assay_name, "SHM Frequencies (Binned)"),
                       reduc = reduction,
                       meta_col = "mu_freq_bins", plot_label = FALSE,
                       legend_label = "Bins", sort_idents = FALSE,
                       order = TRUE, details = details, ...)
      }

      # CDR3 amino acid length
      if ("cdr3_aa_length" %in% comparisons) {
        plots_overview[[paste0("cdr3_aa_length_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$cdr3,
                       plot_title = paste(assay_name, "CDR3 Length"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "cdr3_aa_length",
                       legend_label = "CDR3 Length", idents_char = FALSE,
                       details = details, ...)
      }

      if ("weight_assay" %in% comparisons) {
        plots_overview[[paste0("weight_assay_", obj_name)]] <-
          plot_dimplot(seurat_obj = seurat_obj,
                       data_source = obj_data_source, pt_size = pt_size,
                       clrs_specific = named_colors$weight_assay,
                       plot_title = paste(assay_name, "Modality Weights"),
                       reduc = reduction, plot_label = FALSE,
                       meta_col = "weight_assay",
                       legend_label = "Chosen Assay", details = details, ...)
      }

      # catch other possible comparisons
      if (!all(comparisons %in% possible_comps)) {
        other_comps <- comparisons[!comparisons %in% possible_comps]
        cli::cli_alert("{length(other_comps)} comparison{?s} {?was/were} \\
                       not recognized and will be plotted if {?it/they} \\
                       {?is a/are} valid metadata column{?s}: {other_comps}")
        for (comp in other_comps) {
          plot_title <-
            stringr::str_to_title(stringr::str_replace_all(comp, "_", " "))

          if (comp %in% names(named_colors)) {
            plots_overview[[paste0(comp, "_", obj_name)]] <-
              plot_dimplot(seurat_obj = seurat_obj,
                           data_source = obj_data_source,
                           clrs_specific = named_colors[[comp]],
                           pt_size = pt_size,
                           plot_title = paste(assay_name, plot_title), reduc = reduction,
                           meta_col = comp, plot_label = FALSE,
                           details = details, ...)
          } else {
            plots_overview[[paste0(comp, "_", obj_name)]] <-
              plot_dimplot(seurat_obj = seurat_obj,
                           data_source = obj_data_source, use_hues = TRUE,
                           pt_size = pt_size,
                           plot_title = paste(assay_name, plot_title), reduc = reduction,
                           meta_col = comp, plot_label = FALSE,
                           details = details, ...)
          }
        }
      }
    }
  }

  # reorder plots to match the order of comparisons given
  plot_order <- unlist(lapply(comparisons, function(comp) {
                  paste0(comp, "_", names(seurat_objs))
                }))
  plot_order <- plot_order[plot_order %in% names(plots_overview)]
  plots_overview <- plots_overview[plot_order]

  # combine all of the plots
  if (!missing(ncol)) {
    wrap_plots(plots_overview, ncol = ncol, byrow = TRUE) +
      plot_anno + plot_layout(guides = "collect")
  } else if (length(seurat_objs) == 1) {
    wrap_plots(plots_overview, ncol = min(length(comparisons), 5)) + plot_anno
  } else {
    wrap_plots(plots_overview, nrow = length(comparisons), byrow = TRUE) +
      plot_anno + plot_layout(guides = "collect")
  }
}


#' Plot percentages in a stacked bar plot
#'
#' @description
#' This function plots a stacked bar plot of percentages calculated using [calc_pcts()] with percentages
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
#' @param data_source Dataset description.
#' @param clrs_specific The specific color palette (should be named).
#' @param plot_type One of `All`, `Binary`.
#' @param plot_value What is being plotted.
#' @param x_axis What to put along the x axis.
#' @param x_axis_label The label for the x axis.
#' @param fill_type What to group the bar plot by.
#' @param fill_label The description of what you're filling by.
#' @param perc_min The minimum percentage to show in the plot.
#' @param label_size The size of the percentage labels.
#' @param label_fill Add a white background for clarity.
#' @param label_counts Show raw counts instead of percentages as bar segment labels.
#' @param include_counts Plot the counts on top/bottom.
#' @param drop_zeroes Remove percentages of zeroes.
#' @param reverse_order Change the fill order for a stacked plot.
#' @param total_order Rearrange x axis in descending order by totals instead of alphabetically.
#' @param details An optional custom subtitle.
#'
#' @returns A stacked ggplot bar plot
#' @export
plot_pcts <- function(pcts, data_source = "", clrs_specific,
                      plot_type = "All", plot_value = "Cell Type",
                      x_axis = "sample_id", x_axis_label = "Sample",
                      fill_type = "annotated_clusters", fill_label = fill_type,
                      perc_min = 3, label_size = 3, label_fill = FALSE,
                      label_counts = FALSE, include_counts = TRUE, drop_zeroes = TRUE,
                      reverse_order = FALSE, total_order = FALSE,
                      horizontal = FALSE, details = NULL) {
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
    cli::cli_abort("Please enter a valid plot type.")
  }

  # set up the plot
  p <- ggplot(data = pcts,
              aes(x = !!sym(x_axis), y = Percent, fill = !!sym(fill_type))) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
        labs(title = paste(data_source, "Percent of",
                           plot_value, "per", x_axis_label),
             subtitle = details,
             x = x_axis_label, y = paste(plot_value, "Percentage"),
             fill = fill_label)

  # add aesthetics
  p <- p + scale_fill_manual(values = clrs_specific, limits = force) +
           scale_y_continuous(labels = scales::label_percent(scale = 1),
                              expand = expansion(mult = 0.05))

  # add labels
  bar_labels <- if (label_counts) {
    aes(label = ifelse(Percent > perc_min, Count, NA), group = !!sym(fill_type))
  } else {
    aes(label = scales::label_percent(accuracy = 1, scale = 1)(ifelse(Percent > perc_min, Percent, NA)),
        group = !!sym(fill_type))
  }
  if (!label_fill) {
    p <- p + geom_text(bar_labels, size = label_size, position = position_stack(vjust = 0.5))
  } else {
    p <- p + geom_label(bar_labels, color = "black", fill = "white", size = label_size,
                        position = position_stack(vjust = 0.5))
  }

  if (include_counts) {
    # actually add in the total counts
    if (plot_type == "All") {
      # TODO: rotate the text if horizontal?
      p <- p + geom_text(mapping = aes(x = !!sym(x_axis), y = 100,
                                       label = Total, fill = NULL),
                         hjust = if (horizontal) 0 else 0.5,
                         vjust = if (horizontal) -0.5 else -0.5,
                         position = if (horizontal) position_nudge(y = 0.2) else position_nudge(y = -0.2),
                         color = "black", size = label_size)
    } else if (plot_type == "Binary") { # just for TRUE v FALSE plots
      if (horizontal) {
        # place each count just past the right end of its segment
        p <- p +
          geom_text(mapping = aes(label = Count),
                    position = position_stack(), hjust = -0.2, vjust = 0.5,
                    color = "black", size = label_size)
      } else {
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
  }

  # add remaining styling
  p <- p + theme_bw_custom + labels_standard

  # remove gridlines parallel to bars
  if (horizontal) {
    p <- p + coord_flip()
    p <- p + theme(panel.grid.major.x = element_blank())
  } else {
    p <- p + theme(panel.grid.major.y = element_blank())
  }

  # fix the legend if needed
  if (reverse_order) {
    p <- p + guides(fill = guide_legend(reverse = TRUE))
  }

  # rotate (typically) long names
  if (!horizontal && x_axis %in% c("sample_id", "annotated_clusters",
                                   "annotated_subclusters")) {
    p <- p + labels_rotate_x
  }

  return(p)
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
#' @param clrs_specific The specific color palette (should be named).
#' @param feature The feature of interest.
#' @param reduc The reduction to use for plotting e.g. "bpca" or wnn.umap".
#' @param meta_col What to group by (uses the Idents by default).
#' @param rotate Rotate the labels or not.
#'
#' @returns Two patchworked Seurat plots.
#' @export
plot_vln_feat <- function(seurat_obj, clrs_specific, feature,
                          reduc = "umap", meta_col = NULL, rotate = FALSE, ...) {
  # TODO: pass additional parameters, make labelling optional (be careful with additional)
  # TODO: check if setting the assay is even necessary
  # TODO: deal with how clean_dimplot messing up the heights
  # TODO: give an option for how to sort the violin plot (or not)
  # TODO: include the assay to search for the feature in

  # set the assay and idents
  if (!is.null(meta_col)) Idents(seurat_obj) <- meta_col
  # DefaultAssay(seurat_obj) <- assay

  p1 <- VlnPlot(object = seurat_obj, features = feature, cols = clrs_specific,
                pt.size = 0.1, raster = FALSE) +
          NoLegend()

  if (rotate) p1 <- p1 + labels_standard_vln
  else p1 <- p1 + labels_standard_vln_rotate

  p2 <- FeaturePlot(object = seurat_obj, features = feature, pt.size = 0.1,
                    order = TRUE, min.cutoff = 0, reduction = reduc,
                    raster = FALSE, ...) +
          scale_color_viridis_c(option = "G", direction = -1) +
          labels_standard + # clean_dimplot
          theme(axis.ticks = element_blank(), axis.text = element_blank())

  (p1 | p2) & plot_layout(nrow = 1, widths = c(2, 1))
}
