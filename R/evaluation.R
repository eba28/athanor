#' Calculate the correlation between each cell's expression and the mean of
#' its neighbors' expression
#'
#' @details
#' The Seurat object must have `FindNeighbors()` already run and an ADT assay.
#' Correlations are calculated across all neighbor graphs present in the
#' object.
#'
#' @param seurat_obj The Seurat object.
#' @param nn_names Character vector of neighbor graph names to evaluate.
#'   Defaults to all graphs in the object.
#' @param features_adt Name of the ADT features to evaluate (e.g. "CD27.1").
#' @param adt_assay Name of the assay containing ADT data.
#' @param cor_method Correlation method to use (e.g. "pearson", "spearman").
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#' @param verbose Logical indicating whether or not to print messages.
#'
#' @returns A data frame with columns: Graph, Assay, Feature, Method, Score
#'   (the canonical schema shared by the `calc_*` evaluation functions).
#'   `Method` records the correlation method (e.g. "Spearman").
#' @export
calc_adt_correlation <- function(seurat_obj,
                                 nn_names = names(seurat_obj@neighbors),
                                 features_adt, adt_assay = "ADT",
                                 cor_method = "spearman", drop_self = TRUE,
                                 verbose = FALSE) {
  # TODO: only calculate on a subset of the neighbors if desired

  if (rlang::is_missing(features_adt)) {
    features_adt <- rownames(seurat_obj@assays[[adt_assay]])
  }

  if (any(!(features_adt %in% rownames(seurat_obj@assays[[adt_assay]])))) {
    available_features <- rownames(seurat_obj@assays[[adt_assay]])
    unavailable_features <- setdiff(features_adt, available_features)

    # only use features that are present in the object
    features_adt <- intersect(features_adt, available_features)

    # just for formatting
    available_features <- paste(sort(available_features), collapse = ", ")
    unavailable_features <- paste(sort(unavailable_features), collapse = ", ")

    # pass the pre-formatted strings to cli_inform
    if (verbose) {
      cli::cli_inform(c("i" = "Available features: {available_features}",
                        "i" = "Unavailable features: {unavailable_features}",
                        "v" = "Proceeding with: {features_adt}\n"))
    }
  }

  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  adt_data <- GetAssayData(seurat_obj, assay = adt_assay, layer = "data")
  metrics_df <- c()

  for (nn_name in nn_names) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx

    # a cell can appear as its own neighbor (e.g. distance-0 self match);
    # exclude it so the neighbor mean isn't partly computed from the cell
    # itself, which would inflate the correlation. drop_self_neighbors()
    # always drops one neighbor per cell (self if present, else the
    # farthest) so k stays fixed across cells regardless of whether self
    # was present for all, some, or none of them.
    if (drop_self) {
      if (verbose) {
        cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                        neighbor graph (self-match if present, else the \\
                        farthest neighbor); k reduced from {ncol(nn_idx)} to \\
                        {ncol(nn_idx) - 1}.")
      }
      nn_idx <- drop_self_neighbors(nn_idx)
    }

    neighbor_adt_mean <- vapply(seq_len(ncol(adt_data)), function(i) {
      rowMeans(adt_data[, nn_idx[i, ]])
    }, FUN.VALUE = numeric(nrow(adt_data)))
    colnames(neighbor_adt_mean) <- colnames(seurat_obj)

    for (feat in features_adt) {
      cell_expr <- adt_data[feat, ]
      neighbor_expr <- neighbor_adt_mean[feat, ]
      correlation <- cor(cell_expr, neighbor_expr, method = cor_method)
      metrics_df <- bind_rows(metrics_df,
                              data.frame(Graph = nn_name, Assay = adt_assay,
                                         Feature = feat,
                                         Method = str_to_title(cor_method),
                                         Score = correlation))
    }
  }

  if (verbose) {
    cli::cli_inform("Calculated {cor_method} correlation for the \\
                    {nn_names} neighbor graphs.")
    # cli::cli_inform("Calculated {cor_method} correlation between each cell's \\
    #                 ADT expression and the mean of its neighbors' expression for \\
    #                 {length(features_adt)} features across \\
    #                 {length(seurat_obj@neighbors)} neighbor graphs.")
  }

  metrics_df
}


#' Compute mean ADT distance to each cell's k nearest neighbors
#'
#' @description
#' Calculates the mean distance of each cell's ADT profile to its k nearest
#' neighbors in a specified neighbor space (e.g. RNA, BCR, WNN) using a chosen
#' distance metric (mean absolute difference, Manhattan, Euclidean, or cosine
#' similarity).
#'
#' @details
#' Euclidean distance measures straight-line distance between vectors.
#' Manhattan distance sums absolute differences across features.
#' Mean absolute difference is the per-feature mean of absolute differences.
#' Cosine similarity measures the cosine of the angle between two vectors,
#' with values closer to 1 indicating greater similarity.
#'
#' @param seurat_obj Seurat object that contains neighbor slots.
#' @param base_assay Which neighbor space to use (e.g. "RNA", "BCR", "WNN").
#' @param adt_assay ADT assay name (e.g. "ADT", "ADTnorm").
#' @param layer Data layer to pull from (data, counts, or scale.data).
#' @param features_adt Optional vector of ADT features; if NULL, use all
#'   features present.
#' @param k Number of nearest neighbors.
#' @param distance_metric One of `c("mean_abs", "manhattan", "euclidean",
#'   "cosine")`.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#'
#' @returns A data frame with columns: Graph, Assay, Feature, Method, Score
#'   (the canonical schema shared by the `calc_*` evaluation functions).
#'   Distances are computed jointly across all requested features, so
#'   `Feature` is the "+"-joined set of features used (or "all" when
#'   `features_adt` is NULL) rather than one row per feature. `Method` records
#'   the distance metric. If `return_mean = TRUE`, `Score` is the mean across
#'   all cells; else it contains per-cell values with an additional `cell_id`
#'   column.
#' @export
calc_adt_dists <- function(seurat_obj, base_assay, adt_assay = "ADT",
                           layer = "data", features_adt = NULL, k,
                           distance_metric = "mean_abs",
                           return_mean = TRUE, drop_self = TRUE) {
  adt_mat <- get_adt_matrix(seurat_obj, adt_assay, layer, features_adt)
  nn_idx <- resolve_neighbors(seurat_obj, base_assay)

  if (nrow(nn_idx) != nrow(adt_mat)) {
    cli::cli_abort("Neighbor index rows (", nrow(nn_idx), ") do not match the number ",
         "of cells (", nrow(adt_mat), "). Please ensure cell order matches ",
         "Cells(seurat_obj).")
  }

  scores <- adt_dists_core(adt_mat, nn_idx, distance_metric, drop_self,
                           return_mean = FALSE)

  method_labels <- c(mean_abs = "Mean_Absolute", manhattan = "Manhattan",
                     euclidean = "Euclidean", cosine = "Cosine")
  feature_label <- if (is.null(features_adt)) {
    "all"
  } else {
    paste(sort(colnames(adt_mat)), collapse = "+")
  }

  result <- data.frame(cell_id = rownames(adt_mat),
                       Graph = paste0(map_assay_name(base_assay), ".nn"),
                       Assay = adt_assay, Feature = feature_label,
                       Method = method_labels[[distance_metric]],
                       Score = unname(scores)) %>%
    remove_rownames()

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  result
}


#' Calculate mean absolute ADT distance to nearest neighbors
#'
#' @description
#' For each cell, calculates the mean absolute difference between its ADT
#' expression and that of its nearest neighbors across specified features.
#'
#' @details
#' Use `permute_adt()` to compute a permuted baseline. For
#' the range-based score, see `calc_adt_nn_within_range()`.
#'
#' @param seurat_obj A Seurat object.
#' @param nn_names Character vector of neighbor graph names to evaluate.
#'   Defaults to all graphs in the object.
#' @param features_adt Character vector of ADT feature names to evaluate.
#' @param adt_assay Name of the ADT assay.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame with columns: Graph, Assay, Feature, Method, Score.
#'   If `return_mean = TRUE`, `Score` is the mean across all cells; else it
#'   contains per-cell values with an additional `cell_id` column.
#' @export
calc_adt_mean_absolute <- function(seurat_obj,
                                   nn_names = names(seurat_obj@neighbors),
                                   features_adt, adt_assay = "ADT",
                                   return_mean = TRUE, drop_self = TRUE,
                                   verbose = FALSE) {
  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  invalid_nn <- setdiff(nn_names, names(seurat_obj@neighbors))
  if (length(invalid_nn) > 0) {
    cli::cli_abort("Neighbor graph(s) not found in object: {invalid_nn}")
  }

  features_adt <- intersect(features_adt, rownames(seurat_obj@assays[[adt_assay]]))

  all_graph_results <- vector("list", length(nn_names))
  names(all_graph_results) <- nn_names

  for (nn_name in nn_names) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx

    # a cell can appear as its own neighbor (e.g. distance-0 self match);
    # drop one neighbor per cell (self if present, else the farthest) so k
    # stays fixed across cells regardless of whether self was present.
    if (drop_self) {
      if (verbose) {
        cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                        neighbor graph (self-match if present, else the \\
                        farthest neighbor); k reduced from {ncol(nn_idx)} to \\
                        {ncol(nn_idx) - 1}.")
      }
      nn_idx <- drop_self_neighbors(nn_idx)
    }
    results_list <- vector("list", length(features_adt))

    for (i in seq_along(features_adt)) {
      feat <- features_adt[[i]]
      adt_mat <- get_adt_matrix(seurat_obj, adt_assay, "data", feat)
      scores_adt <- adt_dists_core(adt_mat, nn_idx,
                                   distance_metric = "mean_abs",
                                   drop_self = FALSE, return_mean = FALSE)

      results_list[[i]] <-
        data.frame(cell_id = seurat_obj$cell_id, Graph = nn_name,
                   Assay = adt_assay, Feature = feat,
                   Method = "Mean_Absolute", Score = unname(scores_adt))
    }

    all_graph_results[[nn_name]] <- do.call(rbind, results_list)
  }

  result <- do.call(rbind, all_graph_results)

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  if (verbose) {
    cli::cli_inform("Calculated ADT mean absolute difference scores for {nn_names}.")
  }

  result <- result %>% remove_rownames()

  result
}


#' Calculate Moran's i
#'
#' @description
#' Calculates the global Moran's i index for an ADT feature given a neighbor
#' graph.
#'
#' @details
#' Row standardization makes sure the resulting score will always be between
#' -1 and 1.
#'
#' @param seurat_obj The Seurat object. Must have `FindNeighbors()` already
#'   run.
#' @param features_adt Name of the ADT feature to evaluate (e.g. "CD27.1").
#' @param adt_assay Name of the assay containing ADT data.
#' @param graph_name Name of the neighbor graph slot to use for the weights
#'   matrix (e.g. "RNA.nn", "BCR.nn", "w.nn").
#' @param row_standardize Whether to row-standardize the weights matrix (i.e.
#'   make each row sum to 1).
#' @param verbose If TRUE, print progress messages.
#'
#' @returns A one-row data frame with columns: Graph, Assay, Feature, Method,
#'   Score (the canonical schema shared by the `calc_*` evaluation functions),
#'   where `Score` is the observed Moran's i and `Method` is "Moran_I".
#' @export
calc_adt_moran <- function(seurat_obj, features_adt, adt_assay = "ADT",
                           graph_name, row_standardize = TRUE,
                           verbose = FALSE) {
  if (!features_adt %in% rownames(seurat_obj@assays[[adt_assay]])) {
    cli::cli_abort("The requested ADT feature is not present in assay '", adt_assay,
         "'. Available features: ",
         paste(sort(rownames(seurat_obj@assays[[adt_assay]])), collapse = ", "))
  }

  if (!graph_name %in% names(seurat_obj@graphs)) {
    cli::cli_abort("Graph '", graph_name, "' not found in object. Available graphs: ",
         paste(names(seurat_obj@graphs), collapse = ", "))
  }

  x <- GetAssayData(seurat_obj, assay = adt_assay, layer = "data")
  x <- x[features_adt, ]

  w <- as.matrix(seurat_obj@graphs[[graph_name]])
  if (row_standardize) w <- w / rowSums(w)

  if (verbose) {
    cli::cli_inform("Calculating Moran's I for feature '{features_adt}' using \\
                    graph '{graph_name}' (row_standardize = {row_standardize}).")
  }

  moran <- ape::Moran.I(x, w)

  data.frame(Graph = graph_name, Assay = adt_assay, Feature = features_adt,
             Method = "Moran_I", Score = moran[["observed"]])
}


#' Calculate the proportion of neighbors within an ADT marker's expression
#' range
#'
#' @description
#' For each cell, calculates the proportion of its nearest neighbors whose ADT
#' expression falls within a symmetric relative threshold of the cell's own
#' expression. For example, with `range = 0.20`, neighbors within +/-20% of
#' the cell's expression are counted.
#'
#' @param seurat_obj A Seurat object.
#' @param nn_names Character vector of neighbor graph names to evaluate.
#'   Defaults to all graphs in the object.
#' @param features_adt Character vector of ADT feature names to evaluate.
#' @param adt_assay Name of the ADT assay.
#' @param range Relative threshold. A value of 0.20 means neighbors within
#'   +/-20% of the cell's expression are counted as matches.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame with columns: Graph, Assay, Feature, Method, Score.
#'   If `return_mean = TRUE`, `Score` is the mean across all cells; else it
#'   contains per-cell values with an additional `cell_id` column.
#' @export
calc_adt_nn_within_range <- function(seurat_obj,
                                     nn_names = names(seurat_obj@neighbors),
                                     features_adt,
                                     adt_assay = "ADT",
                                     range = 0.20,
                                     return_mean = TRUE, drop_self = TRUE,
                                     verbose = FALSE) {
  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  invalid_nn <- setdiff(nn_names, names(seurat_obj@neighbors))
  if (length(invalid_nn) > 0) {
    cli::cli_abort("Neighbor graph(s) not found in object: {invalid_nn}")
  }

  features_adt <- intersect(features_adt, rownames(seurat_obj@assays[[adt_assay]]))

  all_graph_results <- vector("list", length(nn_names))
  names(all_graph_results) <- nn_names

  for (nn_name in nn_names) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx

    # a cell can appear as its own neighbor (e.g. distance-0 self match);
    # exclude it so a cell isn't counted as an in-range neighbor of itself.
    # drop_self_neighbors() always drops one neighbor per cell (self if
    # present, else the farthest) so k stays fixed across cells regardless
    # of whether self was present for all, some, or none of them.
    if (drop_self) {
      if (verbose) {
        cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                        neighbor graph (self-match if present, else the \\
                        farthest neighbor); k reduced from {ncol(nn_idx)} to \\
                        {ncol(nn_idx) - 1}.")
      }
      nn_idx <- drop_self_neighbors(nn_idx)
    }
    n_cells <- nrow(nn_idx)
    k <- ncol(nn_idx)
    results_list <- vector("list", length(features_adt))

    for (i in seq_along(features_adt)) {
      feat <- features_adt[[i]]
      adt_expr <- GetAssayData(seurat_obj, assay = adt_assay, layer = "data")[feat, ]

      scores <- vapply(seq_len(n_cells), function(j) {
        cell_expr <- adt_expr[[j]]
        neighbor_expr <- adt_expr[nn_idx[j, ]]
        sum(neighbor_expr >= cell_expr * (1 - range) &
              neighbor_expr <= cell_expr * (1 + range)) / k
      }, FUN.VALUE = numeric(1))

      results_list[[i]] <-
        data.frame(cell_id = seurat_obj$cell_id, Graph = nn_name,
                   Assay = adt_assay, Feature = feat, Method = "Range",
                   Score = scores)
    }

    all_graph_results[[nn_name]] <- do.call(rbind, results_list)
  }

  result <- do.call(rbind, all_graph_results)

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  if (verbose) {
    cli::cli_inform("Calculated ADT range scores for {nn_names}.")
  }

  result
}


#' Calculate the proportion of neighbors within an ADT marker's quantile by
#' expression
#'
#' @description
#' For each cell in a Seurat object, calculates how many of its k nearest
#' neighbors have ADT expression within the same quantile as the cell's own
#' ADT expression for a given feature.
#' Note that the output of this function is a vector and not a data.frame like most of the other evaluation functions.
#'
#' @param seurat_obj A Seurat object containing ADT data and computed neighbor
#'   graphs.
#' @param adt_assay Name of the assay containing ADT data.
#' @param features_adt Name(s) of the ADT feature(s) to evaluate (e.g.
#'   "CD27.1"); one row group per feature is returned.
#' @param base_assay The assay used to compute neighbors. One of
#'   "RNA", "GEX", "BCR", or "WNN".
#' @param k Number of nearest neighbors to evaluate. Must match the k used
#'   when computing the neighbor graph.
#' @param n_quantile The number of bins to divide ADT expression into.
#'   Neighbors count as "within range" when they fall in the same bin as the
#'   cell.
#' @param method One of `c("quantile", "percentile_diff")`. `"quantile"`
#'   compares discrete bins. `"percentile_diff"` returns the mean absolute
#'   percentile rank difference per cell.
#' @param binning How the `n_quantile` bins are defined for
#'   `method = "quantile"`. `"frequency"` (the default) makes bins holding
#'   roughly equal numbers of cells (true quantiles, via
#'   [ggplot2::cut_number()]); `"width"` makes bins of equal expression range
#'   (via [ggplot2::cut_interval()]). Ignored for `method = "percentile_diff"`.
#' @param return_counts If TRUE, `Score` is the count of neighbors in the same
#'   bin (only for `method = "quantile"`); if FALSE, the proportion
#'   (count / k). Ignored for `method = "percentile_diff"`.
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. When TRUE,
#'   `k` is reduced by one to match. See [drop_self_neighbors()] for details.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame with columns: Graph, Assay, Feature, Method, Score (the
#'   canonical schema shared by the `calc_*` evaluation functions). `Method`
#'   is "Quantile" or "Percentile_Diff". If `return_mean = TRUE`, `Score` is
#'   the mean across all cells; else it contains per-cell values with an
#'   additional `cell_id` column.
#' @export
calc_adt_quantile <- function(seurat_obj, adt_assay = "ADT", features_adt,
                              base_assay, k = 20, n_quantile = 10,
                              method = c("quantile", "percentile_diff"),
                              binning = c("frequency", "width"),
                              drop_self = TRUE, return_counts = FALSE,
                              return_mean = TRUE, verbose = FALSE) {
  method <- match.arg(method)
  binning <- match.arg(binning)

  if (any(!(features_adt %in% rownames(seurat_obj@assays[[adt_assay]])))) {
    available_features <- rownames(seurat_obj@assays[[adt_assay]])
    unavailable_features <- setdiff(features_adt, available_features)

    # only use features that are present in the object
    features_adt <- intersect(features_adt, available_features)

    if (verbose) {
      available_features <- paste(sort(available_features), collapse = ", ")
      unavailable_features <- paste(sort(unavailable_features), collapse = ", ")
      cli::cli_inform(c(
        "i" = "Available features: {available_features}",
        "i" = "Unavailable features: {unavailable_features}",
        "v" = "Proceeding with: {features_adt}\n"
      ))
    }
  }

  if (return_counts && method == "percentile_diff") {
    warning("return_counts = TRUE is not applicable for ",
            "method = 'percentile_diff'. Returning mean absolute percentile ",
            "rank difference instead.")
  }

  if (rlang::is_missing(base_assay)) base_assay <- DefaultAssay(seurat_obj)
  nn_idx <- resolve_neighbors(seurat_obj, base_assay)
  nn_name <- paste0(map_assay_name(base_assay), ".nn")

  # a cell can appear as its own neighbor (e.g. distance-0 self match);
  # exclude it so a cell isn't counted as sharing its own bin.
  # drop_self_neighbors() always drops one neighbor per cell (self if
  # present, else the farthest), so reduce k to match the new width.
  if (drop_self) {
    if (verbose) {
      cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                      neighbor graph; k reduced from {ncol(nn_idx)} to \\
                      {ncol(nn_idx) - 1}.")
    }
    nn_idx <- drop_self_neighbors(nn_idx)
    k <- ncol(nn_idx)
  }

  adt_data <- GetAssayData(seurat_obj, assay = adt_assay, layer = "data")
  method_label <- if (method == "quantile") "Quantile" else "Percentile_Diff"

  results_list <- vector("list", length(features_adt))
  for (fi in seq_along(features_adt)) {
    feat <- features_adt[[fi]]
    adt_expr <- as.numeric(adt_data[feat, ])

    if (method == "quantile") {
      # "frequency" -> equal-count bins (true quantiles); "width" -> equal-range
      bin_fun <- if (binning == "frequency") cut_number else cut_interval
      bins <- as.numeric(bin_fun(adt_expr, n = n_quantile))

      scores <- vapply(seq_len(nrow(nn_idx)), function(i) {
        n_same <- sum(bins[nn_idx[i, ]] == bins[i], na.rm = TRUE)
        if (return_counts) n_same else n_same / k
      }, FUN.VALUE = numeric(1))
    } else {
      # percentile_diff: empirical percentile rank in [0, 1]
      adt_percentiles <- rep(NA_real_, length(adt_expr))
      valid <- !is.na(adt_expr)
      n_valid <- sum(valid)
      if (n_valid <= 1) {
        adt_percentiles[valid] <- 0
      } else {
        ranked <- rank(adt_expr[valid], ties.method = "average")
        adt_percentiles[valid] <- (ranked - 1) / (n_valid - 1)
      }

      scores <- vapply(seq_len(nrow(nn_idx)), function(i) {
        cell_pct <- adt_percentiles[i]
        if (is.na(cell_pct)) return(NA_real_)
        mean(abs(adt_percentiles[nn_idx[i, ]] - cell_pct), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
    }

    results_list[[fi]] <-
      data.frame(cell_id = seurat_obj$cell_id, Graph = nn_name,
                 Assay = adt_assay, Feature = feat, Method = method_label,
                 Score = scores)
  }

  result <- do.call(rbind, results_list) %>% remove_rownames()

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  if (verbose) {
    cli::cli_inform("Calculated ADT {method_label} scores for '{nn_name}'.")
  }

  result
}


#' Calculate neighbor matching scores across metadata columns
#'
#' @description
#' For each cell, calculates the proportion of its nearest neighbors that share
#' the same value for each metadata column.
#'
#' @details
#' `cell_id`s are saved for subsetting later if desired. Use
#' `permute_neighbor_matches()` to compute a permuted baseline.
#'
#' @param seurat_obj A Seurat object.
#' @param nn_names Character vector of neighbor graph names to evaluate.
#'   Defaults to all graphs in the object.
#' @param meta_cols Character vector of metadata columns to evaluate.
#' @param cdr3_length_range Integer range for `cdr3_aa_length` matching.
#'   Neighbors within this many amino acids of the query cell are counted as
#'   matches.
#' @param clonal_cols Character vector of `meta_cols` (typically
#'   `"clone_id_unique"`) that should be scored as clone-size-weighted matches
#'   instead of the usual proportion-of-k. For these columns, the score is the
#'   number of matching neighbors divided by `min(clone_size - 1, k)`, i.e.
#'   normalized by the number of same-clone cells actually available to be
#'   found (capped at k) rather than by k itself. This avoids penalizing
#'   small clones that can never fill all k neighbor slots. Cells whose clone
#'   is smaller than `min_clone_size` get a score of NA.
#' @param min_clone_size Minimum clone size (inclusive) for a cell to receive
#'   a score on a `clonal_cols` column; smaller clones get NA. Ignored for
#'   columns not listed in `clonal_cols`.
#' @param gex_cols Character vector of `meta_cols` that describe GEX-derived
#'   groupings (e.g. transcriptomic clusters). These get `Assay = "GEX"` in
#'   the output; every other column gets `Assay = "BCR"`. Defaults to the
#'   known GEX cluster columns plus `"seurat_clusters"`.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame with columns: Graph, Assay, Feature, Method, Score.
#'   If `return_mean = TRUE`, `Score` is the mean across all cells; else it
#'   contains per-cell values with an additional `cell_id` column.
#' @export
calc_neighbor_matches <- function(seurat_obj,
                                  nn_names = names(seurat_obj@neighbors),
                                  meta_cols =
                                    c("annotated_clusters_bcr",
                                      "annotated_clusters_binary",
                                      "annotated_clusters_gex_bcr",
                                      "annotated_clusters_simpler",
                                      "cdr3_aa_length", "clone_id_unique",
                                      "isotype_stage", "locus_light",
                                      "mu_freq_bins_binary", "v_call_family"),
                                  cdr3_length_range = 0,
                                  clonal_cols = character(0),
                                  min_clone_size = 2,
                                  gex_cols = c("annotated_clusters_binary",
                                               "annotated_clusters_gex_bcr",
                                               "annotated_clusters_simpler",
                                               "seurat_clusters"),
                                  return_mean = TRUE, drop_self = TRUE,
                                  verbose = FALSE) {
  # input validation
  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  # only use metadata columns available in the object
  meta_cols <- intersect(meta_cols, colnames(seurat_obj[[]]))

  invalid_nn <- setdiff(nn_names, names(seurat_obj@neighbors))
  if (length(invalid_nn) > 0) {
    cli::cli_abort("Neighbor graph(s) not found in object: {invalid_nn}")
  }

  all_graph_results <- vector("list", length(nn_names))
  names(all_graph_results) <- nn_names

  for (nn_name in nn_names) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx

    # a cell can appear as its own neighbor (e.g. distance-0 self match);
    # exclude it so it doesn't trivially "match" itself and inflate the
    # mixing score. drop_self_neighbors() always drops one neighbor per
    # cell (self if present, else the farthest) so k stays fixed across
    # cells regardless of whether self was present for all, some, or none
    # of them.
    if (drop_self) {
      if (verbose) {
        cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                        neighbor graph (self-match if present, else the \\
                        farthest neighbor); k reduced from {ncol(nn_idx)} to \\
                        {ncol(nn_idx) - 1}.")
      }
      nn_idx <- drop_self_neighbors(nn_idx)
    }
    n_cells <- nrow(nn_idx)
    k <- ncol(nn_idx)
    results_list <- vector("list", length(meta_cols))

    for (i in seq_along(meta_cols)) {
      meta_col <- meta_cols[[i]]
      meta_group <- seurat_obj@meta.data[[meta_col]]
      assay <- if (meta_col %in% gex_cols) "GEX" else "BCR"

      is_length_col <- meta_col == "cdr3_aa_length"
      is_clonal_col <- meta_col %in% clonal_cols

      # clone size normalizes by how many same-clone cells could possibly be
      # found (capped at k), rather than by k itself, so small clones aren't
      # penalized for having fewer same-clone neighbors to find
      if (is_clonal_col) {
        clone_size <- as.integer(table(meta_group)[as.character(meta_group)])
      }

      scores_mixing <- vapply(seq_len(n_cells), function(j) {
        if (is_length_col) {
          n_matches <- sum(abs(meta_group[nn_idx[j, ]] - meta_group[j]) <=
                             cdr3_length_range, na.rm = TRUE)
        } else {
          n_matches <- sum(meta_group[nn_idx[j, ]] == meta_group[j], na.rm = TRUE)
        }
        denom <- if (is_clonal_col) min(clone_size[j] - 1, k) else k
        score <- n_matches / denom
        if (is.na(score)) NA else score
      }, FUN.VALUE = numeric(1))

      if (is_clonal_col) {
        scores_mixing[is.na(meta_group) | clone_size < min_clone_size] <- NA
      }

      if (any(scores_mixing > 1, na.rm = TRUE)) {
        cli::cli_abort("Scores are too high, please check the \\
                       calculation (the detected k is {k})")
      }

      method <- if (is_length_col) paste0("Tol", cdr3_length_range) else "Exact"
      valid_cells <- !is.na(scores_mixing)
      results_list[[i]] <-
        data.frame(cell_id = seurat_obj$cell_id[valid_cells], Graph = nn_name,
                   Assay = assay, Feature = meta_col, Method = method,
                   Score = scores_mixing[valid_cells])
    }

    all_graph_results[[nn_name]] <- do.call(rbind, results_list)
  }

  result <- do.call(rbind, all_graph_results)
  result <- result %>% remove_rownames()

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  if (verbose) {
    cli::cli_inform("Calculated neighbor matching scores for {nn_names}.")
  }

  result
}


#' Calculate shared neighbors between two neighbor graphs
#'
#' @description
#' For each cell in a Seurat object, counts how many of its nearest neighbors
#' are shared between two specified neighbor graphs.
#'
#' @param seurat_obj A Seurat object.
#' @param nn1 Name of the first neighbor graph (e.g. "RNA.nn").
#' @param nn2 Name of the second neighbor graph (e.g. "BCR.nn").
#' @param return_mean If TRUE, return the mean shared neighbor count across all
#'   cells; else return a named numeric vector of per-cell counts.
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) from each graph
#'   so `k` is fixed across cells regardless of whether self-matches are
#'   present. See [drop_self_neighbors()] for details.
#' @param verbose If TRUE, print a summary message.
#'
#' @return Data frame with columns: Graph, Assay, Feature, Method, Score (the
#'   canonical schema shared by the `calc_*` evaluation functions). `Graph` is
#'   `"{nn1} vs {nn2}"`, `Assay` and `Feature` are NA (no ADT feature is
#'   involved), and `Method` is "Overlap". If `return_mean = TRUE`, `Score` is
#'   the mean shared-neighbor count across all cells; else it contains
#'   per-cell counts with an additional `cell_id` column.
#' @export
calc_neighbor_overlap <- function(seurat_obj, nn1, nn2, drop_self = TRUE,
                                  return_mean = TRUE, verbose = FALSE) {
  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  missing_nn <- setdiff(c(nn1, nn2), names(seurat_obj@neighbors))
  if (length(missing_nn) > 0) {
    cli::cli_abort("Neighbor graph(s) not found in object: {missing_nn}")
  }

  # get who the neighbors are
  nn_idx1 <- seurat_obj@neighbors[[nn1]]@nn.idx
  nn_idx2 <- seurat_obj@neighbors[[nn2]]@nn.idx

  # they should always be in the same order, but check just in case
  cells1 <- seurat_obj@neighbors[[nn1]]@cell.names
  cells2 <- seurat_obj@neighbors[[nn2]]@cell.names
  if (!identical(cells1, cells2)) {
    cli::cli_abort("Cell order differs between {nn1} and {nn2}. \
                   Ensure both neighbor graphs were computed on the same object \
                   without reordering cells.")
  }

  # this should never happen
  if (nrow(nn_idx1) != nrow(nn_idx2)) {
    cli::cli_abort("Neighbor graphs have different numbers of cells \\
                   ({nrow(nn_idx1)} vs {nrow(nn_idx2)}).")
  }

  # a cell can appear as its own neighbor (e.g. distance-0 self match);
  # drop one neighbor per cell (self if present, else the farthest) from
  # each graph so self-matches don't inflate the shared count and k stays
  # fixed across cells regardless of whether self was present.
  if (drop_self) {
    if (verbose) {
      cli::cli_inform("Dropping one neighbor per cell from {nn1} and {nn2} \\
                      (self-match if present, else the farthest neighbor).")
    }
    nn_idx1 <- drop_self_neighbors(nn_idx1)
    nn_idx2 <- drop_self_neighbors(nn_idx2)
  }

  shared_counts <- vapply(seq_len(ncol(seurat_obj)), function(i) {
    length(intersect(nn_idx1[i, ], nn_idx2[i, ]))
  }, FUN.VALUE = numeric(1))

  if (verbose) {
    cli::cli_inform("Calculated shared neighbor counts between {nn1} and {nn2}.")
  }

  result <- data.frame(cell_id = seurat_obj$cell_id,
                       Graph = paste(nn1, "vs", nn2),
                       Assay = NA_character_, Feature = NA_character_,
                       Method = "Overlap", Score = shared_counts) %>%
    remove_rownames()

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  result
}


#' Compute a permuted baseline for ADT neighbor matching scores
#'
#' @description
#' Runs ADT neighbor matching `n_permutations` times with shuffled expression
#' values and returns the mean score across permutations as a random baseline.
#' Supports both the mean absolute and range methods.
#'
#' @param seurat_obj A Seurat object.
#' @param nn_names Character vector of neighbor graph names to evaluate.
#'   Defaults to all graphs in the object.
#' @param features_adt Character vector of ADT feature names to evaluate.
#' @param adt_assay Name of the ADT assay.
#' @param adt_range Relative threshold for the range method.
#' @param methods Character vector of methods. One or both of `"mean_abs"` and
#'   `"range"`.
#' @param n_permutations Number of permutations to run.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values (averaged across permutations).
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame in the same format as `calc_adt_mean_absolute()` and
#'   `calc_adt_nn_within_range()`.
#' @export
permute_adt <- function(seurat_obj, nn_names = names(seurat_obj@neighbors),
                        features_adt, adt_assay = "ADT", adt_range = 0.1,
                        methods = c("mean_abs", "range"),
                        n_permutations = 10, return_mean = TRUE,
                        drop_self = TRUE, verbose = FALSE) {
  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  invalid_nn <- setdiff(nn_names, names(seurat_obj@neighbors))
  if (length(invalid_nn) > 0) {
    cli::cli_abort("Neighbor graph(s) not found in object: {invalid_nn}")
  }

  features_adt <- intersect(features_adt, rownames(seurat_obj@assays[[adt_assay]]))

  all_graph_results <- vector("list", length(nn_names))
  names(all_graph_results) <- nn_names

  for (nn_name in nn_names) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx

    # a cell can appear as its own neighbor (e.g. distance-0 self match);
    # drop one neighbor per cell (self if present, else the farthest) so
    # k stays fixed across cells regardless of whether self was present.
    if (drop_self) {
      if (verbose) {
        cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                        neighbor graph (self-match if present, else the \\
                        farthest neighbor); k reduced from {ncol(nn_idx)} to \\
                        {ncol(nn_idx) - 1}.")
      }
      nn_idx <- drop_self_neighbors(nn_idx)
    }
    n_cells <- nrow(nn_idx)
    k <- ncol(nn_idx)

    perm_results <- vector("list", n_permutations)

    for (perm_idx in seq_len(n_permutations)) {
      results_list <- vector("list", length(features_adt) * length(methods))
      i <- 1

      for (feat in features_adt) {
        adt_expr <- as.numeric(GetAssayData(seurat_obj, assay = adt_assay, layer = "data")[feat, ])
        adt_expr <- adt_expr[sample(n_cells)]

        if ("mean_abs" %in% methods) {
          scores_adt <- vapply(seq_len(n_cells), function(j) {
            idx <- nn_idx[j, ]
            if (length(idx) == 0) return(NA)
            mean(abs(adt_expr[idx] - adt_expr[j]))
          }, FUN.VALUE = numeric(1))

          results_list[[i]] <-
            data.frame(cell_id = seurat_obj$cell_id, Graph = nn_name,
                       Assay = adt_assay, Feature = feat,
                       Method = "Mean_Absolute", Score = scores_adt)
          i <- i + 1
        }

        if ("range" %in% methods) {
          scores_adt <- vapply(seq_len(n_cells), function(j) {
            idx <- nn_idx[j, ]
            if (length(idx) == 0) return(NA)
            sum(adt_expr[idx] >= adt_expr[j] * (1 - adt_range) &
                  adt_expr[idx] <= adt_expr[j] * (1 + adt_range)) / k
          }, FUN.VALUE = numeric(1))

          results_list[[i]] <-
            data.frame(cell_id = seurat_obj$cell_id, Graph = nn_name,
                       Assay = adt_assay, Feature = feat, Method = "Range",
                       Score = scores_adt)
          i <- i + 1
        }
      }

      perm_results[[perm_idx]] <- do.call(rbind, results_list[seq_len(i - 1)])
    }

    all_graph_results[[nn_name]] <- do.call(rbind, perm_results)
  }

  result <- do.call(rbind, all_graph_results) %>%
    group_by(cell_id, Graph, Assay, Feature, Method) %>%
    summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  if (verbose) {
    cli::cli_inform("Computed permuted baseline ({n_permutations} permutations) \\
                    for {nn_names}.")
  }

  result
}


#' Compute a permuted baseline for neighbor matching scores
#'
#' @description
#' Runs `calc_neighbor_matches()` `n_permutations` times with shuffled metadata
#' labels and returns the mean score across permutations as a random baseline.
#'
#' @param seurat_obj A Seurat object.
#' @param nn_names Character vector of neighbor graph names to evaluate.
#'   Defaults to all graphs in the object.
#' @param meta_cols Character vector of metadata columns to evaluate.
#' @param gex_cols Character vector of `meta_cols` that describe GEX-derived
#'   groupings; these get `Assay = "GEX"` in the output and every other column
#'   gets `Assay = "BCR"`. Defaults to the known GEX cluster columns plus
#'   `"seurat_clusters"`.
#' @param n_permutations Number of permutations to run.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values (averaged across permutations).
#' @param drop_self Logical indicating whether to drop one neighbor per cell
#'   (the self-match if present, else the farthest neighbor) so `k` is fixed
#'   across cells regardless of whether self-matches are present. See
#'   [drop_self_neighbors()] for details.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame in the same format as `calc_neighbor_matches()`.
#' @export
permute_neighbor_matches <- function(seurat_obj,
                                     nn_names = names(seurat_obj@neighbors),
                                     meta_cols =
                                       c("annotated_clusters_bcr",
                                         "annotated_clusters_binary",
                                         "annotated_clusters_gex_bcr",
                                         "annotated_clusters_simpler",
                                         "cdr3_aa_length", "clone_id_unique",
                                         "isotype_stage", "locus_light",
                                         "mu_freq_bins_binary", "v_call_family"),
                                     gex_cols = c("annotated_clusters_binary",
                                                  "annotated_clusters_gex_bcr",
                                                  "annotated_clusters_simpler",
                                                  "seurat_clusters"),
                                     n_permutations = 10,
                                     return_mean = TRUE, drop_self = TRUE,
                                     verbose = FALSE) {
  if (length(seurat_obj@neighbors) == 0) {
    cli::cli_abort("No neighbor graphs found in object. \\
                   Please run `FindNeighbors()` first.")
  }

  invalid_nn <- setdiff(nn_names, names(seurat_obj@neighbors))
  if (length(invalid_nn) > 0) {
    cli::cli_abort("Neighbor graph(s) not found in object: {invalid_nn}")
  }

  all_graph_results <- vector("list", length(nn_names))
  names(all_graph_results) <- nn_names

  for (nn_name in nn_names) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx

    # a cell can appear as its own neighbor (e.g. distance-0 self match);
    # exclude it so it doesn't trivially "match" itself and inflate the
    # mixing score. drop_self_neighbors() always drops one neighbor per
    # cell (self if present, else the farthest) so k stays fixed across
    # cells regardless of whether self was present for all, some, or none
    # of them.
    if (drop_self) {
      if (verbose) {
        cli::cli_inform("Dropping one neighbor per cell from the '{nn_name}' \\
                        neighbor graph (self-match if present, else the \\
                        farthest neighbor); k reduced from {ncol(nn_idx)} to \\
                        {ncol(nn_idx) - 1}.")
      }
      nn_idx <- drop_self_neighbors(nn_idx)
    }
    n_cells <- nrow(nn_idx)
    k <- ncol(nn_idx)

    perm_results <- vector("list", n_permutations)

    for (perm_idx in seq_len(n_permutations)) {
      results_list <- vector("list", length(meta_cols))

      for (i in seq_along(meta_cols)) {
        meta_col <- meta_cols[[i]]
        meta_group <- seurat_obj@meta.data[[meta_col]][sample(n_cells)]
        assay <- if (meta_col %in% gex_cols) "GEX" else "BCR"

        scores_mixing <- vapply(seq_len(n_cells), function(j) {
          score <- sum(meta_group[nn_idx[j, ]] == meta_group[j], na.rm = TRUE) / k
          if (is.na(score)) NA else score
        }, FUN.VALUE = numeric(1))

        valid_cells <- !is.na(scores_mixing)
        results_list[[i]] <-
          data.frame(cell_id = seurat_obj$cell_id[valid_cells], Graph = nn_name,
                     Assay = assay, Feature = meta_col, Method = "Exact",
                     Score = scores_mixing[valid_cells])
      }

      perm_results[[perm_idx]] <- do.call(rbind, results_list)
    }

    all_graph_results[[nn_name]] <- do.call(rbind, perm_results)
  }

  result <- do.call(rbind, all_graph_results) %>%
    group_by(cell_id, Graph, Assay, Feature, Method) %>%
    summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")

  if (return_mean) {
    result <- result %>%
      group_by(Graph, Assay, Feature, Method) %>%
      summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  if (verbose) {
    cli::cli_inform("Computed permuted baseline ({n_permutations} permutations) \\
                    for {nn_names}.")
  }

  result
}
