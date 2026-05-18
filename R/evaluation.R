#' Calculate the correlation between each cell's expression and the mean of
#' its neighbors' expression
#'
#' @details
#' The Seurat object must have `FindNeighbors()` already run and an ADT assay.
#' Correlations are calculated across all neighbor graphs present in the
#' object.
#'
#' @param seurat_obj The Seurat object.
#' @param features_adt Name of the ADT features to evaluate (e.g. "CD27.1").
#' @param adt_assay Name of the assay containing ADT data.
#' @param cor_method Correlation method to use (e.g. "pearson", "spearman").
#' @param verbose Logical indicating whether or not to print messages.
#'
#' @returns A data frame with columns: Graph, Feature, Score.
#' @export
calc_adt_correlation <- function(seurat_obj, features_adt, adt_assay = "ADT",
                                 cor_method = "spearman", verbose = FALSE) {
  # TODO: only calculate on a subset of the neighbors if desired
  # TODO: print a list of graphs used
  # TODO: include the airrflow version?

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

  # TODO: throw an error if no RNA graph is found

  adt_data <- GetAssayData(seurat_obj, assay = adt_assay, layer = "data")
  metrics_df <- c()

  for (nn_name in names(seurat_obj@neighbors)) {
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx
    neighbor_adt_mean <- vapply(seq_len(ncol(adt_data)), function(i) {
      rowMeans(adt_data[, nn_idx[i, ]])
    }, FUN.VALUE = numeric(nrow(adt_data)))
    colnames(neighbor_adt_mean) <- colnames(seurat_obj)

    for (feat in features_adt) {
      cell_expr <- adt_data[feat, ]
      neighbor_expr <- neighbor_adt_mean[feat, ]
      correlation <- cor(cell_expr, neighbor_expr, method = cor_method)
      # TODO: include Method?
      metrics_df <- bind_rows(metrics_df,
                              data.frame(Graph = nn_name, Feature = feat,
                                         Score = correlation))
    }
  }

  if (verbose) {
    cli::cli_inform("Calculated {cor_method} correlation for the \\
                    {names(seurat_obj@neighbors)} neighbor graphs.")
    # cli::cli_inform("Calculated {cor_method} correlation between each cell's \\
    #                 ADT expression and the mean of its neighbors' expression for \\
    #                 {length(features_adt)} features across \\
    #                 {length(seurat_obj@neighbors)} neighbor graphs.")
  }

  # TODO: return with a column for Type = "Correlation" # spearman too
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
#' @param exclude_self Drop the cell itself from its neighbor list if present.
#'
#' @returns A single numeric value if `return_mean = TRUE`, or a named numeric
#'   vector of per-cell distances if `return_mean = FALSE`.
#' @export
calc_adt_dists <- function(seurat_obj, base_assay, adt_assay = "ADT",
                           layer = "data", features_adt = NULL, k,
                           distance_metric = "mean_abs",
                           return_mean = TRUE, exclude_self = TRUE) {
  adt_mat <- get_adt_matrix(seurat_obj, adt_assay, layer, features_adt)
  nn_idx <- resolve_neighbors(seurat_obj, base_assay)

  if (nrow(nn_idx) != nrow(adt_mat)) {
    cli::cli_abort("Neighbor index rows (", nrow(nn_idx), ") do not match the number ",
         "of cells (", nrow(adt_mat), "). Please ensure cell order matches ",
         "Cells(seurat_obj).")
  }

  adt_dists_core(adt_mat, nn_idx, distance_metric, exclude_self, return_mean)
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
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame with columns: Graph, Assay, Feature, Method, Score.
#'   If `return_mean = TRUE`, `Score` is the mean across all cells; else it
#'   contains per-cell values with an additional `cell_id` column.
#' @export
calc_adt_mean_absolute <- function(seurat_obj,
                                   nn_names = names(seurat_obj@neighbors),
                                   features_adt, adt_assay = "ADT",
                                   return_mean = TRUE, verbose = FALSE) {
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
    k <- ncol(nn_idx)
    assay_name <- substr(nn_name, 1, nchar(nn_name) - 3)
    results_list <- vector("list", length(features_adt))

    for (i in seq_along(features_adt)) {
      feat <- features_adt[[i]]
      scores_adt <-
        calc_adt_dists(seurat_obj = seurat_obj,
                       base_assay = assay_name, adt_assay = adt_assay,
                       features_adt = feat, k = k,
                       distance_metric = "mean_abs",
                       return_mean = FALSE)

      results_list[[i]] <-
        data.frame(cell_id = seurat_obj$cell_id, Graph = nn_name,
                   Assay = adt_assay, Feature = feat,
                   Method = "Mean_Absolute", Score = scores_adt)
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
#'
#' @returns A single numeric value representing the observed Moran's i.
#' @export
calc_adt_moran <- function(seurat_obj, features_adt, adt_assay = "ADT",
                           graph_name, row_standardize = TRUE) {
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

  moran <- ape::Moran.I(x, w)
  moran[["observed"]]
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
                                     return_mean = TRUE, verbose = FALSE) {
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
#'
#' @param seurat_obj A Seurat object containing ADT data and computed neighbor
#'   graphs.
#' @param adt_assay Name of the assay containing ADT data.
#' @param features_adt Name of the ADT feature to evaluate (e.g. "CD27.1").
#' @param base_assay The assay used to compute neighbors. One of
#'   "RNA", "GEX", "BCR", or "WNN".
#' @param k Number of nearest neighbors to evaluate. Must match the k used
#'   when computing the neighbor graph.
#' @param n_quantile The number of quantiles to divide ADT expression into.
#'   Neighbors are considered "within range" if they fall into the same
#'   quantile as the cell.
#' @param method One of `c("quantile", "percentile_diff")`. `"quantile"`
#'   compares discrete quantile bins. `"percentile_diff"` returns mean
#'   absolute percentile rank difference per cell.
#' @param return_counts If TRUE, returns the count of neighbors within the
#'   same quantile bin (only applies to `method = "quantile"`). If FALSE,
#'   returns the proportion (count/k). Ignored when
#'   `method = "percentile_diff"`.
#'
#' @return A named numeric vector with one value per cell in the Seurat object.
#'   For `method = "quantile"`: count or proportion of neighbors in the same
#'   quantile bin. For `method = "percentile_diff"`: mean absolute percentile
#'   rank difference per cell (0 to 1).
#' @export
calc_adt_quantile <- function(seurat_obj, adt_assay = "ADT", features_adt,
                              base_assay, k = 20, n_quantile = 10,
                              method = c("quantile", "percentile_diff"),
                              return_counts = FALSE) {
  # TODO: update this to match the style and output of adt_correlation
  method <- match.arg(method)

  if (any(!(features_adt %in% rownames(seurat_obj@assays[[adt_assay]])))) {
    available_features <- rownames(seurat_obj@assays[[adt_assay]])
    unavailable_features <- setdiff(features_adt, available_features)

    # only use features that are present in the object
    features_adt <- intersect(features_adt, available_features)

    # just for formatting
    available_features <- paste(sort(available_features), collapse = ", ")
    unavailable_features <- paste(sort(unavailable_features), collapse = ", ")

    # pass the pre-formatted strings to cli_inform
    cli::cli_inform(c(
      "i" = "Available features: {available_features}",
      "i" = "Unavailable features: {unavailable_features}",
      "v" = "Proceeding with: {features_adt}\n"
    ))
    # TODO: make having a print-out optional?
  }

  if (rlang::is_missing(base_assay)) base_assay <- DefaultAssay(seurat_obj)
  nn_idx <- resolve_neighbors(seurat_obj, base_assay)

  adt_expr <- GetAssayData(seurat_obj, assay = adt_assay, layer = "data")[features_adt, ]
  adt_expr <- as.numeric(adt_expr)

  if (method == "quantile") {
    adt_quantiles <- cut_interval(adt_expr, n = n_quantile)
    adt_quantiles <- as.numeric(adt_quantiles)

    results <- vapply(seq_len(nrow(nn_idx)), function(i) {
      cell_quantile <- adt_quantiles[i]
      neighbor_quantile <- adt_quantiles[nn_idx[i, ]]
      within_range <- sum(neighbor_quantile == cell_quantile, na.rm = TRUE)
      if (return_counts) within_range else within_range / k
    }, FUN.VALUE = numeric(1))

    names(results) <- Cells(seurat_obj)
    return(results)
  }

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

  if (return_counts) {
    warning("return_counts = TRUE is not applicable for ",
            "method = 'percentile_diff'. Returning mean absolute percentile ",
            "rank difference instead.")
  }

  results <- vapply(seq_len(nrow(nn_idx)), function(i) {
    cell_pct <- adt_percentiles[i]
    if (is.na(cell_pct)) return(NA_real_)
    neighbor_pct <- adt_percentiles[nn_idx[i, ]]
    mean(abs(neighbor_pct - cell_pct), na.rm = TRUE)
  }, FUN.VALUE = numeric(1))

  names(results) <- Cells(seurat_obj)
  results
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
#'   matches. Defaults to 1.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
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
                                  cdr3_length_range = 1,
                                  return_mean = TRUE, verbose = FALSE) {
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
    n_cells <- nrow(nn_idx)
    k <- ncol(nn_idx)
    results_list <- vector("list", length(meta_cols))

    for (i in seq_along(meta_cols)) {
      meta_col <- meta_cols[[i]]
      meta_group <- seurat_obj@meta.data[[meta_col]]
      assay <- ifelse(meta_col %in% str_c("annotated_clusters_",
                                          c("binary", "gex_bcr", "simpler")),
                      "GEX", "BCR")

      is_length_col <- meta_col == "cdr3_aa_length"
      scores_mixing <- vapply(seq_len(n_cells), function(j) {
        if (is_length_col) {
          score <- sum(abs(meta_group[nn_idx[j, ]] - meta_group[j]) <=
                         cdr3_length_range, na.rm = TRUE) / k
        } else {
          score <- sum(meta_group[nn_idx[j, ]] == meta_group[j], na.rm = TRUE) / k
        }
        if (is.na(score)) NA else score
      }, FUN.VALUE = numeric(1))

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
#' @param verbose If TRUE, print progress messages.
#'
#' @return Data frame in the same format as `calc_adt_mean_absolute()` and
#'   `calc_adt_nn_within_range()`.
#' @export
permute_adt <- function(seurat_obj, nn_names = names(seurat_obj@neighbors),
                        features_adt, adt_assay = "ADT", adt_range = 0.1,
                        methods = c("mean_abs", "range"),
                        n_permutations = 10, return_mean = TRUE,
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
            idx <- nn_idx[j, nn_idx[j, ] != j]
            if (length(idx) == 0) return(NA_real_)
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
            idx <- nn_idx[j, nn_idx[j, ] != j]
            if (length(idx) == 0) return(NA_real_)
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
#' @param n_permutations Number of permutations to run.
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values (averaged across permutations).
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
                                     n_permutations = 10,
                                     return_mean = TRUE, verbose = FALSE) {
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
    n_cells <- nrow(nn_idx)
    k <- ncol(nn_idx)

    perm_results <- vector("list", n_permutations)

    for (perm_idx in seq_len(n_permutations)) {
      results_list <- vector("list", length(meta_cols))

      for (i in seq_along(meta_cols)) {
        meta_col <- meta_cols[[i]]
        meta_group <- seurat_obj@meta.data[[meta_col]][sample(n_cells)]
        assay <- ifelse(meta_col %in% str_c("annotated_clusters_",
                                            c("binary", "gex_bcr", "simpler")),
                        "GEX", "BCR")

        scores_mixing <- vapply(seq_len(n_cells), function(j) {
          score <- sum(meta_group[nn_idx[j, ]] == meta_group[j], na.rm = TRUE) / k
          if (is.na(score)) NA_real_ else score
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
