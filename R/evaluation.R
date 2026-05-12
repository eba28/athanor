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
#' For each cell in a Seurat object, calculates how many of its k nearest
#' neighbors have ADT expression within a specified threshold (default 20%) of
#' the cell's own ADT expression for a given feature.
#'
#' @details
#' The range is symmetric around the cell's expression value. For example,
#' with range = 0.20:
#' \itemize{
#'   \item Lower bound = cell_expr * 0.80
#'   \item Upper bound = cell_expr * 1.20
#' }
#'
#' @param seurat_obj A Seurat object containing ADT data and computed neighbor
#'   graphs.
#' @param adt_assay Name of the assay containing ADT data.
#' @param features_adt Name of the ADT feature to evaluate (e.g. "CD27.1").
#' @param base_assay The assay used to compute neighbors. One of
#'   "RNA", "GEX", "BCR", or "WNN".
#' @param k Number of nearest neighbors to evaluate. Must match the k used
#'   when computing the neighbor graph.
#' @param range The relative threshold for considering neighbors similar. A
#'   value of 0.20 means neighbors within +/-20% of the cell's expression are
#'   counted.
#' @param return_counts If TRUE, returns the count of neighbors within range.
#'   If FALSE, returns the proportion (count/k).
#'
#' @return A named numeric vector with one value per cell in the Seurat object.
#' @export
calc_adt_nn_within_range <- function(seurat_obj, adt_assay = "ADT",
                                     features_adt, base_assay, k = 20,
                                     range = 0.20, return_counts = FALSE) {
  if (rlang::is_missing(base_assay)) base_assay <- DefaultAssay(seurat_obj)
  nn_idx <- resolve_neighbors(seurat_obj, base_assay)

  adt_expr <- seurat_obj@assays[[adt_assay]]@data[features_adt, ]

  results <- vapply(seq_len(nrow(nn_idx)), function(i) {
    cell_expr <- adt_expr[i]
    neighbor_expr <- adt_expr[nn_idx[i, ]]
    lower_bound <- cell_expr * (1 - range)
    upper_bound <- cell_expr * (1 + range)
    within_range <- sum(neighbor_expr >= lower_bound &
                          neighbor_expr <= upper_bound)
    if (return_counts) within_range else within_range / k
  }, FUN.VALUE = numeric(1))

  names(results) <- Cells(seurat_obj)
  results
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

  adt_expr <- seurat_obj@assays[[adt_assay]]@data[features_adt, ]
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
#' Calculates the proportion of neighbors that match each cell's metadata
#' category for specified metadata columns and ADT features. Can also perform
#' permutations to generate a random baseline. Results can be saved to a
#' specified path.
#'
#' @details
#' `cell_id`s are saved for subsetting later if desired.
#'
#' @param seurat_obj The Seurat object, with details added to the Misc() slot.
#' @param nn_name Name of the nearest neighbor graph slot.
#' @param meta_cols Character vector of metadata columns to evaluate.
#' @param adt_features Character vector of ADT feature names to evaluate.
#' @param adt_range Numeric vector of ADT expression range threshold(s).
#' @param adt_methods How to calculate the ADT metric(s).
#' @param permute Shuffle labels for each meta column and ADT expression per
#'   feature before computing matches to get a random baseline.
#' @param n_permutations The number of times to permute labels.
#' @param previous_matches Data frame of previous matches to combine with new
#'   results (optional).
#' @param return_mean If TRUE, return the mean across all cells; else return
#'   per-cell values.
#' @param path_save Where to save the results.
#'
#' @return Data frame with columns: Full_Name, Category, Category_Details,
#'   Assay, Meta_Col, Method, Matches. If `return_mean = TRUE`, `Matches` will be the mean across all cells; else it will contain per-cell values.
#' @export
calc_neighbor_matches <- function(seurat_obj, nn_name,
                                  meta_cols =
                                    c("annotated_clusters_bcr",
                                      "annotated_clusters_binary",
                                      "annotated_clusters_gex_bcr",
                                      "annotated_clusters_simpler",
                                      "cdr3_aa_length", "clone_id_unique",
                                      "isotype_stage", "locus_light",
                                      "mu_freq_bins_binary", "v_call_family"),
                                  adt_features = NULL, adt_range = 0.1,
                                  adt_methods = c("mean_abs", "range"),
                                  permute = FALSE, n_permutations = 10,
                                  previous_matches, return_mean = TRUE,
                                  path_save) {
  # TODO: improve the assay that is returned
  # TODO: get rid of using category and category_details??
  # TODO: rename matches to score for consistency

  # input validation
  if (!rlang::is_missing(path_save)) {
    if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)
  }

  if (!is.null(adt_features)) {
    adt_features <- intersect(adt_features, rownames(seurat_obj@assays$ADT))
  }

  if (!all(c("category", "category_details") %in% names(seurat_obj@misc))) {
    seurat_obj@misc$category <- ""
    seurat_obj@misc$category_details <- ""
    cli::cli_inform(c("i" = "The 'category' and 'category_details' fields were not found in the Seurat object's Misc() slot.",
                      "i" = "Filling with empty strings. Please populate these fields for better organization of results."))
  }

  neighbors <- seurat_obj@neighbors[[nn_name]]
  nn_idx <- neighbors@nn.idx
  n_cells <- nrow(nn_idx)
  k <- ncol(nn_idx)

  n_perm <- if (permute) n_permutations else 1
  all_permutation_results <- vector("list", n_perm)

  for (perm_idx in seq_len(n_perm)) {
    results_list <- vector("list", length(meta_cols) +
                             length(adt_features) * length(adt_methods))
    i <- 1

    # calculate matches for categorical/discrete metadata
    for (meta_col in meta_cols) {
      meta_group <- seurat_obj@meta.data[[meta_col]]
      assay <- ifelse(meta_col %in% str_c("annotated_clusters_",
                                          c("binary", "gex_bcr", "simpler")),
                      "GEX", "BCR")

      if (permute) {
        meta_group <- meta_group[sample(n_cells, replace = FALSE)]
      }

      scores_mixing <- vapply(seq_len(n_cells), function(j) {
        current_group <- meta_group[j]
        neighbor_groups <- meta_group[nn_idx[j, ]]
        score <- sum(neighbor_groups == current_group, na.rm = TRUE) / k
        if (is.na(score)) NA_real_ else score
      }, FUN.VALUE = numeric(1))

      valid_cells <- !is.na(scores_mixing)

      if (any(scores_mixing > 1, na.rm = TRUE)) {
        cli::cli_abort(paste0("The non-ADT scores are too high, please check the ",
                    "calculation (the detected k is ", k, ")"))
      }

      results_list[[i]] <-
        data.frame(cell_id = seurat_obj$cell_id[valid_cells], Assay = assay,
                   Meta_Col = meta_col, Method = "Exact",
                   Matches = scores_mixing[valid_cells])
      i <- i + 1
    }

    # calculate matches for ADT features (if provided)
    if (!is.null(adt_features)) {
      assay_name <- substr(nn_name, 1, nchar(nn_name) - 3)
      adt_features <- intersect(adt_features, rownames(seurat_obj@assays$ADT))

      for (feat in adt_features) {
        if ("mean_abs" %in% adt_methods) {
          if (permute) {
            adt_expr <- as.numeric(seurat_obj@assays$ADT@data[feat, ])
            adt_expr <- adt_expr[sample(n_cells, replace = FALSE)]

            scores_adt <- vapply(seq_len(n_cells), function(j) {
              idx <- nn_idx[j, ]
              idx <- idx[idx != j]
              if (length(idx) == 0) return(NA_real_)
              mean(abs(adt_expr[idx] - adt_expr[j]))
            }, FUN.VALUE = numeric(1))
          } else {
            scores_adt <-
              calc_adt_dists(seurat_obj = seurat_obj,
                             base_assay = assay_name, adt_assay = "ADT",
                             features_adt = feat, k = k,
                             distance_metric = "mean_abs",
                             return_mean = FALSE)
          }

          results_list[[i]] <-
            data.frame(cell_id = seurat_obj$cell_id, Assay = "ADT",
                       Meta_Col = feat, Method = "Mean_Absolute",
                       Matches = scores_adt)
          i <- i + 1
        }

        if ("range" %in% adt_methods) {
          if (permute) {
            adt_expr <- as.numeric(seurat_obj@assays$ADT@data[feat, ])
            adt_expr <- adt_expr[sample(n_cells, replace = FALSE)]

            scores_adt <- vapply(seq_len(n_cells), function(j) {
              idx <- nn_idx[j, ]
              idx <- idx[idx != j]
              if (length(idx) == 0) return(NA_real_)
              lower_bound <- adt_expr[j] * (1 - adt_range)
              upper_bound <- adt_expr[j] * (1 + adt_range)
              sum(adt_expr[idx] >= lower_bound &
                    adt_expr[idx] <= upper_bound) / k
            }, FUN.VALUE = numeric(1))
          } else {
            scores_adt <-
              calc_adt_nn_within_range(seurat_obj = seurat_obj,
                                       features_adt = feat,
                                       base_assay = assay_name, k = k,
                                       range = adt_range)
          }

          results_list[[i]] <-
            data.frame(cell_id = seurat_obj$cell_id, Assay = "ADT",
                       Meta_Col = feat, Method = "Range", Matches = scores_adt)
          i <- i + 1
        }
      }
    }

    all_permutation_results[[perm_idx]] <-
      do.call(rbind, results_list[seq_len(i - 1)])
  }

  neighbor_matches <- do.call(rbind, all_permutation_results)

  if (permute) {
    neighbor_matches <- neighbor_matches %>%
      group_by(cell_id, Assay, Meta_Col, Method) %>%
      summarize(Matches = mean(Matches, na.rm = TRUE), .groups = "drop")
  }

  if (!rlang::is_missing(previous_matches)) {
    neighbor_matches <- bind_rows(previous_matches, neighbor_matches) %>%
      arrange(Meta_Col)
  }

  # TODO: make sure it is always a character (could be double)
  category <- seurat_obj@misc$category
  category_details <- seurat_obj@misc$category_details

  neighbor_matches <- neighbor_matches %>%
    mutate(Category = category, Category_Details = category_details,
           .before = Assay) %>%
    mutate(Full_Name = str_c(Category, Category_Details, sep = "_"),
           .before = Category)

  neighbor_matches <- neighbor_matches %>%
    mutate(Full_Name =
             str_replace_all(Full_Name, c("Concatenated_" = "", "_" = " ")),
           Meta_Col = str_replace_all(Meta_Col, c("_" = " ")),
           Category = case_when(str_detect(Category_Details, "Before_PCA") ~
                                  paste0(Category, "_Before_PCA"),
                                str_detect(Category_Details, "After_PCA") ~
                                  paste0(Category, "_After_PCA"),
                                TRUE ~ Category),
           Category = str_replace_all(Category, "_", " "))

  neighbor_matches <-
    neighbor_matches %>%
    mutate(Meta_Col = ifelse(Assay == "ADT", toupper(Meta_Col), Meta_Col))

  if (return_mean) {
    neighbor_matches <-
      neighbor_matches %>%
      group_by(Full_Name, Category, Category_Details, Assay, Meta_Col,
               Method) %>%
      summarize(Matches = mean(Matches, na.rm = TRUE), .groups = "drop")
  }

  # don't need the rownames
  neighbor_matches <- neighbor_matches %>% tibble::remove_rownames()

  # TODO: add the graphs too, or calculate on all graphs
  cli::cli_inform("Scores calculated on {k} neighbors for \\
                  {str_replace_all(seurat_obj@misc$category, '_', ' ')} \\
                  {str_replace_all(seurat_obj@misc$category_details, '_', ' ')}.")

  if (!rlang::is_missing(path_save)) {
    file_name <- tolower(str_c(category, category_details, sep = "_"))
    path_save <- file.path(path_save, paste0(file_name, ".qdata"))
    qd_save(neighbor_matches, path_save)

    cli::cli_inform(c("v" = "Saved neighbor matches to {path_save}"))
  } else {
    return(neighbor_matches)
  }
}
