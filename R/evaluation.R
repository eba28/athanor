#' Compute mean ADT distance to each cell's k nearest neighbors
#'
#' @description
#' This function calculates the mean distance (or similarity) of each cell's ADT profile to its k nearest neighbors in a specified neighbor space (e.g. RNA, BCR, WNN) using a chosen distance metric (e.g. mean absolute difference, Manhattan distance, Euclidean distance, or cosine similarity).
#'
#' @details
#' Pearson correlation coefficient measures linear relationships.
#' Euclidean distance measures the "straight-line" distance.
#' Cosine similarity measures the cosine of the angle between two vectors, with values closer to 1 meaning that they are more similar.
#'
#' @param seurat_obj Seurat object that contains neighbor slots for different k's.
#' @param base_assay Which neighbor space to use (e.g. "RNA", "BCR", "w").
#' @param adt_assay ADT assay name (ADT, ADTnorm).
#' @param layer Data layer to pull from (data, counts, or scale.data).
#' @param feature Optional vector of ADT features; if missing, use all features present.
#' @param k Number of nearest neighbors.
#' @param multiple_k Whether to look for a neighbor slot specific to the provided k (e.g. "RNA.nn_20") or just use the generic one (e.g. "RNA.nn"). The former allows you to have multiple neighbor graphs with different k's, while the latter assumes you only have one neighbor graph per assay.
#' @param distance_metric One of `c(mean_abs, manhattan, euclidean, cosine)`.
#' @param return_mean If TRUE, return the mean across all cells; else return per-cell.
#' @param exclude_self Drop the cell itself from neighbors if present.
#'
#' @returns A single numeric value if `return_mean = TRUE`, or a named numeric vector of per-cell distances if `return_mean = FALSE`.
#' @export
calc_adt_dists <- function(seurat_obj, base_assay, adt_assay = "ADT",
                           layer = "data", feature, k, multiple_k = TRUE,
                           distance_metric = "mean_abs", return_mean = TRUE,
                           exclude_self = TRUE) {
  # use the proper names for the base assays
  if (base_assay == "GEX") base_assay <- "RNA"
  if (base_assay == "WNN") base_assay <- "w"

  # pull the ADT matrix (features by cells) from the specified assay
  # convert to cells by features
  adt_data <- GetAssayData(seurat_obj, assay = adt_assay, layer = layer)

  # make sure that the requested features are available, then get the ADT data
  if (!is_missing(feature)) {
    features_keep <- intersect(feature, rownames(adt_data))

    if (length(features_keep) == 0) {
      stop("The requested ADT feature is not present in assay '",
           adt_assay, "'.")
    }

    # just use the provided feature
    adt_data <- adt_data[features_keep, , drop = FALSE]
    adt_data <- t(as.matrix(adt_data)) # cells x features
  }

  # get the kNNs from the object for the requested base assay and k
  if (multiple_k) nn_name <- paste0(base_assay, ".nn_", k)
  else nn_name <- paste0(base_assay, ".nn")
  if (!nn_name %in% names(seurat_obj@neighbors)) {
    stop("Neighbors slot '", nn_name, "' not found in object. Available: ",
         paste(names(seurat_obj@neighbors), collapse = ", "))
  }

  nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx
  if (nrow(nn_idx) != nrow(adt_data)) {
    stop("Neighbor index rows (", nrow(nn_idx),
         ") do not match with the number of cells (", nrow(adt_data),
         "). Please ensure that the cell order aligns with Cells(seurat_obj).")
  }

  # "Proper per-neighbor distances: compute Euclidean distance per neighbor and average them, not a single distance to all neighbors combined."
  # for each cell, compute Euclidean distance to each neighbor, then average
  k_dist <- vapply(seq_len(nrow(adt_data)), function(cell_i) {
    idx <- nn_idx[cell_i, ]

    # so actually it's the k-1 neighbors (except for WNN)
    if (exclude_self) {
      idx <- idx[idx != cell_i]
      if (length(idx) == 0) return(NA)
    }

    adt_cell <- as.numeric(adt_data[cell_i, feature, drop = FALSE])
    adt_neighbors <- as.numeric(adt_data[idx, feature, drop = FALSE])

    # for Minkowski distances
    # `length(adt_neighbors) should be k or k - 1
    adt_cell_vector <- rep(adt_cell, length(adt_neighbors))

    if (distance_metric == "mean_abs") {
      # I can use absolute value since it's only one feature at a time
      mean(abs(adt_neighbors - adt_cell))
      # diffs <- adt_neighbors - adt_cell
      # mean(sqrt(diffs * diffs))
    } else if (distance_metric == "manhattan") {
      # same as sum(abs_diff)
      dist(rbind(adt_cell_vector, adt_neighbors), method = "minkowski", p = 1)
    } else if (distance_metric == "euclidean") {
      # good for continuous data, overemphasizes big diffs
      # same as sqrt(sum(abs_diff * abs_diff))
      dist(rbind(adt_cell_vector, adt_neighbors), method = "minkowski", p = 2)
    } else if (distance_metric == "cosine") {
      lsa::cosine(x = adt_cell_vector, y = adt_neighbors)
    } else {
      stop("Please choose a valid distance metric.")
    }
  }, FUN.VALUE = numeric(1))

  # return the final distances
  if (return_mean) {
    return(mean(k_dist, na.rm = TRUE))
  } else {
    names(k_dist) <- rownames(adt_data)
    return(k_dist)
  }
}


#' Compute mean ADT distance to each cell's k nearest neighbors (faster version)
#'
#' @description
#' This function calculates the mean distance of each cell's ADT profile to its k nearest neighbors in a specified neighbor space (e.g. RNA, BCR, WNN) using Euclidean distance. It returns a named numeric vector of per-cell mean distances.
#'
#' @details
#' Only returns the mean.
#'
#' @param adt_data ADT data matrix (features by cells) to use for distance calculations.
#' @param features Optional vector of ADT features; if missing, use all features present.
#' @param neighbors The kNN slot.
#' @param exclude_self Drop the cell itself from neighbors if present.
#' @export
calc_adt_dists_fast <- function(adt_data, features, neighbors,
                                exclude_self = TRUE) {
  if (!is_missing(features)) {
    keep <- intersect(features, rownames(adt_data))
    adt_data <- adt_data[keep, , drop = FALSE]
  }
  adt_data <- t(as.matrix(adt_data))

  nn_idx <- neighbors@nn.idx

  mean_k_dist <- vapply(seq_len(nrow(adt_data)), function(i) {
    idx <- nn_idx[i, ]

    if (exclude_self) {
      idx <- idx[idx != i]
      if (length(idx) == 0) return(NA_real_)
    }

    Xi <- adt_data[i, , drop = FALSE]
    nbrs <- adt_data[idx, , drop = FALSE]
    diffs <- sweep(nbrs, 2, Xi, FUN = "-")
    sqrt_rowSums <- sqrt(rowSums(diffs * diffs))
    mean(sqrt_rowSums)
  }, numeric(1))

  names(mean_k_dist) <- colnames(adt_mat)
  return(mean_k_dist)
}


#' Calculate the proportion of neighbors within an ADT marker's expression range
#'
#' @description
#' For each cell in a Seurat object, this function calculates how many of its
#' k nearest neighbors have ADT expression within a specified threshold (default
#' 20%) of the cell's own ADT expression for a given feature.
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
#' @param feature Name of the ADT feature to evaluate (e.g. "CD27.1", "CD38").
#' @param base_assay The assay used to compute neighbors. One of
#'   "RNA", "GEX", "BCR", or "WNN".
#' @param k Number of nearest neighbors to evaluate. Must match the k
#'   used when computing the neighbor graph.
#' @param use_k Whether to look for a neighbor slot specific to the provided k (e.g. "RNA.nn_20") or just use the generic one (e.g. "RNA.nn"). The former allows you to have multiple neighbor graphs with different k's, while the latter assumes you only have one neighbor graph per assay.
#' @param range The relative threshold for considering neighbors
#'   similar. A value of 0.20 means neighbors within ±20% of the cell's
#'   expression are counted.
#' @param return_counts If TRUE, returns the count of neighbors within
#'   range. If FALSE, returns the proportion (count/k).
#'
#' @return A named numeric vector with one value per cell in the Seurat object.
#'   If \code{return_counts = TRUE}, returns the count of neighbors within
#'   range. If \code{return_counts = FALSE}, returns the proportion of
#'   neighbors within range (ranging from 0 to 1). Vector names are cell ids.
#' @export
calc_adt_nn_within_range <- function(seurat_obj, adt_assay = "ADT", feature,
                                     base_assay, k = 20, use_k = TRUE,
                                     range = 0.20, return_counts = FALSE) {
  # TODO: standardize "feature" vs. "adt_features"

  # get the neighbors for the specified assay
  if (rlang::is_missing(base_assay)) base_assay <- DefaultAssay(seurat_obj)
  assay_name <- recode_values(base_assay, "GEX" ~ "RNA", "WNN" ~ "w",
                              default = base_assay)
  assay_name <- paste0(assay_name, ".nn")

  # just don't include k
  if (use_k) assay_name <- paste0(assay_name, "_", k)

  # pull the nn information from the object
  neighbors <- seurat_obj@neighbors[[assay_name]]
  nn_idx <- Indices(neighbors)

  if (is.null(neighbors)) {
    stop(paste("No neighbors found for assay", base_assay, "with k =", k))
  }

  # get ADT expression for the given feature
  adt_expr <- seurat_obj@assays[[adt_assay]]@data[feature, ]

  # calculate the counts or proportion for each cell in the object
  results <- sapply(1:nrow(nn_idx), function(i) {
    cell_expr <- adt_expr[i]
    neighbor_expr <- adt_expr[nn_idx[i, ]]

    # calculate the range bounds
    lower_bound <- cell_expr * (1 - range)
    upper_bound <- cell_expr * (1 + range)

    # count neighbors within the range
    within_range <-
      sum(neighbor_expr >= lower_bound & neighbor_expr <= upper_bound)

    if (return_counts) {
      return(within_range)
    } else {
      return(within_range / k)
    }
  })

  names(results) <- Cells(seurat_obj)
  return(results)
}


#' Calculate the proportion of neighbors within an ADT marker's quantile by expression
#'
#' @description
#' For each cell in a Seurat object, this function calculates how many of its
#' k nearest neighbors have ADT expression within the same quantile as the cell's own ADT expression for a given feature.
#'
#' @param seurat_obj A Seurat object containing ADT data and computed neighbor
#'   graphs.
#' @param adt_assay Name of the assay containing ADT data.
#' @param feature Name of the ADT feature to evaluate (e.g. "CD27.1", "CD38").
#' @param base_assay The assay used to compute neighbors. One of
#'   "RNA", "GEX", "BCR", or "WNN".
#' @param k Number of nearest neighbors to evaluate. Must match the k
#'   used when computing the neighbor graph.
#' @param use_k Whether to look for a neighbor slot specific to the provided k (e.g. "RNA.nn_20") or just use the generic one (e.g. "RNA.nn"). The former allows you to have multiple neighbor graphs with different k's, while the latter assumes you only have one neighbor graph per assay.
#' @param n_quantile The number of quantiles to divide the ADT expression into. Neighbors are considered "within range" if they fall into the same quantile as the cell.
#' @param method One of `c("quantile", "percentile_diff")`.
#'   `"quantile"` compares discrete quantile bins (current behavior).
#'   `"percentile_diff"` returns mean absolute percentile rank difference per cell.
#' @param return_counts If TRUE, returns the count of neighbors within
#'   the same quantile (only applies to `method = "quantile"`). If FALSE, returns
#'   the proportion (count/k). Ignored when `method = "percentile_diff"`.
#'
#' @return A named numeric vector with one value per cell in the Seurat object.
#'   For `method = "quantile"`: If \code{return_counts = TRUE}, returns the count
#'   of neighbors within the same quantile bin. If \code{return_counts = FALSE},
#'   returns the proportion (ranging from 0 to 1). For `method = "percentile_diff"`:
#'   returns mean absolute percentile rank difference per cell (ranging from 0 to 1).
#'   Vector names are cell ids.
#' @export
calc_adt_quantile <- function(seurat_obj, adt_assay = "ADT", feature,
                              base_assay, k = 20, use_k = TRUE, n_quantile = 10,
                              method = c("quantile", "percentile_diff"),
                              return_counts = FALSE) {
  method <- match.arg(method)

  # get the neighbors for the specified assay
  if (rlang::is_missing(base_assay)) base_assay <- DefaultAssay(seurat_obj)
  assay_name <- recode_values(base_assay, "GEX" ~ "RNA", "WNN" ~ "w",
                              default = base_assay)
  assay_name <- paste0(assay_name, ".nn")

  # just don't include k
  if (use_k) assay_name <- paste0(assay_name, "_", k)

  # pull the nn information from the object
  neighbors <- seurat_obj@neighbors[[assay_name]]

  if (is.null(neighbors)) {
    stop(paste("No neighbors found for assay", base_assay, "with k =", k))
  }
  nn_idx <- Indices(neighbors)

  # get ADT expression for the given feature
  adt_expr <- seurat_obj@assays[[adt_assay]]@data[feature, ]
  adt_expr <- as.numeric(adt_expr)

  if (method == "quantile") {
    adt_quantiles <- cut_interval(adt_expr, n = n_quantile)
    adt_quantiles <- as.numeric(adt_quantiles)

    # calculate the counts or proportion for each cell in the object
    results <- sapply(1:nrow(nn_idx), function(i) {
      cell_quantile <- adt_quantiles[i]
      neighbor_quantile <- adt_quantiles[nn_idx[i, ]]

      # count neighbors within the same quantile bin
      within_range <- sum(neighbor_quantile == cell_quantile, na.rm = TRUE)

      if (return_counts) {
        return(within_range)
      } else {
        return(within_range / k)
      }
    })
  } else {
    # empirical percentile rank in [0, 1]
    adt_percentiles <- rep(NA, length(adt_expr))
    valid <- !is.na(adt_expr)
    n_valid <- sum(valid)

    if (n_valid <= 1) {
      adt_percentiles[valid] <- 0
    } else {
      ranked <- rank(adt_expr[valid], ties.method = "average")
      adt_percentiles[valid] <- (ranked - 1) / (n_valid - 1)
    }

    results <- sapply(1:nrow(nn_idx), function(i) {
      cell_percentile <- adt_percentiles[i]
      neighbor_percentile <- adt_percentiles[nn_idx[i, ]]

      if (is.na(cell_percentile)) return(NA)

      pct_diff <- abs(neighbor_percentile - cell_percentile)
      mean(pct_diff, na.rm = TRUE)
    })

    if (return_counts) {
      warning("return_counts = TRUE is not applicable for method = 'percentile_diff'. Returning mean absolute percentile rank difference instead.")
    }

    # convert to similarity by subtracting from 1
    # results <- 1 - results
    return(results)
  }
}


# cluster_distances <- function(embeddings, labels) {
#   clusters <- split(as.data.frame(embeddings), labels)
#
#   # Intra-cluster distance (average pairwise distance within each cluster)
#   intra_distances <- sapply(clusters, function(cluster_data) {
#     mean(dist(cluster_data))
#   })
#
#   # Inter-cluster distance (average distance between cluster centroids)
#   centroids <- sapply(clusters, colMeans)
#   inter_distances <- dist(t(centroids))  # Distance between centroids
#
#   return(list(intra = intra_distances, inter = inter_distances))
# }
#
# # Apply to Seurat data
# umap_embeddings <- Embeddings(seurat_obj, reduction = reduction_name)
# predicted_labels <- seurat_obj$annotated_clusters_simpler
# distances <- cluster_distances(umap_embeddings, predicted_labels)


#' Calculate the correlation between each cell's expression and the mean of its neighbors' expression
#'
#' @details
#' The Seurat object must have `FindNeighbors()` already run at least one time and an assay named "ADT".
#'
#' @param seurat_obj The Seurat object.
#' @param features_adt Name of the ADT features to evaluate on (e.g. "CD27.1").
#' @param cor_method Correlation method to use (e.g. "pearson", "spearman").
#'
#' @returns A data frame with columns: Graph, Feature, Score (correlation value).
calc_correlation <- function(seurat_obj, features_adt, cor_method = "spearman") {
  # TODO: return values or a df
  # TODO: run for each feature or all features
  # TODO: run in parallel

  # make sure that the ADT feature is in the object
  if (rlang::is_missing(features_adt)) {
    features_adt <- rownames(seurat_obj@assays$ADT)
  }

  if (any(!(features_adt %in% rownames(seurat_obj@assays$ADT)))) {
    stop("The requested ADT feature is not present in assay 'ADT'. Available features: ",
         paste(sort(rownames(seurat_obj@assays$ADT)), collapse = ", "))
  }
  # make sure that the object has at least one neighbor graph
  if (length(seurat_obj@neighbors) == 0) {
    stop("No neighbor graphs found in object. Please run FindNeighbors() first.")
  }

  # get the (normalized) ADT expression
  adt_data <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")
  metrics_df <- c()

  for (nn_name in names(seurat_obj@neighbors)) {
    # create a matrix for the neighbor-averaged ADT expression
    nn_idx <- seurat_obj@neighbors[[nn_name]]@nn.idx
    neighbor_adt_mean <- sapply(1:nrow(nn_idx), function(i) {
      rowMeans(adt_data[, nn_idx[i, ]])
    })
    colnames(neighbor_adt_mean) <- colnames(seurat_obj)

    # calculate correlations for each feature across all of the cells
    for (feature_adt in features_adt) {
      cell_expr <- adt_data[feature_adt, ]
      neighbor_expr <- neighbor_adt_mean[feature_adt, ]

      correlation <- cor(cell_expr, neighbor_expr, method = cor_method)

      metrics_df <- bind_rows(metrics_df,
                              data.frame(Graph = nn_name, Feature = feature_adt,
                                         Score = correlation))
    }
  }

  return(metrics_df)
}


#' Calculate Moran's i for a Seurat object
#'
#' @description
#' This function calculates the global Moran's i index.
#'
#' @details
#' We are using `MERINGUE`'s implementation instead of `ape`'s because it runs faster.
#' However, `MERINGUE` is not on CRAN, which means this package could not be published on CRAN.
#' Row standardization makes sure that the resulting score will always be between -1 and 1.
#'
#' @param seurat_obj The Seurat object. Must have `FindNeighbors()` already run and an assay named "ADT".
#' @param feature Name of the ADT feature to evaluate (e.g. "CD27.1", "CD38").
#' @param graph_name Name of the neighbor graph slot to use for the weights matrix (e.g. "RNA.nn", "BCR.nn", "w.nn").
#' @param row_standardize Whether or not to row-standardize the weights matrix (i.e. make each row sum to 1).
#'
#' @returns A single numeric value representing the observed Moran's i index for the specified feature and neighbor graph.
#' @export
calc_moran <- function(seurat_obj, feature, graph_name, row_standardize = TRUE) {
  # TODO: provide multiple features
  # TODO: provide multiple graphs

  # make sure that the ADT feature is in the object
  if (!feature %in% rownames(seurat_obj@assays$ADT)) {
    stop("The requested ADT feature is not present in assay 'ADT'. Available features: ",
         paste(sort(rownames(seurat_obj@assays$ADT)), collapse = ", "))
  }
  # make sure that the graph is in the object
  if (!graph_name %in% names(seurat_obj@graphs)) {
    stop("Graph '", graph_name, "' not found in object. Available graphs: ",
         paste(names(seurat_obj@graphs), collapse = ", "))
  }

  # get the (normalized) expression
  x <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")
  x <- x[feature, ]

  # get the weights matrix
  w <- seurat_obj@graphs[[graph_name]] # as.matrix

  # row-standardize the weights (so each row sums to 1)
  # this shouldn't be necessary because we are comparing graphs with the same k
  if (row_standardize) w <- w / rowSums(w)

  # calculate Moran's i
  moran <- ape::Moran.I(x, w)
  # moran <- MERINGUE::moranTest(x, w)

  # if you want to calculate the local value:
  # list_w <- mat2listw(as.matrix(w)) # spdep
  # local_m <- localmoran(marker_vector, list_w) # produces a value for each cell
  # add to Seurat metadata for plotting
  # seurat_obj@meta.data[[paste0("local_i_", feature)]] <-
  #   local_m[, "Ii"] # 'Ii' is the local Moran statistic

  # just return the observed Moran's i value, not the expected value, sd, or p-value
  return(moran[["observed"]])
}


#' Calculate neighbor matching scores across metadata columns
#'
#' @description
#' This function calculates the proportion of neighbors that match each cell's metadata category for specified metadata columns and ADT features.
#' It can also perform permutations to generate a random baseline for comparison. The results can be saved to a specified path.
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
#' @param permute Shuffle labels for each meta column and ADT expression per feature before computing matches to get a random baseline.
#' @param n_permutations The number of times to permute labels.
#' @param previous_matches Data frame of previous matches to combine with the new results (optional).
#' @param path_save Where to save the results.
#'
#' @return Data frame with these columns: Full_Name, Category, Category_Details, Assay, Meta_Col, Method, Matches
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
                                  previous_matches, path_save) {
  # make sure that the directory exists first if saving
  if (!rlang::is_missing(path_save)) {
    if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)
  }

  # make sure that all of the ADT features are in the object
  if (!is.null(adt_features)) {
    adt_features <- intersect(adt_features, rownames(seurat_obj@assays$ADT))
  }

  # make sure that the category is in the object
  if (!all(c("category", "category_details") %in% names(seurat_obj@misc))) {
    stop("Please fill the 'category' and 'category_details' miscellaneous slots.")
  }

  # get neighbor indices
  # nn_name <- paste(DefaultAssay(seurat_obj))
  neighbors <- seurat_obj@neighbors[[nn_name]]
  nn_idx <- Indices(neighbors)
  n_cells <- nrow(nn_idx)
  k <- ncol(neighbors)

  # determine number of permutations to run and store results
  n_perm <- if (permute) n_permutations else 1
  all_permutation_results <- vector("list", n_perm)

  # run multiple permutations if permute = TRUE
  for (perm_idx in seq_len(n_perm)) {
    # pre-allocate a list to store results (more memory efficient than bind_rows)
    results_list <- vector("list", length(meta_cols) +
                             length(adt_features) * length(adt_methods))
    i <- 1

    # calculate matches for categorical/discrete metadata
    # you could also do a range for values like cdr3_aa_length
    for (meta_col in meta_cols) {
      meta_group <- seurat_obj[[]][[meta_col]]
      assay <- ifelse(meta_col %in% str_c("annotated_clusters_",
                                          c("binary", "gex_bcr", "simpler")),
                      "GEX", "BCR")

      # I already set the random seed during setup
      if (permute) {
        meta_group <- meta_group[sample(n_cells, replace = FALSE)]
      }

      # faster than a for loop and takes up less memory
      scores_mixing <- vapply(1:n_cells, function(i) {
        current_group <- meta_group[i]
        neighbor_groups <- meta_group[nn_idx[i, ]]

        # proportion of neighbors that match the current cell
        score <- sum(neighbor_groups == current_group, na.rm = TRUE) / k
        if (is.na(score)) NA else score}, numeric(1))

      # remove NAs
      valid_cells <- !is.na(scores_mixing)

      # sanity check
      if (any(scores_mixing > 1)) {
        stop(paste0("The non-ADT scores are too high, please check the calculation",
                    " (the detected k is ", k, ")"))
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
        # TODO: include Euclidean and Manhattan distances too

        if ("mean_abs" %in% adt_methods) {
          # TODO: check this
          if (permute) {
            adt_expr <- as.numeric(seurat_obj@assays$ADT@data[feat, ])
            adt_expr <- adt_expr[sample(n_cells, replace = FALSE)]

            scores_adt <-
              vapply(seq_len(n_cells), function(i) {
                idx <- nn_idx[i, ]
                idx <- idx[idx != i]
                if (length(idx) == 0) return(NA)

                cell_expr <- adt_expr[i]
                neighbor_expr <- adt_expr[idx]
                mean(abs(neighbor_expr - cell_expr))}, FUN.VALUE = numeric(1))
          } else {
            scores_adt <-
              calc_adt_dists(seurat_obj = seurat_obj, base_assay = assay_name,
                             adt_assay = "ADT", feature = feat, k = k,
                             multiple_k = FALSE, distance_metric = "mean_abs",
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

            scores_adt <-
              vapply(seq_len(n_cells), function(i) {
                idx <- nn_idx[i, ]
                idx <- idx[idx != i]
                if (length(idx) == 0) return(NA)

                cell_expr <- adt_expr[i]
                neighbor_expr <- adt_expr[idx]
                lower_bound <- cell_expr * (1 - adt_range)
                upper_bound <- cell_expr * (1 + adt_range)

                sum(neighbor_expr >= lower_bound & neighbor_expr <= upper_bound) / k},
                FUN.VALUE = numeric(1))
          } else {
            scores_adt <-
              calc_adt_nn_within_range(seurat_obj = seurat_obj, feature = feat,
                                       base_assay = assay_name, k = k,
                                       use_k = FALSE, range = adt_range)
          }

          results_list[[i]] <-
            data.frame(cell_id = seurat_obj$cell_id, Assay = "ADT",
                       Meta_Col = feat, Method = "Range", Matches = scores_adt)
          i <- i + 1
        }
      }
    }

    # combine results from this permutation
    all_permutation_results[[perm_idx]] <- do.call(rbind, results_list[1:(i - 1)])
  }

  # combine all permutations (if applicable)
  neighbor_matches <- do.call(rbind, all_permutation_results)

  if (permute) {
    # average across permutations
    neighbor_matches <- neighbor_matches %>%
      group_by(cell_id, Assay, Meta_Col, Method) %>%
      summarize(Matches = mean(Matches, na.rm = TRUE), .groups = "drop")
  }

  if (!rlang::is_missing(previous_matches)) {
    neighbor_matches <- bind_rows(previous_matches, neighbor_matches) %>%
      arrange(Meta_Col)
  }

  # add category columns
  category <- seurat_obj@misc$category
  category_details <- seurat_obj@misc$category_details
  neighbor_matches <-
    neighbor_matches %>%
    mutate(Category = category, Category_Details = category_details,
           .before = Assay) %>%
    mutate(Full_Name = str_c(Category, Category_Details, sep = "_"),
           .before = Category)

  # simplify and standardize the category names
  neighbor_matches <-
    neighbor_matches %>%
    mutate(Full_Name =
             str_replace_all(Full_Name, c("Concatenated_" = "", "_" = " ")),
           Meta_Col = str_replace_all(Meta_Col, c("_" = " ")),
           Category = case_when(str_detect(Category_Details, "Before_PCA") ~
                                  paste0(Category, "_Before_PCA"),
                                str_detect(Category_Details, "After_PCA") ~
                                  paste0(Category, "_After_PCA"),
                                TRUE ~ Category),
           Category = str_replace_all(Category, "_", " ")) # for plotting

  # standardize the ADT feature names
  neighbor_matches <-
    neighbor_matches %>%
    mutate(Meta_Col = ifelse(Assay == "ADT", toupper(Meta_Col), Meta_Col))

  # neighbor_matches <-
  #   neighbor_matches %>%
  #   mutate(Assay =
  #          case_when(Meta_Col %in% c("annotated clusters binary",
  #                                    "annotated clusters simpler") ~ "GEX",
  #                    Meta_Col == "annotated clusters gex bcr" ~ "GEX_BCR",
  #                    Assay == "BCR_GEX" ~ "BCR",
  #                    Assay == "ADT" ~ paste0(Assay, "_", Method),
  #                    .default = Assay))

  # useful printout
  cat(paste("Scores calculated on", k, "neighbors for",
            str_replace_all(obj@misc$category, "_", " "),
            str_replace_all(obj@misc$category_details, "_", " "), "\n"))

  # save or return the results
  # TODO: address balm_paired vs. BALM-paired (can mess up csv reading)
  if (!rlang::is_missing(path_save)) {
    file_name <- tolower(str_c(category, category_details, sep = "_"))
    qd_save(neighbor_matches,
            file.path(path_save, paste0(file_name, ".qdata")))
  } else {
    return(neighbor_matches)
  }
}
