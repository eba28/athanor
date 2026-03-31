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

#' Calculate cluster distances for a Seurat object
#'
#' @description
#' This function calculates various cluster distance metrics using the `fpc::cluster.stats` function, which provides a comprehensive set of clustering statistics based on a distance matrix and cluster assignments.
#' The function computes both single-value metrics (e.g. Calinski-Harabasz index) and cluster-wise metrics (e.g. average within-cluster distance) depending on the specified criteria.
#' The results are returned in a tidy data frame format for easy plotting and comparison across different embeddings and reductions.
#'
#' @details
#' The `as.factor()` is needed in case you give a categorical cluster col.
#' Assumes you're using the embeddings approach.
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param reduction_name The name of the reduction to use for distances.
#' @param criteria One of: Between_Mean, Between_Min, Calinski_Harabasz, Within_Between, Within_Max, Within_Mean, Within_Median.
#' @param labels_true The name of the column in the metadata that contains the cluster labels to evaluate.
#' @param labels_name A more descriptive name for the labels to use in plotting (optional).
#'
#' @returns A data.frame with a row per metric and cluster containing the score.
#' @export
calc_distances <- function(seurat_obj, reduction_name, criteria = "Within_Max",
                           labels_true = "annotated_clusters", labels_name) {
  # set up the inputs and output
  embeddings <- Embeddings(seurat_obj, reduction = reduction_name)
  dist_matrix <- stats::dist(embeddings) # method = "euclidean"
  clusters <- seurat_obj[[]][[labels_true]]
  clusters_fact <- as.factor(clusters)
  clusters_int <- as.integer(clusters_fact)

  # some metrics are per label group and some are for everything
  distances_single <- c()
  distances_multi <- c()

  # check for issues before calling cluster.stats
  if (any(is.na(dist_matrix))) {
    cat("Warning: Distance matrix contains NAs.\n")
    # Option: Remove NAs or use complete cases only
  }
  if (any(is.na(clusters_int))) {
    cat("Warning: Clustering vector contains NAs.\n") # should never happen
    valid_indices <- !is.na(clusters_int)
    clusters_int <- clusters_int[valid_indices]
    # may also need to subset the distance matrix accordingly
  }
  if (xfun::attr2(dist_matrix, "Size") != length(clusters_int)) {
    cat("Length mismatch between the distance matrix and the clustering.\n")
  }

  # call cluster.stats with error handling
  tryCatch({
    cluster_stats <- fpc::cluster.stats(d = dist_matrix, clustering = clusters_int)
  }, error = function(e) {
    cat("Error in cluster.stats:", e$message, "\n")
    return(NULL)
  })

  # calculate cluster_stats
  # cluster_stats <- fpc::cluster.stats(d = dist_matrix, clustering = clusters_int)

  if ("Between_Mean" %in% criteria) {
    # "vector of clusterwise average distances of a point in the cluster to the points of other clusters"
    distances_multi <-
      bind_rows(distances_multi,
                data.frame(Metric = "Between_Mean",
                           Score = cluster_stats$average.toother))
  }

  if ("Between_Min" %in% criteria) {
    # "vector of clusterwise minimum distances of a point in the cluster to a point of another cluster"
    distances_multi <-
      bind_rows(distances_multi,
                data.frame(Metric = "Between_Min",
                           Score = cluster_stats$separation))
  }

  if ("Calinski_Harabasz" %in% criteria) {
    # "Calinski and Harabasz index (Calinski and Harabasz 1974, optimal in Milligan and Cooper 1985; generalised for dissimilarites in Hennig and Liao 2013)."
    distances_single <-
      bind_rows(distances_single,
                data.frame(Metric = "Calinski_Harabasz",
                           Score = cluster_stats$ch))
  }

  if ("Within_Between" %in% criteria) {
    # "average.within/average.between"
    distances_single <-
      bind_rows(distances_single,
                data.frame(Metric = "Within_Between",
                           Score = cluster_stats$wb.ratio))
  }

  if ("Within_Max" %in% criteria) {
    # "vector of cluster diameters (maximum within cluster distances)"
    distances_multi <-
      bind_rows(distances_multi,
                data.frame(Metric = "Within_Max",
                           Score = cluster_stats$diameter))
  }

  if ("Within_Mean" %in% criteria) {
    # "vector of clusterwise within cluster average distances"
    distances_multi <-
      bind_rows(distances_multi,
                data.frame(Metric = "Within_Mean",
                           Score = cluster_stats$average.distance))
  }

  if ("Within_Median" %in% criteria) {
    # "vector of clusterwise within cluster distance medians"
    distances_multi <-
      bind_rows(distances_multi,
                data.frame(Metric = "Within_Median",
                           Score = cluster_stats$median.distance))
  }

  # add a few more columns (the order might not be right)
  # criteria could include single and multi
  if (!is.null(distances_multi)) {
    distances_multi <-
      distances_multi %>%
      mutate(Clusters = rep(levels(clusters_fact),
                            nrow(distances_multi) / nlevels(clusters_fact)))
  }
  distances <-
    bind_rows(distances_single, distances_multi) %>%
    mutate(Embedding = seurat_obj@misc$embedding_type,
           Reduction = reduction_name, Labeling = labels_true,
           .before = Metric)
  if (!rlang::is_missing(labels_name)) distances$Labeling <- labels_name

  return(distances)
}


#' Calculate internal clustering metrics for a Seurat object
#'
#' @description
#' This function calculates various internal clustering metrics using the `cluster` package, which provides a comprehensive set of clustering evaluation statistics based on a distance matrix and cluster assignments.
#'
#' @details
#' The Satija lab used `cluster` for their analyses.
#' The `as.factor()` is needed in case you give a categorical cluster col.
#' Use the `bluster` package for speed if desired e.g.
#' `bluster::approxSilhouette(x = embeddings, clusters = clusters)`
#' Assumes you're using the embeddings approach.
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param reduction_name The name of the reduction to use for distances.
#' @param criteria One of: DB, Dunn, Intra_Complete, Silhouette
#' @param labels_true The name of the column in the metadata that contains the cluster labels to evaluate.
#' @param labels_name A more descriptive name for the labels to use in plotting (optional).
#' @param return_full If TRUE, return the full silhouette object instead of just the mean silhouette width.
#'
#' @returns A data.frame with a row per metric containing the combined score.
#' @export
calc_int_metrics <- function(seurat_obj, reduction_name,
                             criteria = "Silhouette",
                             labels_true = "annotated_clusters", labels_name,
                             return_full = FALSE) {
  # set up the inputs and output
  embeddings <- Embeddings(seurat_obj, reduction = reduction_name)
  dist_matrix <- dist(embeddings) # method = "euclidean"
  clusters <- seurat_obj[[]][[labels_true]]
  clusters_int <- as.integer(as.factor(clusters))

  metrics <- c()

  if ("Davies-Bouldin" %in% criteria) {
    score_db <- clusterSim::index.DB(embeddings, clusters_int)$DB

    metrics <- bind_rows(metrics,
                         data.frame(Metric = "Davies-Bouldin", Score = score_db))
  }

  # could also be done with fpc::cluster.stats
  if ("Dunn" %in% criteria) {
    score_dunn <- clValid::dunn(distance = dist_matrix, clusters = clusters_int)

    metrics <- bind_rows(metrics,
                         data.frame(Metric = "Dunn", Score = score_dunn))
  }

  if ("Intra_Complete" %in% criteria) {
    # clv is no longer on CRAN
    # cls.scatt.diss.mx
    # score_intra_complete <- clv::cls.scatt.data(data = embeddings,
    #                                             clust = clusters_int,
    #                                             dist = "euclidean")
    score_intra_complete <- mean(score_intra_complete$intracls.complete)

    metrics <- bind_rows(metrics,
                         data.frame(Metric = "Intra_Complete",
                                    Score = score_intra_complete))
  }

  if ("Silhouette" %in% criteria) {
    score_sil <- cluster::silhouette(x = clusters_int, dist_matrix)

    if (!return_full) score_sil <- mean(score_sil[, "sil_width"])

    metrics <- bind_rows(metrics,
                         data.frame(Metric = "Silhouette", Score = score_sil))
  }

  # add a few more columns
  metrics <- metrics %>%
    mutate(Embedding = seurat_obj@misc$embedding_type,
           Reduction = reduction_name, Labeling = labels_true,
           .before = Metric)
  if (!rlang::is_missing(labels_name)) metrics$Labeling <- labels_name

  return(metrics)
}


#' Calculate external clustering metrics for a Seurat object
#'
#' @description
#' This function calculates various external clustering metrics using the `mclust` and `clevr` packages, which provide a comprehensive set of clustering evaluation statistics based on true cluster labels and predicted cluster assignments.
#'
#' @details
#' The Satija lab used `cluster` for their analyses.
#' The `as.factor()` is needed in case you give a categorical cluster col.
#' `sklearn.metrics.cluster` has `completeness_score` and `homogeneity_score`
#' Assumes you're using the embeddings approach.
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param reduction_name The name of the reduction to use for distances.
#' @param criteria One of: ARI, Completeness, Homogeneity
#' @param labels_true The name of the column in the metadata that contains the true cluster labels to evaluate.
#' @param labels_pred The name of the column in the metadata that contains the predicted cluster labels to evaluate.
#'
#' @returns A data.frame with a row per metric containing the combined score.
#' @export
calc_ext_metrics <- function(seurat_obj, reduction_name,
                             criteria = c("Completeness", "Homogeneity"),
                             labels_true = "annotated_clusters",
                             labels_pred = "seurat_clusters") {
  # set up the inputs and output
  clusters_true <- seurat_obj[[]][[labels_true]]
  clusters_pred <- seurat_obj[[]][[labels_pred]]

  metrics <- c()

  # Adjusted Rand Index
  if ("ARI" %in% criteria) {
    metrics <-
      bind_rows(metrics,
                data.frame(Metric = "ARI",
                           Score = mclust::adjustedRandIndex(clusters_pred,
                                                             clusters_true)))
  }

  if ("Completeness" %in% criteria) {
    metrics <-
      bind_rows(metrics,
                data.frame(Metric = "Completeness",
                           Score = clevr::completeness(true = clusters_true,
                                                       pred = clusters_pred)))
  }

  if ("Homogeneity" %in% criteria) {
    metrics <-
      bind_rows(metrics,
                data.frame(Metric = "Homogeneity",
                           Score = clevr::homogeneity(true = clusters_true,
                                                      pred = clusters_pred)))
  }

  # add a few more columns
  metrics <- metrics %>%
    mutate(Embedding = seurat_obj@misc$embedding_type,
           Reduction = reduction_name, Labeling = labels_true,
           .before = Metric)

  return(metrics)
}


#' Plot the internal or external clustering metrics
#'
#' @description
#' This function creates a ggplot visualization of the internal or external clustering metrics across different embeddings and reductions.
#' It uses a color scale to represent the metric scores and outlines the best scores for each label and metric combination, as well as the best scores across reductions.
#' The plot is faceted by metric and reduction for easy comparison.
#'
#' @details
#' Make sure to check that the best score is consistent across your plotting metrics.
#'
#' @param metrics Data frame of metrics to plot, with columns: Embedding, Reduction, Labeling, Metric, Score.
#' @param plot_title Title to use for the plot.
#' @param best_score One of "higher" or "lower" to indicate whether higher or lower scores are better for the metrics being plotted. This is used to determine which scores to outline in the plot.
#' @param type One of "Internal" or "External" to indicate the type of metrics being plotted, used for the plot title.
#' @param y_axis The name of the column in the metrics data frame to use for the y-axis (e.g. "Labeling").
#' @param round_to Number of decimal places to round the score labels to in the plot.
#' @param details Additional details to include in the plot title (optional).
#'
#' @returns A ggplot showing the metrics across embeddings and reductions, with the best scores outlined.
#' @export
plot_metrics <- function(metrics, plot_title = "", best_score = "higher",
                         type = "Internal", y_axis = "Labeling", round_to = 2,
                         details = "") {
  # calculate the top scores
  # could have ties
  top_scores <-
    metrics %>%
    mutate(id = row_number()) %>%
    group_by(Reduction, !!sym(y_axis), Metric) %>%
    group_modify(~ {
      if (best_score == "higher") {
        slice_max(.x, order_by = Score, n = 1, with_ties = FALSE)
      } else {
        slice_min(.x, order_by = Score, n = 1, with_ties = FALSE)
      }}) %>%
    ungroup(Reduction) %>%
    mutate(best_id = ifelse(best_score == "higher",
                            id[which.max(Score)], id[which.min(Score)]))

  # pick the color palette
  palette <- rev(pals::brewer.rdbu(n = 7)) # the color scale function's default
  if (best_score != "higher") palette <- rev(palette)

  ggplot(metrics, aes(x = Embedding, y = !!sym(y_axis), fill = Score)) %>%
    plot_color_scale(palette = palette, val_col = "Score", fill_by = "fill") +
    geom_tile(linewidth = 0.4, color = "white") +
    geom_text(aes(label = round(Score, digits = round_to)), size = 3) +
    # outline the highest score by label, reduction and metric
    geom_tile(data = metrics[top_scores$id, ], fill = NA,
              color = "black", linewidth = 0.4) +
    # outline the highest score across reductions
    geom_tile(data = metrics[top_scores$best_id, ], fill = NA,
              color = "black", linewidth = 0.8) +
    labs(title = paste(plot_title, type, "Clustering", details),
         subtitle = "Black outlines = best scores per label and metric") +
    scale_y_discrete(limits = rev) +
    facet_grid(rows = vars(Metric), cols = vars(Reduction),
               scales = "free_x", space = "free_x") +
    theme_bw + labels_standard +
    theme(panel.grid.major.y = element_blank())
}


#' Calculate homogeneity scores for binary ADT features across embeddings
#'
#' @description
#' This function calculates homogeneity scores for binary ADT features across different embeddings and reductions in a Seurat object.
#'
#' @param seurat_objs List of Seurat objects for each embedding type.
#' @param meta_res Named list of cluster columns for each reduction.
#' @param metric_type One of: Internal, External
#' @param metrics List of metric names.
#' @param adt_features ADT feature name (e.g. "CD27.1"). You can provide multiple.
#' @param adt_cutoff Numeric cutoff for binary classification.
#'
#' @return Data frame of scores for each embedding/reduction.
#' @export
calc_adt_scores <- function(seurat_objs, meta_res, metric_type, metrics,
                            adt_features = "CD27.1", adt_cutoff = 1) {
  scores <- c()

  for (adt_feat in adt_features) {
    # embedding-dependent (BCR and WNN only)
    # TODO: replace with map or map2
    scores_feat <- map_dfr(names(seurat_objs), function(embedding_type) {
      meta_res_type <- meta_res[[embedding_type]]
      obj <- seurat_objs[[embedding_type]]
      obj$adt_feature <-
        case_when(obj@assays$ADT@data[adt_feat, ] > adt_cutoff ~ TRUE, TRUE ~ FALSE)

      reduction_combinations <-
        list(list(reduction = "bcr.umap",
                  labels_pred = meta_res_type[["bcr.umap"]]),
             list(reduction = "rna.umap",
                  labels_pred = meta_res_type[["rna.umap"]]),
             list(reduction = "wnn.umap",
                  labels_pred = meta_res_type[["wnn.umap"]]))

      # TODO: replace with map or map2
      if (metric_type == "Internal") {
        map_dfr(reduction_combinations, function(combo) {
          calc_int_metrics(seurat_obj = obj, reduction_name = combo$reduction,
                           criteria = metrics, labels_true = "adt_feature")
        })
      } else {
        map_dfr(reduction_combinations, function(combo) {
          calc_ext_metrics(seurat_obj = obj, reduction_name = combo$reduction,
                           criteria = metrics, labels_true = "adt_feature",
                           labels_pred = combo$labels_pred)
        })
      }
    })

    # combine results and sort by score
    scores_feat <-
      scores_feat %>%
      mutate(method_combo = paste(Embedding, Reduction, sep = " + "),
             Feature = adt_feat,
             Reduction = factor(Reduction, levels = reductions)) %>%
      arrange(desc(Score))

    scores <- bind_rows(scores, scores_feat)
  }

  return(scores)
}


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
