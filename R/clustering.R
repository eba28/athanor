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
