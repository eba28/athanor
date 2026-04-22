#' Compute per-cell mean distances to nearest neighbors
#'
#' @param adt_mat A cells by features numeric matrix.
#' @param nn_idx An integer matrix of nearest-neighbor indices (cells by k).
#' @param distance_metric One of `"euclidean"`, `"manhattan"`, `"mean_abs"`,
#'   or `"cosine"`.
#' @param exclude_self Logical; whether to drop a cell from its own neighbor
#'   list before computing distances.
#' @param return_mean Logical; if `TRUE` returns a single grand mean instead of
#'   a per-cell vector.
#'
#' @returns A named numeric vector of per-cell mean distances, or a scalar if
#'   `return_mean = TRUE`.
#' @keywords internal
adt_dists_core <- function(adt_mat, nn_idx, distance_metric = "euclidean",
                           exclude_self = TRUE, return_mean = FALSE) {
  k_dist <- vapply(seq_len(nrow(adt_mat)), function(cell_i) {
    idx <- nn_idx[cell_i, ]

    if (exclude_self) {
      idx <- idx[idx != cell_i]
      if (length(idx) == 0) return(NA_real_)
    }

    xi <- adt_mat[cell_i, , drop = FALSE]
    nbrs <- adt_mat[idx, , drop = FALSE]
    diffs <- sweep(nbrs, 2, as.numeric(xi), FUN = "-")

    per_nbr <- if (distance_metric == "mean_abs") {
      rowMeans(abs(diffs))
    } else if (distance_metric == "manhattan") {
      rowSums(abs(diffs))
    } else if (distance_metric == "euclidean") {
      sqrt(rowSums(diffs * diffs))
    } else if (distance_metric == "cosine") {
      xi_vec <- as.numeric(xi)
      vapply(seq_len(nrow(nbrs)), function(j) {
        lsa::cosine(xi_vec, as.numeric(nbrs[j, ]))
      }, numeric(1))
    } else {
      stop("Please choose a valid distance metric: mean_abs, manhattan, ",
           "euclidean, cosine.")
    }

    mean(per_nbr, na.rm = TRUE)
  }, FUN.VALUE = numeric(1))

  if (return_mean) {
    mean(k_dist, na.rm = TRUE)
  } else {
    names(k_dist) <- rownames(adt_mat)
    k_dist
  }
}

#' Extract and validate an ADT matrix from a Seurat object
#'
#' @details
#' The Seurat object must contain an assay with the specified `adt_assay` name, and that assay must contain a layer with the specified `layer` name.
#' If `features_adt` is provided, it must be a character vector of feature names that are present in the assay.
#' The function will return a matrix of dimension cells by features, where the rows are named with cell IDs and the columns are named with feature names.
#'
#' @param seurat_obj A Seurat object.
#' @param adt_assay Name of the ADT assay.
#' @param layer Layer to retrieve from the assay.
#' @param features_adt Optional character vector of ADT features to subset.
#'   If `NULL`, all features are returned.
#'
#' @returns A numeric matrix of dimension cells by features.
#' @keywords internal
get_adt_matrix <- function(seurat_obj, adt_assay = "ADT", layer = "data",
                           features_adt = NULL) {
  # check the arguments
  if (!adt_assay %in% names(seurat_obj@assays)) {
    stop("Assay '", adt_assay, "' not found in Seurat object. Available assays: ",
         paste(names(seurat_obj@assays), collapse = ", "))
  }
  if (!layer %in% names(seurat_obj@assays[[adt_assay]]@data)) {
    stop("Layer '", layer, "' not found in assay '", adt_assay, "'. Available layers: ",
         paste(names(seurat_obj@assays[[adt_assay]]@data), collapse = ", "))
  }

  mat <- GetAssayData(seurat_obj, assay = adt_assay, layer = layer)

  if (!is.null(features_adt)) {
    keep <- intersect(features_adt, rownames(mat))
    if (length(keep) == 0) {
      stop("None of the requested ADT features are present in assay '",
           adt_assay, "'.")
    }
    mat <- mat[keep, , drop = FALSE]
  }

  t(as.matrix(mat))
}
