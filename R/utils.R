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


#' Map user-facing assay names to internal neighbor slot prefixes
#'
#' @param base_assay Character string; one of `"GEX"`, `"WNN"`, or an assay
#'   name passed through unchanged.
#'
#' @returns A character scalar prefix used to look up the neighbor slot.
#' @keywords internal
map_assay_name <- function(base_assay) {
  switch(base_assay, GEX = "RNA", WNN = "w", base_assay)
}


#' Reduce a Seurat object's size
#'
#' @description
#' This function reduces a Seurat object by removing count matrices and keeping only
#' specified dimensionality reductions. This is especially useful for creating
#' lightweight objects for Shiny apps or sharing.
#'
#' @details
#' Modify this as needed if your object is built differently (e.g. tSNE instead).
#' Uses `DietSeurat` to remove counts while preserving reductions and metadata.
#' This is especially useful if you are making a Shiny app or just visualizing the data.
#' This works well in conjunction with [seurat_pipeline()].
#'
#' @param seurat_obj Processed Seurat object to reduce.
#' @param dim_reducs Vector of dimensionality reductions to keep.
#' @param meta_cols Vector of metadata column names to keep. If unspecified, keeps all metadata columns.
#' @param remove_neighbors Whether or not to remove neighbor graphs from the object.
#' @param print_size Whether to print info about how much the object was reduced.
#'
#' @returns A reduced Seurat object with specified reductions kept.
#' @export
reduce_object <- function(seurat_obj, dim_reducs = "rna.umap", meta_cols,
                          remove_neighbors = TRUE, print_size = TRUE, ...) {
  cli::cli_inform("Currently reducing: {deparse(substitute(seurat_obj))}")

  # modify this as desired
  obj_reduced <- DietSeurat(object = seurat_obj, dimreducs = dim_reducs, ...)

  # the metadata can take up a lot of memory, so you can filter it down to just the columns you need
  if (!rlang::is_missing(meta_cols)) {
    if (!all(meta_cols %in% colnames(obj_reduced[[]]))) {
      stop("Please check that all of the specified metadata columns are present in the Seurat object.")
    }

    obj_reduced@meta.data <- obj_reduced@meta.data %>% select(all_of(meta_cols))
  }

  # you might not need the graphs here if you're just using the object for plotting
  if (remove_neighbors) {
    obj_reduced@neighbors <- list()
  }

  # print the before and after sizes
  if (print_size) {
    cli::cli_inform(c(
      "i" = "Original object size: {format(object.size(seurat_obj), units = 'auto')}",
      "i" = "Reduced object size: {format(object.size(obj_reduced), units = 'auto')}"
    ))
  }

  return(obj_reduced)
}


#' Resolve the nearest-neighbor index matrix from a Seurat object
#'
#' @param seurat_obj A Seurat object containing a neighbors slot.
#' @param base_assay Character string passed to `map_assay_name()` to determine
#'   which neighbor slot to retrieve.
#'
#' @returns An integer matrix of nearest-neighbor indices (cells by k).
#' @keywords internal
resolve_neighbors <- function(seurat_obj, base_assay) {
  prefix <- map_assay_name(base_assay)
  nn_name <- paste0(prefix, ".nn")

  if (!nn_name %in% names(seurat_obj@neighbors)) {
    stop("Neighbors slot '", nn_name, "' not found in object. Available: ",
         paste(names(seurat_obj@neighbors), collapse = ", "))
  }

  seurat_obj@neighbors[[nn_name]]@nn.idx
}