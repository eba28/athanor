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
#'
#' @param seurat_obj Processed Seurat object to reduce.
#' @param dim_reducs Vector of dimensionality reductions to keep.
#' @param meta_cols Vector of metadata column names to keep. If unspecified, keeps all metadata columns.
#' @param remove_neighbors Whether or not to remove neighbor graphs from the object.
#' @param print_size Whether to print info about how much the object was reduced.
#'
#' @returns A reduced Seurat object with specified reductions kept.
#' @export
reduce_object <- function(seurat_obj, dim_reducs = "umap", meta_cols,
                          remove_neighbors = TRUE, print_size = TRUE, ...) {
  cat(paste("Currently reducing:", deparse(substitute(seurat_obj)), "\n"))

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
    cat(paste("Original object size:",
              format(object.size(seurat_obj), units = "auto"), "\n"))
    cat(paste("Reduced object size:",
              format(object.size(obj_reduced), units = "auto"), "\n"))
  }

  return(obj_reduced)
}
