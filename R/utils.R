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
      cli::cli_abort"Please check that all of the specified metadata columns are present in the Seurat object.")
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
    cli::cli_abort"Neighbors slot '", nn_name, "' not found in object. Available: ",
         paste(names(seurat_obj@neighbors), collapse = ", "))
  }

  seurat_obj@neighbors[[nn_name]]@nn.idx
}
