#' Regenerate neighbor graphs and UMAPs.
#'
#' @description
#' Regenerates PCA, neighbor graphs, and UMAP for a Seurat object.
#' This is useful if you subset a Seurat object and need to recompute these reductions for the new object.
#' It can also be used to regenerate these reductions if you have modified the object in a way that requires them to be redone (e.g. filtering out cells or genes).
#'
#' @details
#'
#' @param seurat_obj The Seurat object.
#' @param pca_name Name of the PCA reduction to use for neighbor finding and UMAP. This should be the name of an existing PCA reduction in the Seurat object (e.g. "pca" or "rpca").
#' @param assay Name of the assay to use for neighbor finding and UMAP.
#' @param num_dims Number of PCA dimensions to use for neighbor finding.
#' @param k_param Number of nearest neighbors.
#' @param verbose Print out Seurat's progress messages.
#'
#' @returns A processed Seurat object with the `graphs`, `neighbors`, and `reductions` slots filled in or updated.
#' @export
regen_reduc <- function(seurat_obj, pca_name = "rpca", assay = "RNA",
                        num_dims = 20, k_param = 20, verbose = TRUE) {
  # TODO: regenerate PCA too?
  # TODO: integrate this into other functions e.g. seurat_pipeline, run_wnn?

  # parameter checks
  if (!pca_name %in% names(seurat_obj@reductions)) {
    cli::cli_abort("{pca_name} is not a valid PCA name. Please select one of: {names(seurat_obj@reductions)}.")
  }

  # if assay was not provided, then fall back to the assay that was used to
  # originally compute the PCA reduction
  # note that this could cause undesired behavior e.g. new reductions made by batch effect integration methods will typically use a different name
  if (missing(assay)) {
    assay <- seurat_obj@reductions[[pca_name]]@assay.used
    cli::cli_inform("Assay not provided, so using assay {assay} that was used to compute the PCA reduction.")
  }

  # if dims and k are not provided, get them from the object
  nn_name <- stringr::str_c(assay, ".nn")
  nn_command <- stringr::str_c("FindNeighbors", assay, pca_name, sep = ".")
  if (missing(num_dims) | missing(k_param)) {
    # try to use the existing neighbors slot if possible
    if (stringr::str_c(assay, ".nn") %in% names(seurat_obj@neighbors)) {
      nn <- seurat_obj@neighbors[[nn_name]]

      if (missing(num_dims)) num_dims <- nn@alg.info$ndim
      if (missing(k_param)) k_param <- ncol(nn@nn.idx)

      cli::cli_inform("Using existing neighbor graph for assay {assay} to determine num_dims ({num_dims}) and k_param ({k_param}).")
    } else if (nn_command %in% names(seurat_obj@commands)) {
      cmd <- seurat_obj@commands[[nn_command]]

      if (missing(num_dims)) num_dims <- length(cmd$dims)
      if (missing(k_param)) k_param <- cmd$k.param

      cli::cli_inform("Using existing command {nn_command} for assay {assay} to determine num_dims ({num_dims}) and k_param ({k_param}).")
    }
      else {
      cli::cli_warn("No existing neighbor graph found for assay {assay}, \\
                    so using default values for num_dims (20) and k_param (20).")
      if (missing(num_dims)) num_dims <- 20
      if (missing(k_param)) k_param <- 20
    }
  }

  # fill in the "graphs" slot
  seurat_obj <- FindNeighbors(seurat_obj, reduction = pca_name,
                              dims = 1:num_dims, k.param = k_param,
                              graph.name =
                                stringr::str_c(assay, "_", c("", "s"), "nn"),
                              verbose = verbose)
  # fill in the "neighbors" slot
  seurat_obj <- FindNeighbors(seurat_obj, reduction = pca_name,
                              dims = 1:num_dims, k.param = k_param,
                              return.neighbor = TRUE, graph.name = nn_name,
                              verbose = verbose)
  # fill in the "reductions" slot
  assay <- tolower(assay)
  seurat_obj <- RunUMAP(seurat_obj, reduction = pca_name,
                        n.neighbors = k_param, nn.name = nn_name,
                        reduction.name = stringr::str_c(assay, ".umap"),
                        reduction.key = stringr::str_c(assay, "MAP_"),
                        verbose = verbose)

  seurat_obj
}
