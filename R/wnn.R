#' Run Weighted Nearest Neighbors (WNN) analysis on combined GEX and BCR data
#'
#' @description
#' This function takes a Seurat object containing gene expression (GEX) data and a matrix of BCR embeddings, and performs WNN analysis to integrate the two modalities.
#' It processes each assay separately (normalization, variable feature selection, scaling, PCA, neighbor finding, UMAP), then finds multimodal neighbors and runs clustering if specified.
#' The function also adds metadata about modality weights and run information to the resulting Seurat object.
#'
#' @details
#' * Currently only works for the embeddings approach and BCR data.
#' * The GEX object must have a `cell_id` metadata column.
#' * If I end up combining multiple embeddings into one object, then I should use
#'   something like `bcr_assay_name <- "BCR" # paste0("BCR_", embedding_type)`
#' * The neighbors for the assays can be saved in both the `graphs` slot (`compute.SNN` for clustering) and the `neighbors` slot (`return.neighbors` for distance calculations later)
#' The `compute.SNN` option constructs a shared nearest neighbor graph using Jaccard index.
#' Perhaps the GEX and BCR sections should be run if `modality_weights` is provided.
#'
#' @param seurat_obj A Seurat object containing GEX data (at the least). If this
#'   is already a merged object produced by [merge_gex_bcr()] (i.e., it contains
#'   a BCR assay and a `bpca` reduction), `embeddings` and `embedding_type` may
#'   be omitted and the existing BCR infrastructure will be reused.
#' @param embeddings Matrix of BCR embeddings (genes by cells). Optional if
#'   `seurat_obj` already contains a BCR assay and `bpca` reduction.
#' @param embedding_type The embeddings method. Optional for merged objects;
#'   falls back to `seurat_obj@misc$embedding_type` if not provided.
#' @param pc_gex The number of PCs for the GEX assay.
#' @param pc_bcr The number of PCs for the BCR assay.
#' @param k The number of neighbors to use for each modality.
#' @param cluster Whether or not to perform clustering.
#' @param cluster_res Named list of clustering resolutions for GEX, BCR, and WNN.
#' @param modality_weights Named vector of modality weights. If NULL, Seurat will calculate automatically.
#' @param verbose Whether or not to show output from Seurat functions.
#'
#' @returns A Seurat object with WNN run.
#' @export
run_wnn <- function(seurat_obj, embeddings, embedding_type, pc_gex = 20,
                    pc_bcr = 20, k_param = 20, cluster = FALSE,
                    cluster_res = list("GEX" = 1, "BCR" = 1, "WNN" = 1),
                    modality_weights = NULL, verbose = FALSE) {
  # TODO: update this to be able to run on other omics e.g. GEX & ADT
  # TODO: add the option to filter out genes again???
  # TODO: give the option to use an integrated GEX assay

  # input validation
  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("seurat_obj must be a Seurat object.")
  }
  if (!"cell_id" %in% colnames(seurat_obj[[]])) {
    cli::cli_abort("Cell ID column not found in metadata.")
  }
  if (missing(pc_gex)) {
    pc_gex <- if ("rpca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["rpca"]])
    } else if ("pca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["pca"]])
    } else {
      20
    }
    cli::cli_inform(c("i" = "Using pc_gex = {pc_gex}", " from existing reductions."))
  }
  if (missing(pc_bcr)) {
    pc_bcr <- if ("bpca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["bpca"]])
    } else {
      20
    }
    cli::cli_inform(c("i" = "Using pc_bcr = {pc_bcr}", " from existing reductions."))
  }
  if (missing(k_param)) {
    # use the k from the GEX neighbors if it exists, otherwise default to 20
    if ("RNA.nn" %in% names(seurat_obj@neighbors)) {
      nn <- seurat_obj@neighbors[["RNA.nn"]]
      k_param <- ncol(nn@nn.idx)
      cli::cli_inform(c("i" = "Using k = {k_param}", " from RNA neighbors."))
    } else {
      k_param <- 20
      cli::cli_inform(c("i" = "Using default k = {k_param}."))
    }
  }

  # detect if the object is already merged (has BCR assay + bpca from merge_gex_bcr)
  is_merged <- "BCR" %in% names(seurat_obj@assays) &&
    "bpca" %in% names(seurat_obj@reductions)

  if (is_merged) {
    cli::cli_inform(c("i" = "Merged object detected: using existing GEX and BCR assays and reductions."))
    if (missing(embedding_type)) {
      embedding_type <- seurat_obj@misc$embedding_type # should exist
    }
  } else {
    if (missing(embeddings)) {
      cli::cli_abort(c(
        "Must provide {.arg embeddings} or pass a merged object.",
        "i" = "Run {.fn merge_gex_bcr} first, or provide an embeddings matrix."
      ))
    }
    if (ncol(embeddings) == 0) {
      cli::cli_abort("No cells found in the provided embeddings matrix.")
    }
  }

  # remove any clustering that might exist and reset the factor levels
  # assumes a UMAP might have been run for the RNA data (instead of tSNE)
  cols_to_remove <- c()
  if (any(grepl("^snn_res", names(seurat_obj[[]])))) {
    cols_to_remove <-
      c(cols_to_remove,
        grep("^snn_res", names(seurat_obj[[]]), value = TRUE))
  }
  if ("seurat_clusters" %in% names(seurat_obj[[]])) {
    cols_to_remove <- c(cols_to_remove, "seurat_clusters")
  }
  seurat_obj@meta.data <-
    seurat_obj[[]] %>%
    {if (length(cols_to_remove) > 0)
      select(., -all_of(cols_to_remove)) else .} %>%
    droplevels()

  if (!is_merged) {
    # TODO: double check this part

    # build BCR object (subset to cells present in the GEX object)
    bcr_obj <- bcr_embeddings_pipeline(
      embeddings[, intersect(colnames(embeddings), seurat_obj$cell_id)],
      embedding_type = embedding_type,
      num_pcs = pc_bcr, num_dims = pc_bcr, k_param = k_param, verbose = verbose)

    # add BCR assay; Seurat v5 does not propagate nCount/nFeature automatically
    suppressWarnings(seurat_obj[["BCR"]] <- bcr_obj[["BCR"]])
    seurat_obj$nCount_BCR <- bcr_obj[[]][Cells(seurat_obj), "nCount_BCR"]
    seurat_obj$nFeature_BCR <- bcr_obj[[]][Cells(seurat_obj), "nFeature_BCR"]
    if ("nFeature_ADT" %in% names(seurat_obj[[]])) {
      seurat_obj@meta.data <-
        seurat_obj[[]] %>%
        relocate(nCount_BCR, nFeature_BCR, .after = nFeature_ADT)
    }

    if (ncol(seurat_obj@assays$RNA) != ncol(seurat_obj@assays$BCR)) {
      cli::cli_abort("The number of cells in the RNA and BCR assays do not match.")
    }

    # transfer reductions and graphs
    seurat_obj[["bpca"]] <- bcr_obj[["bpca"]]
    seurat_obj[["bcr.umap"]] <- bcr_obj[["bcr.umap"]]
    seurat_obj@neighbors[["BCR.nn"]] <- bcr_obj@neighbors[["BCR.nn"]]
    seurat_obj@graphs[["BCR_nn"]] <- bcr_obj@graphs[["BCR_nn"]]
    seurat_obj@graphs[["BCR_snn"]] <- bcr_obj@graphs[["BCR_snn"]]

    # add run info
    Misc(seurat_obj, slot = "embedding_type") <- embedding_type
    Misc(seurat_obj, slot = "embedding_dims") <- nrow(seurat_obj@assays[["BCR"]])
  }

  # find multimodal neighbors, then do clustering and make a UMAP
  seurat_obj <-
    Seurat::FindMultiModalNeighbors(object = seurat_obj,
                                    reduction.list = list("rpca", "bpca"),
                                    dims.list = list(1:pc_gex, 1:pc_bcr),
                                    k.nn = k_param,
                                    # match the RNA and BCR style
                                    knn.graph.name = "w_nn",
                                    snn.graph.name = "w_snn",
                                    weighted.nn.name = "w.nn",
                                    modality.weight.name =
                                      str_c(c("RNA", "BCR"), ".weight"),
                                    return.intermediate = TRUE,
                                    modality.weight = modality_weights,
                                    verbose = verbose)

  # check for NA values
  # TODO: add this check to other UMAPs
  na_nn <- is.na(seurat_obj$RNA.weight)
  if (any(na_nn)) {
    bad_cells <- seurat_obj[[]]$cell_id[na_nn]
    bcr_mat <- GetAssayData(seurat_obj, assay = "BCR", layer = "data")
    bcr_mat <- as.matrix(bcr_mat[, bad_cells])
    # add 1 because the first column is not included in the duplicate count
    dup_count <- sum(duplicated(t(bcr_mat))) + 1

    # TODO: record all duplicates, not just the ones that failed??
    Misc(seurat_obj, slot = "embedding_dups") <- bad_cells

    cli::cli_inform(c("x" = "{sum(na_nn)} NA rows found in w.nn neighbor \\
                            indices (out of {ncol(seurat_obj)} total cells). \\
                            This means that these cells do not have any valid \\
                            neighbors in the WNN space.",
                      "i" = "This is most likely due to these cells having \\
                             identical BCR embeddings. Out of the {sum(na_nn)} \\
                             cells with NA neighbors, {dup_count} have \\
                             identical BCR embeddings.",
                      ">" = "UMAP will not be run for the WNN reduction because
                      of the NA values."))
  } else {
    seurat_obj <- Seurat::RunUMAP(object = seurat_obj, nn.name = "w.nn",
                                  n.neighbors = k_param, # might not be needed
                                  reduction.name = "wnn.umap",
                                  reduction.key = "wnnUMAP_",
                                  verbose = verbose)
  }

  # the Leiden algorithm (4) has been shown to be better than Louvain (1), but
  # in order to use it you have to install the 'leidenbase' package
  if (cluster) {
    algo <- 1

    # cluster the BCR assay
    seurat_obj <- FindClusters(object = seurat_obj,
                               graph.name = "BCR_snn",
                               resolution = cluster_res[["BCR"]],
                               algorithm = algo, verbose = verbose)
    # cluster the GEX assay
    # TODO: do this in seurat_pipeline??
    seurat_obj <- FindClusters(object = seurat_obj,
                               graph.name = "RNA_snn",
                               resolution = cluster_res[["GEX"]],
                               algorithm = algo, verbose = verbose)
    # cluster the WNN assay
    seurat_obj <- FindClusters(object = seurat_obj,
                               graph.name = "w_snn",
                               resolution = cluster_res[["WNN"]],
                               algorithm = algo, verbose = verbose)

    # set the cluster identities and fix the order (RNA and BCR are fine)
    meta_res_wnn <- paste0("w_snn_res.", cluster_res[["WNN"]])
    seurat_obj[[]][[meta_res_wnn]] <- fct_inseq(seurat_obj[[]][[meta_res_wnn]])
    Idents(seurat_obj) <- meta_res_wnn
  }

  # add a metadata column listing which assay was chosen for each cell
  seurat_obj@meta.data <-
    seurat_obj[[]] %>%
    mutate(weight_assay =
             case_when(seurat_obj[[]][["RNA.weight"]] > 0.5 ~ "RNA",
                       seurat_obj[[]][["RNA.weight"]] < 0.5 ~ "BCR",
                       is.na(seurat_obj[[]][["RNA.weight"]]) ~ NA,
                       .default = "Tie"))

  if (verbose) {
    if (!any(na_nn)) {
      cli::cli_inform(c("v" = "WNN neighbors calculated and UMAP run.",
                        "i" = "Use {.fn Seurat::DimPlot} with \\
                               {.arg reduction = 'wnn.umap'} or \\
                               {.fn athanor::plot_dimplot} with \\
                               {.arg reduc = 'wnn.umap'} to visualize the results."))
    }

    # explain the weighting approach
    if (is.null(modality_weights)) {
      cli::cli_inform(c("v" = "Modality weights were automatically calculated based on the provided assays.",
                        "i" = "Check {.fn Seurat::FindMultiModalNeighbors} for more details."))
      # summarize the resulting weight distribution with a total count
      # weight_summary <- seurat_obj[[]] %>%
      #   group_by(weight_assay) %>%
      #   summarise(count = n()) %>%
      #   ungroup()
      # cli::cli_inform(c("v" = "Summary of modality weights across cells:",
      #                   "i" = "{.code print(weight_summary)}"))
    } else {
      cli::cli_inform(c("v" = "Custom modality weights were used:",
                        "i" = "{.arg modality_weights}"))
    }
  }

  return(seurat_obj)
}


#' Give a summary of a Seurat object post-WNN
#'
#' @description
#' This function generates a summary message about the post-WNN Seurat object, including the number of cells, details about the assays (e.g. number of genes, markers, embedding dimensions), information about the reductions used for WNN, and the number of clusters identified in each modality (RNA, BCR, and WNN) based on the largest resolutions.
#'
#' @details
#' Assumes that embeddings were used (for now) and that the object has RNA.
#'
#' @param seurat_obj The post-WNN Seurat object
#' @param gex_pca The name of the GEX PCA reduction.
#' @param other_pca The name of the BCR/ADT/etc. PCA reduction.
#' @param other_type The second assay.
#'
#' @returns A text message.
#' @export
extract_wnn_vars <- function(seurat_obj, gex_pca = "rpca",
                             other_pca = "bpca", other_type = "BCR") {
  # TODO: update this to work on a concatenated object too
  # base message (after AIRR integration)
  message <- paste("This object has", ncol(seurat_obj), "cells, ")

  # assay details
  assay_details <- c()
  bcr_embeddings_name <- grep("^BCR", names(seurat_obj@assays), value = TRUE)

  # GEX info
  if ("RNA" %in% names(seurat_obj@assays)) {
    assay_details <-
      c(assay_details, paste(nrow(seurat_obj@assays$RNA), "genes"))
  }

  # ADT info
  if ("ADT" %in% names(seurat_obj@assays)) {
    assay_details <-
      c(assay_details,
        paste(nrow(seurat_obj@assays$ADT), "cell surface protein markers"))
  }

  # BCR info
  if (length(bcr_embeddings_name) > 0) {
    assay_details <-
      c(assay_details,
        paste(nrow(seurat_obj@assays[[bcr_embeddings_name]]),
              "embedding dimensions"))
  }

  message <- paste0(message, str_c(assay_details, collapse = ", "), ".")

  # reductions info
  if ("weighted.nn" %in% names(seurat_obj@neighbors)) {
    message <- paste(message, "WNN was run with",
                     ncol(seurat_obj@reductions[[gex_pca]]), "GEX PCs and",
                     ncol(seurat_obj@reductions[[other_pca]]), other_type, "PCs.")
  }

  # TODO: remove this?
  # clustering info
  # does not account for multiple clustering resolutions
  if ("seurat_clusters" %in% colnames(seurat_obj[[]])) {
    nclusters_rna <-
      select(seurat_obj[[]], starts_with("RNA_snn")) %>% select(last_col())
    nclusters_bcr <-
      select(seurat_obj[[]], starts_with("BCR_snn")) %>% select(last_col())
    nclusters_wnn <-
      select(seurat_obj[[]], starts_with("w_snn")) %>% select(last_col())

    message <- paste(message,
                     "There were", n_distinct(nclusters_rna), "RNA clusters,",
                     n_distinct(nclusters_bcr), "BCR clusters, and",
                     n_distinct(nclusters_wnn),
                     "WNN clusters identified (in the largest resolutions).")
    # cat(str_glue("There were {nclusters_rna} RNA clusters, {nclusters_bcr} BCR clusters, and {nclusters_wnn} WNN clusters identified."))
  }

  # output the message
  cli::cli_inform(message)
}
