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
#' @param seurat_obj A Seurat object containing GEX data (at the least).
#' @param embeddings Matrix of BCR embeddings (genes x cells).
#' @param embedding_type The embeddings method.
#' @param pc_gex The number of PCs for the GEX assay.
#' @param pc_bcr The number of PCs for the BCR assay.
#' @param k_param The number of neighbors to use for each modality. Can be a single value or a vector of values to test.
#' @param k_main The main number of neighbors to use for the final WNN UMAP and clustering.
#' @param cluster Whether or not to perform clustering.
#' @param cluster_res Named list of clustering resolutions for GEX, BCR, and WNN.
#' @param modality_weights Named vector of modality weights. If NULL, Seurat will calculate automatically.
#' @param show_output Whether or not to show verbose output from Seurat functions.
#'
#' @returns A Seurat object with WNN run.
#' @export
run_wnn <- function(seurat_obj, embeddings, embedding_type, pc_gex = 20,
                    pc_bcr = 20, k_param = 20, k_main = 20, cluster = FALSE,
                    cluster_res = list("GEX" = 1, "BCR" = 1, "WNN" = 1),
                    modality_weights = NULL, show_output = FALSE) {
  # input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!"cell_id" %in% colnames(seurat_obj[[]])) {
    stop("Cell ID column not found in metadata.")
  }
  if (ncol(embeddings) == 0) {
    stop("No cells found in embeddings matrix.")
  }

  # add k to the names only if there is more than one
  add_k <- ifelse(length(k_param) > 1, TRUE, FALSE)

  # remove any GEX clustering that might exist and reset the factor levels
  # assumes a UMAP might have been run for the RNA data (instead of tSNE)
  cols_to_remove <- c()
  if (any(grepl("^RNA_snn_res", names(seurat_obj[[]])))) {
    cols_to_remove <-
      c(cols_to_remove,
        grep("^RNA_snn_res", names(seurat_obj[[]]), value = TRUE))
  }
  if ("seurat_clusters" %in% names(seurat_obj[[]])) {
    cols_to_remove <- c(cols_to_remove, "seurat_clusters")
  }
  seurat_obj@meta.data <- seurat_obj[[]] %>%
                            {if (length(cols_to_remove) > 0)
                            select(., -all_of(cols_to_remove)) else .} %>%
                          droplevels()
  # keep wiping the slate clean
  Idents(seurat_obj) <- "orig.ident"
  seurat_obj@graphs <- list() # RNA.nn, RNA.snn
  seurat_obj@reductions <- list() # pca, umap

  # create and format the BCR assay
  bcr_assay <- CreateAssayObject(counts = embeddings) # CreateAssay5Object
  # subset down to cells present in both assays
  bcr_assay <- subset(bcr_assay, cells = seurat_obj$cell_id)

  # add BCR data as a new assay to the GEX object
  suppressWarnings(seurat_obj[["BCR"]] <- bcr_assay)
  # this will add nCount_BCR and nFeature_BCR to the metadata

  # change the order of the metadata just for consistency
  seurat_obj@meta.data <-
    seurat_obj[[]] %>%
    relocate(nCount_BCR, nFeature_BCR, .after = nFeature_ADT)

  # make sure that the number of cells match
  if (ncol(seurat_obj@assays$RNA) != ncol(seurat_obj@assays$BCR)) {
    stop("The number of cells in the RNA and BCR assays do not match.")
  }

  # BCR processing
  DefaultAssay(seurat_obj) <- "BCR"
  VariableFeatures(seurat_obj) <- rownames(seurat_obj[["BCR"]])
  # seurat_obj <- NormalizeData(object = seurat_obj)
  seurat_obj <- ScaleData(object = seurat_obj, verbose = show_output)
  seurat_obj <- RunPCA(object = seurat_obj,
                       npcs = pc_bcr, reduction.name = "bpca",
                       reduction.key = "bpca_", verbose = show_output)
  for (k in k_param) {
    # fill the neighbors slot
    neighbor_name <- ifelse(add_k, paste0("BCR.nn_", k), "BCR.nn")
    seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "bpca",
                                assay = "BCR", k.param = k,
                                return.neighbor = TRUE, verbose = show_output,
                                graph.name = neighbor_name)

    # fill the graphs slot
    graph_name <- str_c("BCR_", c("", "s"), "nn")
    if (add_k) graph_name <- str_c(graph_name, "_", k_main)
    seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "bpca",
                                dims = 1:pc_bcr, k.param = k_main,
                                compute.SNN = TRUE, verbose = show_output,
                                graph.name = graph_name)
  }
  nn_name <- ifelse(add_k, paste0("BCR.nn_", k_main), "BCR.nn")
  seurat_obj <- RunUMAP(object = seurat_obj, # dims = 1:pc_bcr,
                        reduction = "bpca",
                        n.neighbors = k_main, nn.name = nn_name,
                        reduction.name = "bcr.umap", reduction.key = "bcrUMAP_",
                        verbose = show_output)

  # RNA processing (a lot of this was already done)
  # TODO: use seurat_pipeline()
  # TODO: filter out genes
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(object = seurat_obj, verbose = show_output)
  seurat_obj <- FindVariableFeatures(object = seurat_obj, verbose = show_output)
  seurat_obj <- ScaleData(object = seurat_obj, verbose = show_output)
  seurat_obj <- RunPCA(object = seurat_obj,
                       npcs = pc_gex, reduction.name = "rpca",
                       reduction.key = "rpca_", verbose = show_output)
  for (k in k_param) {
    # fill the neighbors slot
    neighbor_name <- ifelse(add_k, paste0("RNA.nn_", k), "RNA.nn")
    seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rpca",
                                assay = "RNA", k.param = k,
                                return.neighbor = TRUE, verbose = show_output,
                                graph.name = neighbor_name)

    # fill the graphs slot
    graph_name <- str_c("RNA_", c("", "s"), "nn")
    if (add_k) graph_name <- str_c(graph_name, "_", k_main)
    seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rpca",
                                dims = 1:pc_gex, k.param = k_main,
                                compute.SNN = TRUE, verbose = show_output,
                                graph.name = graph_name)

  }
  nn_name <- ifelse(add_k, paste0("RNA.nn_", k_main), "RNA.nn")
  seurat_obj <- RunUMAP(object = seurat_obj, # dims = 1:pc_gex,
                        reduction = "rpca",
                        n.neighbors = k_main, nn.name = nn_name,
                        reduction.name = "rna.umap", reduction.key = "rnaUMAP_",
                        verbose = show_output)

  # find multimodal neighbors, then do clustering and make a UMAP
  for (k in k_param) {
    mw_name <- str_c(c("RNA", "BCR"), ".weight")
    if (add_k) mw_name <- str_c(mw_name, "_", k)
    seurat_obj <-
      FindMultiModalNeighbors(object = seurat_obj,
                              reduction.list = list("rpca", "bpca"),
                              dims.list = list(1:pc_gex, 1:pc_bcr),
                              k.nn = k,
                              # match the RNA and BCR style
                              knn.graph.name =
                                ifelse(add_k, paste0("w_nn_", k), "w_nn"),
                              snn.graph.name =
                                ifelse(add_k, paste0("w_snn_", k), "w_snn"),
                              weighted.nn.name =
                                ifelse(add_k, paste0("w.nn_", k), "w.nn"),
                              modality.weight.name = mw_name,
                              return.intermediate = TRUE,
                              modality.weight = modality_weights,
                              verbose = show_output)
  }
  nn_name <- ifelse(add_k, paste0("w.nn_", k_main), "w.nn")
  seurat_obj <- RunUMAP(object = seurat_obj, # dims = 1:pc_gex,
                        n.neighbors = k_main, nn.name = nn_name,
                        reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                        verbose = show_output)

  # use the "main" k for clustering
  if (cluster) {
    # cluster the BCR assay
    graph_name <- ifelse(add_k, paste0("BCR_snn_", k_main), "BCR_snn")
    seurat_obj <- FindClusters(object = seurat_obj,
                               graph.name = graph_name,
                               resolution = cluster_res[["BCR"]],
                               algorithm = 1, verbose = show_output)

    # cluster the GEX assay
    graph_name <- ifelse(add_k, paste0("RNA_snn_", k_main), "RNA_snn")
    seurat_obj <- FindClusters(object = seurat_obj,
                               graph.name = graph_name,
                               resolution = cluster_res[["GEX"]],
                               algorithm = 1, verbose = show_output)

    # cluster the WNN assay
    graph_name <- ifelse(add_k, paste0("w_snn_", k_main), "w_snn")
    seurat_obj <- FindClusters(object = seurat_obj,
                               graph.name = graph_name,
                               resolution = cluster_res[["WNN"]],
                               algorithm = 1, verbose = show_output)

    # set the cluster identities and fix the order (RNA and BCR are fine)
    meta_res_wnn <- paste0("w_snn_res.", cluster_res[["WNN"]])
    seurat_obj[[]][[meta_res_wnn]] <- fct_inseq(seurat_obj[[]][[meta_res_wnn]])
    Idents(seurat_obj) <- meta_res_wnn
  }

  # add a metadata column listing which assay was chosen for each cell
  # only for the main k
  weight_col <- ifelse(add_k, paste0("RNA.weight_", k_main), "RNA.weight")
  seurat_obj@meta.data <-
    seurat_obj[[]] %>%
    mutate(weight_assay =
             case_when(seurat_obj[[]][[weight_col]] > 0.5 ~ "RNA",
                       seurat_obj[[]][[weight_col]] < 0.5 ~ "BCR",
                       .default = "Tie"))

  # add run info
  Misc(seurat_obj, slot = "embeddings_dims") <- nrow(seurat_obj@assays[["BCR"]])
  Misc(seurat_obj, slot = "embedding_type") <- embedding_type
  Misc(seurat_obj, slot = "default_k") <- k_main

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
#' @param other_type The second assay. Defaults to "BCR".
#'
#' @returns A text message.
#' @export
extract_wnn_vars <- function(seurat_obj, gex_pca = "rpca",
                             other_pca = "bpca", other_type = "BCR") {
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
  cat(message)
}
