#' Manually simulate gene expression data
#'
#' @description
#' This function simulates gene expression data by generating a matrix of counts using a Poisson distribution, where the rate parameter is set to 0.5.
#' The resulting matrix is then converted to a sparse matrix format and formatted with gene and cell names.
#' Finally, a Seurat object is created from the counts matrix, and metadata is added with cell identifiers.
#'
#' @details
#' Partially based off of [Single Cell workshop: Chapter 6](http://gateway.training.ncgr.org/single-cell-workshop/seurat.html)
#'
#' @param num_genes Number of genes to simulate.
#' @param num_cells Number of cells to simulate.
#' @param separator Separator for gene and cell names.
#'
#' @returns A Seurat object
#' @export
sim_gex_manual <- function(num_genes = 1000, num_cells = 2000,
                           separator = "-") {
  # set up the counts matrix
  raw_counts <- as.integer(rexp(num_genes * num_cells, rate = 0.5))
  raw_counts <- Matrix(data = raw_counts,
                       nrow = num_genes, ncol = num_cells,
                       sparse = TRUE)

  # format the genes and cells
  rownames(raw_counts) <- paste("Gene", 1:num_genes, sep = separator)
  colnames(raw_counts) <- paste("Cell", 1:num_cells, sep = separator)

  # create the object
  manual <- CreateSeuratObject(counts = raw_counts, project = "Manual")

  # add metadata
  manual$cell_id <- Cells(manual)

  return(manual)
}


#' Simulate gene expression data using Splatter
#'
#' @description
#' This function simulates gene expression data using the Splatter package, which generates synthetic single-cell RNA-seq data based on a specified set of parameters.
#' The function allows for customization of the number of genes, cells, simulation method, group probabilities, and random seed for reproducibility.
#' The resulting simulated data is then converted into a Seurat object with appropriate metadata and cell identifiers.
#'
#' @details
#' Splatter doesn't use a separator for its fake names.
#'
#' @param num_genes Number of genes to simulate. Defaults to 1000.
#' @param num_cells Number of cells to simulate. Defaults to 2000.
#' @param splatter_method Splatter simulation method. Defaults to "single".
#' @param splatter_groups Group probabilities for Splatter simulation. Defaults to 1.
#' @param seed Random seed for reproducible simulations. Defaults to 42.
#' @param verbose Whether or not to print verbose output. Defaults to FALSE.
#'
#' @returns A Seurat object containing the simulated gene expression data with
#'   metadata and cell identifiers
#' @export
sim_gex_splatter <- function(num_genes = 1000, num_cells = 2000,
                             splatter_groups = 1, splatter_method = "single",
                             seed = 42, verbose = FALSE) {
  # splatter with defaults - gives a SCE object
  sim <- splatSimulate(params = newSplatParams(nGenes = num_genes,
                                               batchCells = num_cells,
                                               group.prob = splatter_groups,
                                               seed = seed),
                       method = splatter_method,
                       verbose = verbose)

  # we don't need most of what's in the object e.g. the rowData
  # (Gene, BaseGeneMean, OutlierFactor, GeneMean, DEFacGroup1, DEFacGroup2) or
  # most of the colData (Cell, Batch, Group, ExpLibSize)
  sim <- minimiseSCE(sim, colData.keep = "Cell", metadata.keep = TRUE,
                     verbose = verbose)

  # convert to Seurat
  seurat_obj <-
    CreateSeuratObject(counts = sim@assays@data@listData[["counts"]],
                       project = "Splatter", assay = "RNA",
                       meta.data = as.data.frame(sim@colData@listData,
                                                 row.names = "Cell"))
  Misc(seurat_obj, slot = "meta") <- sim@metadata$Params

  # add metadata
  seurat_obj$cell_id <- Cells(seurat_obj)

  return(seurat_obj)
}


#' Simulate an assay made from BCR embeddings independent of any GEX object
#'
#' @description
#' This function simulates an assay of BCR embeddings by generating a matrix of random values between -0.6 and 0.6, with a specified number of cells and embedding dimensions.
#'
#' @details
#' The choice of value ranges was based off of real data, including having 9 decimal points.
#'
#' @param num_cells Number of cells to simulate.
#' @param num_dims Number of embedding dimensions to simulate.
#' @param separator Separator for cell and dimension names.
#'
#' @returns A Seurat Assay
#' @export
sim_bcr_manual <- function(num_cells, num_dims, separator) {
  # simulate a matrix of immune2vec-style embeddings
  bcr_embeddings <- round(runif(n = num_cells * num_dims, -0.6, 0.6), 9)
  bcr_embeddings <- Matrix(data = bcr_embeddings,
                           nrow = num_cells, ncol = num_dims,
                           sparse = TRUE)

  # format the AIRR parameters and cells
  rownames(bcr_embeddings) <- paste("Cell", 1:num_cells, sep = separator)
  colnames(bcr_embeddings) <- paste("Dim", 1:num_dims, sep = separator)

  # match the existing format
  bcr_embeddings <- t(bcr_embeddings)

  return(bcr_embeddings)
}

#' Run manual WNN simulations by varying a specific variable while keeping others constant
#' @export
run_wnn_sims <- function(count_range, sim_var, other_vars,
                         show_progress = FALSE) {
  # setup
  manual_wnn_test <- data.frame()

  if (show_progress) progress_bar <- txtProgressBar(min = min(count_range),
                                                    max = max(count_range),
                                                    char = "=", style = 3)

  # go through each count
  for (i in count_range) {
    no_error <- 1
    other_vars[[sim_var]] <- i

    # show the current progress
    if (show_progress) setTxtProgressBar(progress_bar, value = i)

    tryCatch({
      run_wnn(seurat_obj =
                sim_gex_manual(num_genes = other_vars[["Genes"]],
                               num_cells = other_vars[["Cells"]],
                               separator = "-"),
              embeddings =
                sim_bcr_manual(num_cells = other_vars[["Cells"]],
                               num_dims = other_vars[["Dimensions"]],
                               separator = "-"),
              pc_gex = other_vars[["GEX PCs"]],
              pc_bcr = other_vars[["BCR PCs"]])
    },
    error = function(e) {no_error <<- 0},
    warning = function(w) {no_error <<- 0.5}
    )

    manual_wnn_test <- bind_rows(manual_wnn_test,
                                 c("Count" = i, "Passed" = no_error))
  }

  if (show_progress) close(progress_bar)
  return(manual_wnn_test)
}


#' Plot the results of manual WNN testing
#'
#' @description
#' This function creates a step plot to visualize the results of manual WNN testing, showing how the success of Seurat's WNN computation changes as the count of a specific variable is varied while keeping other variables constant.
#'
#' @param manual_wnn_test A data frame containing the results of the manual WNN testing, with columns "Count" and "Passed".
#' @param sim_var The variable that was varied in the manual WNN testing (e.g. "Genes", "Cells", "Dimensions", "GEX PCs", or "BCR PCs").
#' @param other_vars A named list of the other variables that were held constant in the manual WNN testing, with names corresponding to the variable names (e.g. "Genes", "Cells", "Dimensions", "GEX PCs", or "BCR PCs") and values corresponding to the constant values used in the testing.
#' @param count_range A numeric vector specifying the range of counts that were tested in the manual WNN testing, used for setting the x-axis breaks in the plot.
#'
#' @returns A ggplot object showing the step plot of the WNN testing results, with the x-axis representing the count of the varied variable and the y-axis representing whether Seurat's WNN computation passed (1), failed (0), or gave a warning (0.5). The plot includes a title indicating the variable that was varied and a subtitle listing the other variables that were held constant.
#' @export
plot_wnn_testing <- function(manual_wnn_test, sim_var,
                             other_vars, count_range) {
  # don't include the test var in the subtitle info
  other_vars <- other_vars[names(other_vars) != sim_var]

  # step plot
  ggplot(manual_wnn_test, aes(x = Count, y = Passed)) +
    geom_step(aes(color = Passed), linewidth = 2,
              color = named_colors$sim_vars[[sim_var]]) +
    labs(title = paste("Number of", sim_var, "Needed for Manual WNN"),
         subtitle = paste(other_vars, names(other_vars),
                          sep = " ", collapse = " | "),
         x = paste(sim_var, "Count"), y = "Seurat's WNN Computation") +
    scale_x_continuous(breaks = count_range) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5),
                       labels = c("Failed", "Warning", "Passed")) +
    theme_bw + labels_standard_larger +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "none")
}


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


#' Plot UMAPs of a Seurat object post-WNN
#'
#' @description
#' This function creates a combined plot of the GEX, BCR, and WNN UMAPs from a post-WNN Seurat object, colored by a specified metadata column (e.g. clusters or cell types).
#' The function allows for customization of the plot title, point size, color palette, and whether to display metadata labels on the UMAPs.
#' It uses the `plot_dimplot` function to create individual UMAP plots for each assay and then combines them using `patchwork` for a cohesive visualization.
#'
#' @details
#' Should be able to plot ADT instead of BCR too.
#' Only plots one clustering.
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param data_source The source of the data for the plot title.
#' @param airr_type The type of AIRR data.
#' @param airr_processing The type of AIRR processing; one of `c("Embeddings", "Features)`
#' @param reducs_list Should be in order (GEX, AIRR, WNN).
#' @param clusters_col The metadata column to color the UMAPs by.
#' @param plot_label Whether or not to plot the metadata labels on the UMAPs.
#' @param pt_size The point size for the UMAPs.
#' @param clrs_specific A specific (must have names) color palette for the clusters. If not provided, the default Seurat colors will be used.
#'
#' @returns A combined plot of the GEX, BCR, and WNN UMAPs colored by the specified metadata column.
#' @export
plot_wnn_umaps <- function(seurat_obj, data_source = "Manual",
                           airr_type = "BCR",
                           airr_processing = "Embeddings",
                           reducs_list = c("rna.umap", "bcr.umap", "wnn.umap"),
                           clusters_col = "seurat_clusters",
                           plot_label = TRUE, pt_size = 0.8, clrs_specific) {
  # set the groupings
  Idents(seurat_obj) <- clusters_col

  # could use nlevels(seurat_obj), but maybe the annotations column isn't a factor
  num_clusts <- n_distinct(seurat_obj[[clusters_col]])

  # if you want to use the default Seurat colors
  if (rlang::is_missing(clrs_specific)) clrs_specific <- hue_pal()(num_clusts)

  if (!rlang::is_missing(airr_processing)) {
    airr_assay <- paste0(airr_type, " (", airr_processing, ")")
  }
  else{
    airr_assay <- airr_type
  }

  # reducs_list <- names(seurat_obj@reductions)

  # RNA UMAP
  p1 <- plot_dimplot(seurat_obj = seurat_obj, data_source = data_source,
                     clrs_specific = clrs_specific, pt_size = pt_size,
                     assay = "GEX", reduc = reducs_list[1],
                     plot_label = plot_label, clusters_col = clusters_col)

  # BCR UMAP
  p2 <- plot_dimplot(seurat_obj = seurat_obj, data_source = data_source,
                     clrs_specific = clrs_specific, pt_size = pt_size,
                     assay = airr_assay, reduc = reducs_list[2],
                     plot_label = plot_label, clusters_col = clusters_col)

  # WNN UMAP
  p3 <- plot_dimplot(seurat_obj = seurat_obj, data_source = data_source,
                     clrs_specific = clrs_specific, pt_size = pt_size,
                     assay = paste("GEX &", airr_type, "WNN"), reduc = reducs_list[3],
                     plot_label = plot_label, clusters_col = clusters_col)

  if (clusters_col == "seurat_clusters") {
    p1 <- p1 + labs(color = "WNN")
    p2 <- p2 + labs(color = "WNN")
    p3 <- p3 + labs(color = "WNN")
  } else {
    p1 <- p1 + labs(color = clusters_col)
    p2 <- p2 + labs(color = clusters_col)
    p3 <- p3 + labs(color = clusters_col)
  }

  # all UMAPs
  (p1 | p2 | p3) +
    plot_layout(guides = "collect") & plot_anno &
    guides(color = guide_legend(nrow = 2)) & # optional
    legend_bottom
}


#' Plot a box plot of modality weights per cell type
#'
#' @description
#' This function creates a box plot to visualize the distribution of modality weights (e.g. RNA vs. BCR) across different cell types or clusters in a post-WNN Seurat object.
#'
#' @details
#' Assumes annotated_clusters is a column
#'
#' @param seurat_obj The post-WNN Seurat object.
#' @param details Details to add to the plot title.
#' @param second_assay List of other assays run through WNN in order.
#' @param clrs_specific A specific (must have names) color palette.
#' @param split_by A meta.data column to split the box plots up by.
#' @param y_axis_label Label for the y-axis.
#'
#' @returns A ggplot with the distribution of weights
#' @export
plot_mws <- function(seurat_obj, details = "", second_assay = "BCR",
                     clrs_specific = named_colors$mu_freq_bins,
                     split_by = "mu_freq_bins",
                     y_axis_label = "SHM Frequency Bins") {
  main_assay <- ifelse(length(second_assay) > 1, second_assay[-1], second_assay)

  weight <- grep(paste0("^", main_assay, ".*\\.weight.*$"),
                 colnames(seurat_obj[[]]), value = TRUE)
  # in case of multiple k's
  if (length(weight) > 1) {
    weight <- grep(seurat_obj@misc$default_k, weight, value = TRUE)
  }

  n_assay <- length(second_assay) + 1

  if (!is.factor(seurat_obj[[]][split_by])) {
    seurat_obj[[]][split_by] <- factor(seurat_obj[[]][[split_by]])
  }

  # don't require using the embeddings approach
  if ("embedding_type" %in% names(seurat_obj@misc)) {
    subtitle <- embedding_types[[seurat_obj@misc$embedding_type]]
  } else {
    subtitle <- NULL
  }

  p <- ggplot(seurat_obj[[]],
              aes(x = !!sym(weight), y = !!sym(split_by),
                  fill = !!sym(split_by))) +
         geom_boxplot(outlier.size = 0.5) +
         geom_jitter(size = 0.2) +
         labs(title = paste(details, "Weights by Cell Type"),
              subtitle = subtitle, x = "Weights", y = y_axis_label) +
         scale_fill_manual(values = clrs_specific) +
         facet_wrap(vars(annotated_clusters), scales = "fixed") +
         theme_bw + labels_standard + theme(legend.position = "none")

  if (n_assay == 2) {
      p <- p +
             geom_vline(xintercept = 0.50, linetype = "dashed") +
             scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                                labels = c("0 [GEX]", "0.25", "0.50", "0.75",
                                           paste0("1 [", second_assay, "]")),
                                limits = c(0, 1))
  } else if (n_assay == 3) {
      p <- p +
             geom_vline(xintercept = 1/3, linetype = "dashed") +
             geom_vline(xintercept = 2/3, linetype = "dashed") +
             scale_x_continuous(breaks = seq(0, 1, by = 0.1),
                                labels = c("0.0 [GEX]", "0.1", "0.2",
                                           paste0("0.3 [GEX:", second_assay[1], "]"),
                                           "0.4", "0.5",
                                           paste0("0.6 [", second_assay[1], ":",
                                                  second_assay[2], "]"),
                                           "0.7", "0.8", "0.9",
                                           paste0("1 [", second_assay[2], "]")),
                                limits = c(0, 1))
  }

  p
}
