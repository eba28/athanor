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
