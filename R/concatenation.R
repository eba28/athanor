#' Concatenate GEX and BCR data in a Seurat object
#'
#' @description
#' Creates a new assay combining gene expression (GEX) and B-cell receptor (BCR)
#' data by concatenating processed BCR features with RNA data. Runs the standard
#' Seurat pipeline on the combined data.
#'
#' @details
#' The function:
#'   1. Extracts specified BCR metadata columns from the Seurat object
#'   2. Processes BCR features using `process_bcr_features()`
#'   3. Creates a new assay by row-binding RNA and BCR data
#'   4. Optionally filters out IG/TR genes from variable features
#'   5. Runs standard Seurat workflow: normalize, scale, PCA, neighbors, UMAP
#'
#' @note
#' Currently assumes BCR data is already integrated into object metadata.
#' See `integrate_gex_airr()` for adding AIRR data to a Seurat object.
#'
#' @param seurat_obj A Seurat object containing RNA assay and BCR metadata.
#' @param pca_stage Add BCR information before PCA or after PCA.
#' @param cols_to_include Character vector of BCR metadata column names to include
#'   in the concatenated assay e.g. c("mu_freq", "isotype") or embedded dimensions.
#' @param var_features Logical; if TRUE, run FindVariableFeatures on the combined
#'   assay. If FALSE, concatenate BCR features onto existing variable features.
#' @param normalize Logical; if TRUE, normalize the combined assay using
#'   LogNormalize. If FALSE, skip normalization.
#' @param num_dims Integer; number of dimensions to use for PCA and neighbor finding.
#' @param filter_genes Logical; if TRUE, filter out immunoglobulin (IG) and
#'   T-cell receptor (TR) genes from variable features. Requires `remove_genes`
#'   to be defined in the environment.
#'
#' @return
#' A Seurat object with:
#'   - New `RNA_BCR` assay containing concatenated GEX and BCR features
#'   - PCA reduction (`rna_bcr.pca`)
#'   - UMAP reduction (`rna_bcr.umap`)
#'   - Neighbor graphs computed on the combined data
#' @export
concatenate_gex_bcr <- function(seurat_obj, pca_stage = "Before",
                                cols_to_include, var_features = FALSE,
                                normalize = TRUE, num_dims = 20, # k
                                filter_genes = TRUE) {
  # check that we have everything we need
  if (filter_genes) {
    if (!exists("remove_genes")) {
      stop("Make sure that `remove_genes` has been defined.")
    }
  }

  if (pca_stage == "Before") {
    # TODO: check if the RNA assay contains the names of the BCR assay
    # intersect(cols_to_include, rownames(GetAssayData(seurat_obj, assay = "RNA")))

    # select features
    bcr_features <- seurat_obj[[]] %>% select(all_of(cols_to_include))

    # TODO: add an option to do weighting to influence the effect of the BCRs (post-normalization)

    # process BCR features (rename, convert, normalize)
    # TODO: check if the embeddings output has already been normalized by amulety
    bcr_features <- process_bcr_features(bcr_features)
    # bcr_features[, 1:5] # check them visually

    # make a new assay with both genes and BCR data
    # TODO: try making a BCR assay first?

    # need counts if we're going to normalize later
    seurat_layer <- ifelse(normalize, "counts", "data")
    seurat_obj@assays$RNA_BCR <-
      CreateAssayObject(counts = rbind(GetAssayData(seurat_obj, assay = "RNA",
                                                    layer = seurat_layer),
                                       bcr_features))
    DefaultAssay(seurat_obj) <- "RNA_BCR"

    # make sure that the BCR data are included as variable features
    # TODO: try other selection methods (not just vst)
    if (var_features) {
      seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
      ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

      ## Warning: Layer counts isn't present in the assay object; returning NULL
      ## Warning message:
      ## In FindVariableFeatures.Assay(object = object[[assay]], selection.method = selection.method,  :
      ##   selection.method set to 'vst' but count slot is empty; will use data slot instead
    } else {
      # just concatenate onto the existing ones
      # which already had IG features removed, so it's less than 2000
      VariableFeatures(seurat_obj) <-
        c(VariableFeatures(seurat_obj, assay = "RNA"), rownames(bcr_features))
      VariableFeatures(seurat_obj) <- unique(VariableFeatures(seurat_obj))
    }
    # TODO: run the finder, then concatenate anyways?

    ##   Warning messages:
    ## 1: In eval(predvars, data, env) : NaNs produced
    ## 2: In hvf.info$variance.expected[not.const] <- 10^fit$fitted :
    ##   number of items to replace is not a multiple of replacement length

    # TODO: add the option to not filter out the IG genes
    # so I'd have to rerun FindVariableFeatures on the RNA alone

    # filter out the IG and TR genes
    # (assumes that features_meta is already loaded and remove_genes defined)
    # see the "Features to be filtered out" section
    if (filter_genes) {
      remove_feats <- VariableFeatures(seurat_obj) %in% remove_genes
      VariableFeatures(seurat_obj) <- VariableFeatures(seurat_obj)[!remove_feats]

      cat(paste0("After removing IG/TR genes, the total number of variable features is: ",
                 length(VariableFeatures(seurat_obj)), ". Of these, ",
                 length(setdiff(VariableFeatures(seurat_obj),
                                str_replace_all(rownames(bcr_features), "_", "-"))),
                 " are in GEX alone.\n"))
    }

    # process with the standard Seurat pipeline
    # TODO: update seurat_pipeline() to take an alternate reduction name?
    # seurat_obj <- seurat_pipeline(seurat_obj, num_dims = k,
    #                               filter_genes = filter_genes, verbose = FALSE)

    # don't normalize twice
    if (normalize) {
      seurat_obj <- NormalizeData(seurat_obj,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000, verbose = FALSE)
    }

    # TODO: don't scale twice??
    seurat_obj <- ScaleData(object = seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(object = seurat_obj, npcs = 30,
                         reduction.name = "rna_bcr.pca",
                         reduction.key = "rnabcrpca_")
  } else {
    # check that the necessary PCAs exist
    if (!"rpca" %in% names(seurat_obj@reductions)) {
      stop("Make sure that the Seurat object contains an RNA PCA.")
    }
    if (!"bpca" %in% names(seurat_obj@reductions)) {
      stop("Make sure that the Seurat object contains a BCR PCA.")
    }

    # remove unnecessary WNN information
    seurat_obj@graphs$w_nn <- c()
    seurat_obj@graphs$w_snn <- c()
    seurat_obj@neighbors$w.nn <- c()
    seurat_obj@reductions$wnn.umap <- c()

    # just use the RNA slot by default

    # update the number of dimension
    num_dims <- num_dims * 2

    # assumes that both reductions have the same number of PCs
    combined_pca <-
      cbind(Embeddings(seurat_obj, "rpca"), Embeddings(seurat_obj, "bpca"))
    colnames(combined_pca) <- str_c("pca_", 1:num_dims) # necessary? what about cell ids?
    # what about the feature loadings?

    # new PCA reduction
    seurat_obj[["rna_bcr.pca"]] <-
      CreateDimReducObject(embeddings = combined_pca, key = "rnabcrpca_",
                           assay = DefaultAssay(seurat_obj))

  }

  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rna_bcr.pca",
                              k.param = num_dims, dims = 1:num_dims,
                              verbose = FALSE, graph.name = "RNA_BCR_nn")
  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rna_bcr.pca",
                              k.param = num_dims, dims = 1:num_dims,
                              return.neighbor = TRUE,
                              verbose = FALSE, graph.name = "RNA_BCR.nn")
  seurat_obj <- RunUMAP(object = seurat_obj, # n.neighbors = num_dims,
                        nn.name = "RNA_BCR.nn",
                        reduction = "rna_bcr.pca", # dims = 1:num_dims,
                        reduction.name = "rna_bcr.umap") # , verbose = FALSE
  ## 16:48:25 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30

  return(seurat_obj)
}
