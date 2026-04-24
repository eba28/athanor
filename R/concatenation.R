#' Concatenate GEX and BCR data in a Seurat object
#'
#' @description
#' Creates a new assay combining gene expression (GEX) and B-cell receptor (BCR)
#' data by concatenating processed BCR features with RNA data. Runs the standard
#' Seurat pipeline on the combined data.
#'
#' @details
#' This would typically be used after [seurat_pipeline()] and [gex_add_airr()].
#' The function:
#'   1. Extracts specified BCR metadata columns from the Seurat object
#'   2. Processes BCR features using [process_bcr_features()]
#'   3. Creates a new assay by row-binding RNA and BCR data
#'   4. Optionally filters out IG/TR genes from variable features
#'   5. Runs standard Seurat workflow: normalize, scale, PCA, neighbors, UMAP
#'
#' @note
#' Currently assumes BCR data is already integrated into object metadata.
#'
#' @param seurat_obj A Seurat object containing RNA assay and BCR metadata.
#' @param pca_stage Add BCR information before PCA or after PCA.
#' @param cols_to_include Character vector of BCR metadata column names to include
#'   in the concatenated assay e.g. c("mu_freq", "isotype") or embedded dimensions.
#' @param var_features If TRUE, run FindVariableFeatures on the combined
#'   assay. If FALSE, concatenate BCR features onto existing variable features.
#' @param normalize If TRUE, normalize the combined assay using
#'   LogNormalize. If FALSE, skip normalization.
#' @param num_pcs Number of principal components to compute.
#' @param num_dims Number of PCA dimensions to use for neighbor finding and UMAP.
#' @param k_param Number of nearest neighbors.
#' @param filter_genes If specified, filter out genes from this category (e.g. "IG" and/or "TR")
#' @param ensembl_version If filtering genes, specify the Ensembl version to use for gene annotations (e.g. "GRCh38.104"). If NULL, uses the default version in [get_airr_genes()].
#'
#' @return
#' A Seurat object with:
#'   - New `RNA_BCR` assay containing concatenated GEX and BCR features
#'   - PCA reduction (`rna_bcr.pca`)
#'   - UMAP reduction (`rna_bcr.umap`)
#'   - Neighbor graphs computed on the combined data
#' @export
concatenate_gex_bcr <- function(seurat_obj, pca_stage = c("Before", "After"),
                                cols_to_include, var_features = FALSE,
                                normalize = TRUE, num_pcs = 50, num_dims = 20,
                                k_param = 20, filter_genes,
                                ensembl_version = NULL) {
  # TODO: include k as a parameter?
  pca_stage <- match.arg(pca_stage)
  # TODO: double check the formatting of filter genes
  # TODO: give an option to use a pre-defined list of filter_genes

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
    seurat_obj[["RNA_BCR"]] <-
      CreateAssayObject(counts = rbind(GetAssayData(seurat_obj, assay = "RNA",
                                                    layer = seurat_layer),
                                       bcr_features))
    DefaultAssay(seurat_obj) <- "RNA_BCR"

    # TODO: use seurat_pipeline()

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

    # filter out the IG and/or TR genes
    # TODO: check the object for the species, ensembl version
    if (!rlang::is_missing(filter_genes)) {
      seurat_obj <- filter_variable_features(seurat_obj, filter_genes,
                                             ensembl_version = ensembl_version,
                                             bcr_features = bcr_features)
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
    seurat_obj <- RunPCA(object = seurat_obj, npcs = num_pcs,
                         reduction.name = "rna_bcr.pca",
                         reduction.key = "rnabcrpca_")
  } else {
    # check that the necessary PCAs exist
    if (!"rpca" %in% names(seurat_obj@reductions)) {
      cli::cli_abort(c(
        "{.arg seurat_obj} must contain an RNA PCA reduction.",
        i = "Run {.fn RunPCA} with {.code reduction.name = \"rpca\"} first."
      ))
    }
    if (!"bpca" %in% names(seurat_obj@reductions)) {
      cli::cli_abort(c(
        "{.arg seurat_obj} must contain a BCR PCA reduction.",
        i = "Run {.fn RunPCA} with {.code reduction.name = \"bpca\"} first."
      ))
    }

    # remove unnecessary WNN information
    # seurat_obj@graphs$w_nn <- c()
    # seurat_obj@graphs$w_snn <- c()
    # seurat_obj@neighbors$w.nn <- c()
    # seurat_obj@reductions$wnn.umap <- c()
    seurat_obj[["w_nn"]] <- NULL
    seurat_obj[["w_snn"]] <- NULL
    seurat_obj[["w.nn"]] <- NULL
    seurat_obj[["wnn.umap"]] <- NULL

    # just use the RNA slot by default

    # assumes that both reductions have the same number of PCs
    # TODO: address how in the combined PCA has total_dims columns but FindNeighbors uses dims = 1:num_dims, which only covers half the combined dimensions
    total_dims <- num_dims * 2
    combined_pca <-
      cbind(Embeddings(seurat_obj, "rpca"), Embeddings(seurat_obj, "bpca"))
    colnames(combined_pca) <- str_c("pca_", 1:total_dims) # TODO: check if this is necessary and what happens to cell ids
    # what about the feature loadings?

    # new PCA reduction
    seurat_obj[["rna_bcr.pca"]] <-
      CreateDimReducObject(embeddings = combined_pca, key = "rnabcrpca_",
                           assay = DefaultAssay(seurat_obj))

  }

  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rna_bcr.pca",
                              k.param = k_param, dims = 1:num_dims,
                              verbose = FALSE, graph.name = "RNA_BCR_nn")
  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rna_bcr.pca",
                              k.param = k_param, dims = 1:num_dims,
                              return.neighbor = TRUE,
                              verbose = FALSE, graph.name = "RNA_BCR.nn")
  seurat_obj <- RunUMAP(object = seurat_obj, nn.name = "RNA_BCR.nn",
                        n.neighbors = k_param, reduction.name = "rna_bcr.umap",
                        reduction.key = "rnabcrUMAP_")

  return(seurat_obj)
}
