#' Print an overview of a Seurat object.
#'
#' @description
#' Prints a structured summary of a Seurat object, covering its assays,
#' reductions, neighbor graphs, graphs, WNN configuration (if present), and
#' Misc slot contents.
#'
#' @details
#' Each section is only printed if the corresponding slot is non-empty, so the
#' output degrades gracefully on minimal objects.
#' For RNA, the number of variable features is shown.
#' For ADT, all marker names are listed.
#' The WNN section is printed only when a `w.nn` neighbor graph is detected and the
#' `FindMultiModalNeighbors` command is recorded in the object.
#' Misc slot values are printed inline for short atomic vectors and summarized by type
#' for larger or complex objects.
#'
#' @param seurat_obj The Seurat object.
#'
#' @returns A overview of the Seurat object, and invisibly returns the input object for piping if desired.
#' @export
object_overview <- function(seurat_obj) {
  # --- Cells & Assays ---
  cli::cli_h2("Cells & Assays")
  cli::cli_inform("{ncol(seurat_obj)} cells")

  for (assay_name in names(seurat_obj@assays)) {
    assay_obj <- seurat_obj@assays[[assay_name]]
    n_feats <- nrow(assay_obj)

    if (assay_name == "RNA") {
      n_var <- length(VariableFeatures(seurat_obj))
      cli::cli_inform("RNA: {n_feats} genes, {n_var} variable features")
    } else if (assay_name == "ADT") {
      markers <- rownames(assay_obj)
      cli::cli_inform("ADT: {n_feats} markers ({toString(stringr::str_sort(markers, numeric = TRUE))})")
    } else {
      cli::cli_inform("{assay_name}: {n_feats} features")
    }
  }

  # --- Reductions ---
  if (length(seurat_obj@reductions) > 0) {
    cli::cli_h2("Reductions")
    for (reduc_name in names(seurat_obj@reductions)) {
      n_dims <- ncol(seurat_obj@reductions[[reduc_name]])
      pc_info <- ""
      if (grepl("umap", tolower(reduc_name))) {
        for (cmd_name in names(seurat_obj@commands)) {
          if (!startsWith(cmd_name, "RunUMAP")) next
          cmd <- seurat_obj@commands[[cmd_name]]
          if (!identical(cmd@params$reduction.name, reduc_name)) next
          nn_ref <- cmd@params$nn.name
          if (!is.null(nn_ref) && nn_ref %in% names(seurat_obj@neighbors)) {
            nn_dims <- seurat_obj@neighbors[[nn_ref]]@alg.info$ndim
            if (!is.null(nn_dims)) pc_info <- paste0(" (", nn_dims, " PCs)")
          } else if (!is.null(cmd@params$dims)) {
            nn_dims <- length(cmd@params$dims)
            pc_info <- paste0(" (", nn_dims, " PCs)")
          }
          break
        }
      }
      cli::cli_inform("{reduc_name}: {n_dims} dimensions{pc_info}")
    }
  }

  # --- Neighbors ---
  if (length(seurat_obj@neighbors) > 0) {
    cli::cli_h2("Neighbors")
    for (nn_name in names(seurat_obj@neighbors)) {
      nn <- seurat_obj@neighbors[[nn_name]]
      k <- ncol(nn@nn.idx)
      n_dims <- nn@alg.info$ndim
      if (!is.null(n_dims)) {
        cli::cli_inform("{nn_name}: k = {k}, {n_dims} PCs")
      } else {
        cli::cli_inform("{nn_name}: k = {k}")
      }
    }
  }

  # --- Graphs ---
  if (length(seurat_obj@graphs) > 0) {
    cli::cli_h2("Graphs")
    fnn_cmds <- Filter(function(nm) startsWith(nm, "FindNeighbors"),
                       names(seurat_obj@commands))
    graph_dims <- NULL
    for (cmd_name in fnn_cmds) {
      cmd <- seurat_obj@commands[[cmd_name]]
      if (!is.null(cmd@params$dims) && !isTRUE(cmd@params$return.neighbor)) {
        graph_dims <- length(cmd@params$dims)
        break
      }
    }
    if (!is.null(graph_dims)) {
      cli::cli_inform("{toString(names(seurat_obj@graphs))} ({graph_dims} PCs)")
    } else {
      cli::cli_inform("{toString(names(seurat_obj@graphs))}")
    }
  }

  # --- WNN (if present) ---
  wnn_cmd_name <- "Seurat..FindMultiModalNeighbors"
  if ("w.nn" %in% names(seurat_obj@neighbors) &&
      wnn_cmd_name %in% names(seurat_obj@commands)) {
    cli::cli_h2("WNN")
    dims_list <- seurat_obj@commands[[wnn_cmd_name]]@params$dims.list
    reduction_list <- seurat_obj@commands[[wnn_cmd_name]]@params$reduction.list
    for (i in seq_along(dims_list)) {
      reduc <- if (!is.null(reduction_list)) reduction_list[[i]] else paste("modality", i)
      n_dims <- length(dims_list[[i]])
      cli::cli_inform("{reduc}: {n_dims} dimensions")
    }
  }

  # --- Misc ---
  misc <- Seurat::Misc(seurat_obj)
  if (length(misc) > 0) {
    cli::cli_h2("Misc")
    for (slot_name in names(misc)) {
      val <- misc[[slot_name]]
      if (is.atomic(val) && length(val) <= 10) {
        cli::cli_inform("{slot_name}: {toString(val)}")
      } else if (is.atomic(val)) {
        cli::cli_inform("{slot_name}: [{length(val)}-length {class(val)} vector]")
      } else {
        cli::cli_inform("{slot_name}: [{class(val)[1]}]")
      }
    }
  }

  invisible(seurat_obj)
}


#' Regenerate neighbor graphs and UMAPs.
#'
#' @description
#' Regenerates PCA, neighbor graphs, and UMAP for a Seurat object.
#' This is useful if you subset a Seurat object and need to recompute these reductions for the new object.
#' It can also be used to regenerate these reductions if you have modified the object in a way that requires them to be redone (e.g. filtering out cells or genes).
#'
#' @details
#' If you are providing a reduction based on batch effect integration (e.g. Harmony, RPCA), then you should use the name of that reduction for `pca_name`.
#' Although this function is called "regen_reduc", it can also be run to generate reductions for the first time (e.g. if you have already run PCA and just want to calculate neighbors and run UMAP).
#'
#' @param seurat_obj The Seurat object.
#' @param pca_name Name of the PCA reduction to use for neighbor finding and UMAP. This should be the name of an existing PCA reduction in the Seurat object (e.g. "pca" or "rpca").
#' @param assay Name of the assay to use for neighbor finding and UMAP.
#' @param num_dims Number of PCA dimensions to use for neighbor finding.
#' @param k_param Number of nearest neighbors.
#' @param verbose Logical indicating whether or not to print messages.
#'
#' @returns A processed Seurat object with the `graphs`, `neighbors`, and `reductions` slots filled in or updated.
#' @export
regen_reduc <- function(seurat_obj, pca_name = "rpca", assay = "RNA",
                        num_dims = 20, k_param = 20, verbose = TRUE) {
  # TODO: regenerate PCA too?
  # TODO: integrate this into other functions e.g. seurat_pipeline, run_wnn?
  # TODO: don't use the cli calls if verbose = FALSE
  # TODO: just reuse seurat_pipeline and bcr_embeddings_pipeline???

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
  # TODO: don't use assay to search for a command because Harmony and other batch effect integration methods will typically use a different assay name for their neighbor finding step
  nn_command <- stringr::str_c("FindNeighbors", assay, pca_name, sep = ".")

  # TODO: only print out the actual values that are filled in
  if (missing(num_dims) || missing(k_param)) {
    # try to use the existing neighbors slot if possible
    if (nn_name %in% names(seurat_obj@neighbors)) {
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
      cli::cli_inform(c("i" = "No existing neighbor graph found for assay {assay}, \\
                               so using default values for num_dims (20) and k_param (20)."))
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
  # the function won't use assay since we are providing the nn graph, but we provide it
  # so that the command is recorded properly (it will use the default assay otherwise)
  umap_name <- paste0(tolower(assay), ".umap")
  umap_key <- paste0(gsub("_", "", tolower(assay)), "UMAP_")
  seurat_obj <- RunUMAP(seurat_obj, reduction = pca_name, assay = assay,
                        nn.name = nn_name, n.neighbors = k_param,
                        reduction.name = umap_name, reduction.key = umap_key,
                        verbose = verbose)

  seurat_obj
}


#' Run Seurat's standard pipeline
#'
#' @description
#' Runs the standard Seurat pipeline (normalize, scale, PCA, neighbors, UMAP) on a
#' Seurat object. Optionally filters cells by QC metrics first, filters IG/TR genes
#' from variable features, and clusters. Supports any assay (e.g. `"RNA"`,
#' `"RNA_BCR"`), with reduction names derived automatically from the assay name.
#'
#' @details
#' Cell filtering (`nfeatures_RNA`, `perc_mt`) and ADT normalization are only
#' applied when `assay = "RNA"`. For other assays, pass `normalize = FALSE` if
#' the data has already been normalized.
#'
#' @param seurat_obj The Seurat object.
#' @param assay Name of the assay to run the pipeline on. Defaults to `"RNA"`.
#' @param pca_name Name to give the PCA reduction. If `NULL`, defaults to
#'   `"rpca"` for the RNA assay, or `paste0(tolower(assay), ".pca")` for others
#'   (e.g. `"rna_bcr.pca"` for `assay = "RNA_BCR"`).
#' @param nfeatures_RNA Minimum number of RNA features to retain per cell. If
#'   omitted, cell filtering is skipped. Only applies when `assay = "RNA"`.
#' @param perc_mt Maximum percentage of mitochondrial genes to retain. If
#'   omitted, cell filtering is skipped. Only applies when `assay = "RNA"`.
#' @param num_features Desired number of variable features.
#' @param num_pcs Number of principal components to compute.
#' @param num_dims Number of PCA dimensions to use for neighbor finding.
#' @param k_param Number of nearest neighbors.
#' @param normalize If `TRUE`, normalize the assay using LogNormalize before
#'   scaling. Set to `FALSE` if the assay has already been normalized.
#' @param find_var_features If `TRUE`, run [FindVariableFeatures()] before
#'   scaling. Set to `FALSE` if variable features are already set on the assay
#'   (e.g. when calling this after manually setting them in [concatenate_gex_bcr()]).
#' @param cluster_res Clustering resolution(s). If `NULL`, clustering is skipped.
#' @param filter_genes If specified, filter out genes from this category (e.g. `"IG"` and/or `"TR"`).
#' @param ensembl_version Ensembl version for gene annotations (e.g. `"GRCh38.104"`). If `NULL`, uses the default in [get_airr_genes()].
#' @param cache_file Passed to [get_airr_genes()]. Path to a cached RDS result to use instead of querying Ensembl.
#' @param verbose Logical indicating whether or not to print messages.
#'
#' @returns A processed Seurat object with PCA, neighbor graphs, optional
#'   clusters, and UMAP. Reduction names are derived from `assay` and `pca_name`.
#' @export
seurat_pipeline <- function(seurat_obj, assay = "RNA", pca_name = NULL,
                            nfeatures_RNA, perc_mt,
                            num_features = 2000, num_pcs = 50, num_dims = 20,
                            k_param = 20, normalize = TRUE,
                            find_var_features = TRUE, cluster_res = NULL,
                            filter_genes, ensembl_version = NULL,
                            cache_file = NULL, verbose = TRUE) {
  # derive reduction and graph names from assay
  if (is.null(pca_name)) {
    pca_name <- if (assay == "RNA") "rpca" else paste0(tolower(assay), ".pca")
  }
  pca_key <- paste0(gsub("[._]", "", tolower(pca_name)), "_")
  snn_name <- paste0(assay, "_snn")

  # cell filtering - RNA only
  if (assay == "RNA" && !rlang::is_missing(nfeatures_RNA) && !rlang::is_missing(perc_mt)) {
    if ("percent.mt" %in% names(seurat_obj[[]])) {
      seurat_obj <- subset(seurat_obj,
                           subset = nFeature_RNA > nfeatures_RNA & percent.mt < perc_mt)
    } else {
      warning("No filtration was performed upon this object.")
    }
  }

  # normalization
  if (normalize) {
    seurat_obj <- NormalizeData(seurat_obj, assay = assay,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000, verbose = verbose)
    # ADT normalization - RNA only, normalizes across cells not features
    if (assay == "RNA" && "ADT" %in% names(seurat_obj@assays)) {
      seurat_obj <- NormalizeData(seurat_obj,
                                  normalization.method = "CLR", margin = 2,
                                  assay = "ADT", verbose = verbose)
    }
  }

  # highly variable features
  if (find_var_features) {
    seurat_obj <- FindVariableFeatures(seurat_obj, assay = assay,
                                       selection.method = "vst",
                                       nfeatures = num_features, verbose = verbose)
  }

  # note: `features = rownames(seurat_obj)` can cause crashes
  seurat_obj <- ScaleData(seurat_obj, assay = assay, verbose = verbose)

  if (!rlang::is_missing(filter_genes)) {
    seurat_obj <- filter_variable_features(seurat_obj, filter_genes,
                                           ensembl_version = ensembl_version,
                                           cache_file = cache_file)
  }

  # irlba throws a warning to "use a standard svd instead" when requesting more
  # than 50% of all singular values, so use exact SVD in that case (also faster
  # when the embedding dimension is small)
  scale_data <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "scale.data")
  max_dim <- min(nrow(scale_data), ncol(scale_data))

  # ScaleData silently no-ops when VariableFeatures names don't match assay
  # feature names (Seurat 5 rewrites "_" to "-" at assay creation, so a
  # pre-rename VariableFeatures list finds nothing to scale)
  if (nrow(scale_data) == 0L) {
    cli::cli_abort(c(
      "ScaleData produced an empty {.code scale.data} layer for assay {.val {assay}}.",
      "i" = "Usually means {.fn VariableFeatures} contains names that don't \\
             match the assay's actual feature names.",
      "i" = "Seurat 5 rewrites {.code _} to {.code -} in feature names at \\
             assay creation; pre-rename BCR features (e.g. \\
             {.code gsub('_', '.', x)}) if your VariableFeatures list \\
             contains underscores."
    ))
  }

  # Cap num_pcs at matrix rank - 1 (RunPCA errors or hangs otherwise; common
  # in concatenate_gex_bcr "reduced" where combined_mat has ~26 rows)
  if (num_pcs > max_dim - 1L) {
    if (verbose) {
      cli::cli_warn("num_pcs ({num_pcs}) exceeds matrix rank ({max_dim}); \\
                     capping at {max_dim - 1L}.")
    }
    num_pcs <- max_dim - 1L
  }

  # FindNeighbors errors if dims > computed PCs
  if (num_dims > num_pcs) {
    if (verbose) {
      cli::cli_inform(c("i" = "num_dims ({num_dims}) > num_pcs ({num_pcs}); \\
                                capping num_dims to {num_pcs}."))
    }
    num_dims <- num_pcs
  }

  use_approx <- num_pcs < max_dim / 2
  seurat_obj <- RunPCA(seurat_obj, assay = assay, npcs = num_pcs,
                       reduction.name = pca_name, reduction.key = pca_key,
                       approx = use_approx, verbose = verbose)
  cli::cli_inform(c("v" = "Computed PCA with {num_pcs} dimensions using \\
                    {ifelse(use_approx, 'approximate', 'exact')} SVD."))

  # SNN graph for clustering + neighbor object for UMAP and evaluation
  seurat_obj <- regen_reduc(seurat_obj, pca_name = pca_name, assay = assay,
                            num_dims = num_dims, k_param = k_param,
                            verbose = verbose)

  if (!is.null(cluster_res)) {
    seurat_obj <- FindClusters(seurat_obj, graph.name = snn_name,
                               resolution = cluster_res, verbose = verbose)

    # fix the cluster levels (they sort alphabetically by default)
    for (res in cluster_res) {
      res_col <- paste0(snn_name, "_res.", res)
      seurat_obj[[res_col]] <-
        factor(seurat_obj[[]][[res_col]],
               str_sort(unique(seurat_obj[[]][[res_col]]), numeric = TRUE))
    }
  }

  seurat_obj
}
