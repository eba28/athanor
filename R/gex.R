#' Add a user-specified list of cluster annotations to a Seurat object
#'
#' @description
#' This function adds a user-specified list of cluster annotations to the Seurat object.
#' It can be used for adding both manual and automated annotations, updating cell identities
#' and adding annotation columns to the metadata.
#'
#' @details
#' This would typically be used after [seurat_pipeline()].
#' This assumes that cell typing was done on a cluster level, so the annotations_df should have one row per cluster. If you have cell-level annotations, you can skip the relabeling and just add the metadata column to the Seurat object.
#'
#' @param seurat_obj The Seurat object to annotate.
#' @param annotations_df Data frame containing cluster-to-cell-type mappings, typically with "Cluster" and "CellType" columns.
#' @param cell_types_col The name of the column containing the cell type annotations.
#' @param relabel Whether to update the Seurat object's active identities. Sometimes you just want to add the metadata.
#' @param relocate Whether to relocate the annotation column in metadata.
#' @param alphabetize Whether to alphabetize the cell types.
#' @param clusters_col The name of the column containing cluster IDs.
#' @param annotations_col The name of the new metadata column for the annotations.
#'
#' @returns A Seurat object with added annotation information.
#' @export
add_annotations <- function(seurat_obj, annotations_df,
                            cell_types_col = "CellType",
                            relabel = TRUE, relocate = TRUE, alphabetize = TRUE,
                            clusters_col = "seurat_clusters",
                            annotations_col = "annotated_clusters") {
   # prepare the annotation mapping: cluster -> cell type
   annotations <- annotations_df[[cell_types_col]]
   if (!is.factor(seurat_obj[[]][[clusters_col]])) {
      seurat_obj[[clusters_col]] <-
         factor(seurat_obj[[]][[clusters_col]],
                levels = str_sort(unique(seurat_obj[[]][[clusters_col]])))
   }
   names(annotations) <- levels(seurat_obj[[]][[clusters_col]])

   # TODO: simplify the renaming logic
   current_idents <- Idents(seurat_obj)
   Idents(seurat_obj) <- clusters_col
   seurat_obj <- RenameIdents(seurat_obj, annotations)

   if (alphabetize) {
      Idents(seurat_obj) <- factor(Idents(seurat_obj),
                                   levels = sort(levels(seurat_obj)))
   }

   seurat_obj[[annotations_col]] <- Idents(seurat_obj)

   if (!relabel) Idents(seurat_obj) <- current_idents

   if (relocate) {
      seurat_obj@meta.data <-
         seurat_obj[[]] %>%
         relocate(!!sym(annotations_col), .after = !!sym(clusters_col))
   }

   return(seurat_obj)
}


#' Run automated cell type annotation
#'
#' @description
#' This function runs automated cell type annotation using `CellTypist`.
#'
#' @details
#' Supports CellTypist annotation methods.
#' Assumes that the `Cells()` of `seurat_obj` are properly formatted (i.e. unique).
#' For `CellTypist`, assumes the H5AD file and predictions have already been generated.
#' This would typically be used after [seurat_pipeline()].
#'
#' "Majority voting refines the prediction result in a local cell cluster by
#' choosing the dominant cell type label but may increase the runtime especially
#' for a large dataset due to the over-clustering step. This approach usually
#' improves the cell annotation, as voting is conducted in small subclusters
#' derived from over-clustering (cells belonging to a given cell type will be
#' assigned the same label regardless of potential batch effects separating them)."
#' - https://www.celltypist.org/tutorials/onlineguide
#'
#' @param seurat_obj The Seurat object. Must be the path to a H5AD object if using CellTypist.
#' @param annotation_method Which method to use: CellTypist
#' @param reference Reference or model to use for prediction.
#' @param majority_voting Whether to enable majority voting for CellTypist predictions, which refines predictions based on local cluster information but may increase runtime.
#'
#' @returns A data.frame with the annotations for each cell
#' @export
automated_annotation <- function(seurat_obj, annotation_method,
                                 reference = "pbmcref",
                                 majority_voting = FALSE,
                                 over_clustering = NULL) {
   # TODO: Rename seurat_obj since it can be scanpy too

   # validate input
   valid_methods <- c("CellTypist") # for now
   if (!any(annotation_method %in% valid_methods)) {
     cli::cli_abort("Method must be one of: {paste(valid_methods, collapse = ', ')}")
   }

   # run CellTypist annotation
   if (annotation_method == "CellTypist") {
      cli::cli_inform("Processing CellTypist annotation...")

      if (!requireNamespace("reticulate", quietly = TRUE)) {
         cli::cli_abort("The `reticulate` package is required for running \\
                        CellTypist predictions through R.")
      }

      # load model and data
      model <- celltypist$models$Model$load(model = reference)
      obj_h5ad <- scanpy$read_h5ad(filename = seurat_obj)

      # CellTypist requires log1p-normalized counts in .X (not scaled data).
      # If .X contains scaled data (values < 0 or max >> 15), normalize from raw counts.
      if (min(obj_h5ad$X) < 0) {
         cli::cli_inform("Detected scaled data in .X; normalizing from raw counts for CellTypist.")
         if (!is.null(obj_h5ad$raw)) {
            reticulate::py_set_attr(obj_h5ad, "X", obj_h5ad$raw$X)
         } else {
            cli::cli_abort("CellTypist requires log1p-normalized counts in .X, but .X contains \\
                           scaled data and no raw counts slot is available. \\
                           Re-export the H5AD before scaling (e.g. before ScaleData()).")
         }
         scanpy$pp$normalize_total(obj_h5ad, target_sum = 1e4)
         scanpy$pp$log1p(obj_h5ad)
      }

      # run predictions (on a per cell level)
      if (majority_voting) {
         cli::cli_inform("Majority voting enabled: predictions will be refined based on local cluster information, which may increase runtime.")
      } else {
         cli::cli_inform("Majority voting not enabled: predictions will be made on a per cell basis without refinement.")
      }
      predictions <- celltypist$annotate(filename = obj_h5ad, model = model,
                                         majority_voting = majority_voting,
                                         over_clustering = over_clustering)

      # get the predicted cell types and the confidence score
      annotations <- predictions$to_adata()
      if (majority_voting) {
         annotations <- annotations$obs %>%
                           select(cell_id, majority_voting, conf_score) %>%
                           # for consistency with the non-majority voting output, rename the majority_voting column to predicted_labels
                           rename(predicted_labels = majority_voting) %>%
                           remove_rownames()
      } else {
         annotations <- annotations$obs %>%
                           select(cell_id, predicted_labels, conf_score) %>%
                           remove_rownames()
      }

   }

   return(annotations)
}


#' Map cell types to Seurat clusters
#'
#' @description
#' This function maps cell types to Seurat clusters by counting the number of cells in each cluster that belong to each cell type, and then assigning the most common cell type to each cluster.
#' It returns a data frame with the assigned cell types for each cluster.
#'
#' @details
#' This would typically be used after [seurat_pipeline()] and [add_annotations()].
#'
#' @param seurat_obj The Seurat object.
#' @param clusters_col The metadata column with the Seurat clusters.
#' @param annotations_col The metadata column with the cell types.
#'
#' @returns A data.frame with a row for each cell type.
#' @export
cell_type_clusters <- function(seurat_obj, clusters_col = "seurat_clusters",
                               annotations_col) {
   seurat_obj[[]] %>%
      select(!!sym(clusters_col), !!sym(annotations_col)) %>%
      group_by(!!sym(clusters_col), !!sym(annotations_col)) %>%
      tally() %>% # TODO: replace with count
      slice_max(n, with_ties = FALSE) %>%
      select(-n) %>%
      group_by(!!sym(annotations_col)) %>%
      transmute(Clusters = paste0(!!sym(clusters_col), collapse = ", ")) %>%
      distinct() %>%
      arrange(!!sym(annotations_col))
}


#' Filter AIRR genes from variable features
#'
#' @description
#' Removes IG and/or TR genes from a Seurat object's variable features list.
#' Optionally reports how many remaining features are GEX-only when BCR features are present.
#'
#' @details
#' This will usually just be used as part of [seurat_pipeline()].
#'
#' @param seurat_obj A Seurat object.
#' @param filter_genes Category of genes to remove (e.g. "IG" and/or "TR").
#' @param ensembl_version Ensembl version for gene annotations (e.g. "v114").
#'   If NULL, uses the default in `get_airr_genes()`.
#' @param bcr_features Optional matrix of BCR features (rows = features). If provided, the log message also reports the number of GEX-only features.
#' @param cache_file Passed to [get_airr_genes()]. Path to a cached RDS result to use instead of querying Ensembl.
#'
#' @returns The Seurat object with filtered variable features and Ensembl version saved to `Misc(seurat_obj, "ensembl_version")`.
#' @export
filter_variable_features <- function(seurat_obj, filter_genes,
                                     ensembl_version = NULL,
                                     bcr_features = NULL, cache_file = NULL) {
   airr_genes <- get_airr_genes(category = filter_genes,
                                ensembl_version = ensembl_version,
                                cache_file = cache_file)

   remove_feats <- VariableFeatures(seurat_obj) %in% airr_genes$remove_genes
   VariableFeatures(seurat_obj) <- VariableFeatures(seurat_obj)[!remove_feats]
   Misc(seurat_obj, slot = "ensembl_version") <- airr_genes$ensembl_version

   n_remaining <- length(VariableFeatures(seurat_obj))
   if (!is.null(bcr_features)) {
      n_gex <- length(setdiff(VariableFeatures(seurat_obj),
                              str_replace_all(rownames(bcr_features), "_", "-")))
      cli::cli_inform(c("i" = "After removing {str_c(filter_genes, collapse = '/')} genes: \\
{n_remaining} variable features ({n_gex} GEX-only) remain."))
   } else {
      cli::cli_inform(c("i" = "After removing {str_c(filter_genes, collapse = '/')} genes: \\
{n_remaining} variable features remain."))
   }

   seurat_obj
}


#' Find the right clustering resolution to obtain the desired number of clusters
#'
#' @description
#' This function iteratively tests clustering resolutions in Seurat to find the resolution
#' that yields the desired number of clusters. It uses the specified graph and returns the
#' Seurat object with clusters if successful, or stops if the desired number is exceeded or not found.
#'
#' @details
#' This would typically be used after [seurat_pipeline()].
#'
#' @param seurat_obj The Seurat object.
#' @param graph_name The name of the graph to use for clustering.
#' @param desired_k The desired number of clusters.
#'
#' @returns The Seurat object with clusters at the resolution that matches desired_k.
#' @export
find_k_clusters <- function(seurat_obj, graph_name = "RNA_snn", desired_k) {
   # TODO: increase this range
   for (res in seq(0.1, 2, by = 0.1)) {
      cli::cli_inform("Checking resolution {res}...")

      seurat_obj <-
         suppressWarnings(FindClusters(seurat_obj, resolution = res,
                                       graph.name = graph_name, algorithm = 1, # 4
                                       verbose = FALSE))
      n_clusters <- n_distinct(seurat_obj$seurat_clusters)
      cli::cli_inform("Resolution {res}: {n_clusters} cluster{?s}")

      if (n_clusters == desired_k) {
         # message(paste("Resolution", res, "gives", desired_k, "clusters"))
         return(seurat_obj)
      } else if (n_clusters > desired_k) {
         cli::cli_abort("The number of desired clusters has been exceeded.")
      } else {
         # don't keep the other resolutions
         seurat_obj[[paste0(graph_name, "_res.", res)]] <- c()
      }
   }

   cli::cli_abort("Could not find resolution to match desired clusters.")
}


#' Get IG and TR genes from Ensembl using biomaRt
#'
#' @description
#' This function retrieves IG and TR genes from Ensembl using the biomaRt package. This allows for accurate IG and TR gene names instead of using a search for genes that begin with "IG" or "TR".
#'
#' @details
#' This will usually just be used as part of [seurat_pipeline()].
#'
#' @param genome The genome to use for gene annotation (e.g. "hsapiens" or "mmusculus").
#' @param ensembl_version The Ensembl version to use for gene annotation (e.g. "114").
#' @param category The category of genes to retrieve: "IG" for immunoglobulin genes, "TR" for T cell receptor genes, or both.
#' @param cache_file Optional path to an RDS file. If the file exists, its contents are returned directly without
#'   querying Ensembl. After a successful Ensembl query the result is saved to this path for future offline use.
#'
#' @returns A character vector of IG and/or TR gene names to be filtered out from the most variable features.
#' @export
get_airr_genes <- function(genome = "hsapiens", ensembl_version = NULL,
                           category = c("IG", "TR"), cache_file = NULL) {
   if (!is.null(cache_file) && file.exists(cache_file)) {
      cli::cli_inform(c("i" = "Cache file found: {cache_file}. Loading AIRR gene list from cache."))
      return(readRDS(cache_file))
   }

   retry_ensembl <- function(fn, attempts = 3) {
      for (i in seq_len(attempts)) {
         result <- tryCatch(fn(), error = function(e) {
            if (i < attempts) {
               cli::cli_warn("Ensembl query failed (attempt {i}/{attempts}), retrying...")
               NULL
            } else {
               cli::cli_abort("Ensembl query failed after {attempts} attempts: {conditionMessage(e)}")
            }
         })
         if (!is.null(result)) return(result)
      }
   }

   # if the Ensembl version isn't specified, use the current version
   current_version <-
      retry_ensembl(function() biomaRt::listEnsemblArchives()) %>%
      dplyr::filter(current_release == "*") %>%
      dplyr::pull(name)
   # just the number
   current_version <- stringr::str_split(current_version, " ")[[1]][2]
   if (is.null(ensembl_version)) {
      ensembl_version <- current_version
      cli::cli_inform(c("i" = "Using the current Ensembl version: v{ensembl_version}."))
   }

   # `version` routes to a versioned archive (e.g. e114.ensembl.org) instead of
   # www.ensembl.org, so the two can have different availability
   ensembl <- retry_ensembl(function() biomaRt::useEnsembl(biomart = "genes"))
   dataset <- biomaRt::searchDatasets(mart = ensembl, pattern = genome)$dataset
   ensembl <- tryCatch(
      retry_ensembl(function() biomaRt::useEnsembl(biomart = "genes", dataset = dataset,
                                                    version = ensembl_version)),
      error = function(e) {
         cli::cli_inform(c(
            "Ensembl archive v{ensembl_version} is unreachable; falling back to the current live release (v{current_version}).",
            "i" = "Use {.arg cache_file} to avoid this in future runs."
         ))
         ensembl_version <<- current_version
         retry_ensembl(function() biomaRt::useEnsembl(biomart = "genes", dataset = dataset))
      }
   )

   # filter for the features of interest
   filters <- biomaRt::listFilters(ensembl)
   attributes <- biomaRt::listAttributes(ensembl)
   attribute_names <- c("ensembl_gene_id", "external_gene_name",
                        "gene_biotype", "hgnc_symbol", "description")
   attributes <- dplyr::filter(attributes,
                               name %in% attribute_names &
                                  page == "feature_page")

   # known biotypes
   biotypes <- list(IG = c("IG_V_gene", "IG_V_pseudogene", "IG_LV_gene",
                           "IG_D_gene", "IG_D_pseudogene", "IG_J_gene",
                           "IG_C_gene", "IG_C_pseudogene", "IG_pseudogene"),
                    TR = c("TR_V_gene", "TR_V_pseudogene", "TR_D_gene",
                           "TR_J_gene", "TR_J_pseudogene", "TR_C_gene"))
   biotypes <- as.character(unlist(biotypes[category]))
   biotype_pattern <- str_c("^", category, "_", collapse = "|")

   features_meta <- retry_ensembl(function() biomaRt::getBM(attributes = attribute_names,
                                                             filters = "biotype", values = biotypes,
                                                             mart = ensembl))

   # "biotypes_excl is just unique(features_meta$gene_biotype), which is already exactly biotypes from the BioMart filter. The subsequent filter(gene_biotype %in% biotypes_excl) is a no-op. The commented-out line below suggests the intent was to filter down further — worth revisiting what that logic should be."

   # define the genes to be excluded
   biotypes_excl <- unique(features_meta$gene_biotype)
   # biotypes_excl <- biotypes_excl[grepl(biotype_pattern, biotypes_excl)]
   remove_genes <- features_meta %>%
      dplyr::filter(gene_biotype %in% biotypes_excl) %>%
      dplyr::pull(external_gene_name)
   remove_genes <- remove_genes[remove_genes != ""] # remove empty strings
   remove_genes <- unique(remove_genes)

   result <- list(ensembl_version = ensembl_version, remove_genes = remove_genes)
   cli::cli_inform(c("v" = "Retrieved {length(remove_genes)} AIRR genes from Ensembl v{ensembl_version}."))

   if (!is.null(cache_file)) {
      saveRDS(result, cache_file)
      cli::cli_inform(c("v" = "Saved AIRR gene list to cache file: {cache_file}"))
   } else {
      cli::cli_inform(c("i" = "No cache file specified. To avoid querying Ensembl in future runs, provide a path to save the results using {.arg cache_file}."))
   }

   return(result)
}
