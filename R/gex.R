#' Get IG and TR genes from Ensembl using biomaRt.
#'
#' @description
#' This function retrieves IG and TR genes from Ensembl using the biomaRt package. This allows for accurate IG and TR gene names instead of using a search for genes that begin with "IG" or "TR".
#'
#' @param genome The genome to use for gene annotation (e.g. "hsapiens" or "mmusculus").
#' @param ensembl_version The Ensembl version to use for gene annotation (e.g. "114").
#' @param category The category of genes to retrieve: "IG" for immunoglobulin genes, "TR" for T cell receptor genes, or both.
#'
#' @returns A character vector of IG and/or TR gene names to be filtered out from the most variable features.
#' @export
get_airr_genes <- function(genome = "hsapiens", ensembl_version = NULL,
                           category = c("IG", "TR")) {
   # if the Ensembl version isn't specified, use the current version
   if (is.null(ensembl_version)) {
      ensembl_version <-
         biomaRt::listEnsemblArchives() %>%
         dplyr::filter(current_release == "*") %>%
         dplyr::pull(name)

      # just the number
      ensembl_version <- stringr::str_split(ensembl_version, " ")[[1]][2]
   }

   # search Ensembl's datasets
   ensembl <- biomaRt::useEnsembl(biomart = "ensembl")

   # use listEnsembl() for biomart options (we'll always use genes)
   # you could set `mirror` but you'd have to remove the version arg
   dataset <- biomaRt::searchDatasets(mart = ensembl, pattern = genome)$dataset
   ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset,
                                  version = ensembl_version)
   # could also just do:
   # ensembl <- useMart(biomart = "ensembl", dataset = dataset)

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

   features_meta <- biomaRt::getBM(attributes = attribute_names,
                                   filters = "biotype", values = biotypes,
                                   mart = ensembl)

   # define the genes to be excluded
   biotypes_excl <- unique(features_meta$gene_biotype)
   # biotypes_excl <- biotypes_excl[grepl(biotype_pattern, biotypes_excl)]
   remove_genes <- features_meta %>%
                     dplyr::filter(gene_biotype %in% biotypes_excl) %>%
                     dplyr::pull(external_gene_name)
   remove_genes <- remove_genes[remove_genes != ""] # remove empty strings
   remove_genes <- unique(remove_genes)

   # in case you don't specify the version
   return(list(ensembl_version = ensembl_version, remove_genes = remove_genes))
}


#' Run Seurat's standard pipeline
#'
#' @description
#' This function takes a Seurat object containing combined GEX data and runs the standard Seurat pipeline for normalization, scaling, dimensionality reduction, clustering, and UMAP visualization.
#' It includes optional filtering of IG/TR genes and allows for customization of various parameters such as the number of variable features, principal components, and clustering resolution.
#'
#' @details
#' It is highly recommended to save the resulting object as an RDS or qs file.
#' This pipeline is loosely based on [Seurat's pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial).
#' Unlike previous analyses, `features = rownames(obj)` was removed from the `ScaleData` step since the data is too large and only the top variable features are needed to do `RunPCA`.
#'
#' @param seurat_obj The Seurat object containing combined GEX data.
#' @param nfeatures_RNA Minimum number of RNA features.
#' @param perc_mt Maximum percentage of mitochondrial genes to retain.
#' @param num_features Desired number of variable features.
#' @param num_pcs Number of principal components to compute.
#' @param num_dims Number of dimensions to use for neighbor finding and UMAP.
#' @param cluster_res Clustering resolution parameter.
#' @param filter_genes If specified, filter out genes from this category (e.g. "IG" and/or "TR")
#' @param verbose Print out Seurat's progress messages.
#'
#' @returns A processed Seurat object with normalization, scaling, PCA, clustering, and UMAP.
#' @export
seurat_pipeline <- function(seurat_obj, nfeatures_RNA = 200, perc_mt = 15,
                            num_features = 2000, num_pcs = 30, num_dims = 20,
                            cluster_res = 0.4, filter_genes,
                            ensembl_version = NULL, verbose = TRUE) {
   # filtration
   if ("percent.mt" %in% names(seurat_obj[[]])) {
      seurat_obj <-
         subset(seurat_obj,
                subset = nFeature_RNA > nfeatures_RNA & percent.mt < perc_mt)
   } else {
      warning("No filtration was performed upon this object.")
   }

   # standard normalization
   seurat_obj <- NormalizeData(seurat_obj,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000, verbose = verbose)
   if ("ADT" %in% names(seurat_obj@assays)) {
      # normalize across the cells, not the features
      seurat_obj <- NormalizeData(seurat_obj,
                                  normalization.method = "CLR", margin = 2,
                                  assay = "ADT", verbose = verbose)
   }

   # highly variable features
   seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                      nfeatures = num_features,
                                      verbose = verbose)

   # scaling
   # note: `features = rownames(seurat_obj)` can cause crashes
   seurat_obj <- ScaleData(seurat_obj, verbose = verbose)

   # filter out the IG and/or TR genes
   if (!rlang::is_missing(filter_genes)) {
      # TODO: check the object for the species
      airr_genes <- get_airr_genes(category = filter_genes,
                                   ensembl_version = ensembl_version)

      remove_feats <- VariableFeatures(seurat_obj) %in% airr_genes$remove_genes
      VariableFeatures(seurat_obj) <- VariableFeatures(seurat_obj)[!remove_feats]

      cat(paste("After removing", str_c(filter_genes, collapse = "/"),
                "genes, the total number of variable features is:",
                length(VariableFeatures(seurat_obj)), "\n"))
   }

   # TODO: pull from get_airr_genes, list latest if NULL
   # save the Ensembl version
   Misc(seurat_obj, slot = "ensembl_version") <- airr_genes$ensembl_version

   # TODO: call this rpca instead
   # dimensionality reduction
   seurat_obj <- RunPCA(seurat_obj, npcs = num_pcs, verbose = verbose)

   # neighbor detection
   seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:num_dims,
                               verbose = verbose)

   # might not always want to perform clustering, so make it optional
   if (!rlang::is_missing(cluster_res)) {
      seurat_obj <- FindClusters(seurat_obj, resolution = cluster_res,
                                 verbose = verbose)

      # fix the cluster levels (for some reason they sort alphabetically now)
      for (res in cluster_res) {
         res <- paste0("RNA_snn_res.", res)
         seurat_obj[[res]] <-
            factor(seurat_obj[[]][[res]],
                   str_sort(unique(seurat_obj[[]][[res]]), numeric = TRUE))
      }
   }

   # TODO: call this rna.umap instead
   # for visualization
   seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:num_dims,
                         verbose = verbose)

   return(seurat_obj)
}


#' Visualize how many cells are in each Seurat cluster or cell type
#'
#' @description
#' This function generates a bar plot showing the number of cells in each Seurat cluster or cell type, optionally filled by another annotation (e.g. annotated clusters).
#' It allows for customization of the x-axis, fill colors, and plot titles.
#'
#' @param seurat_obj The Seurat object containing the clusters and/or cell type annotations to plot.
#' @param tissue_type The tissue type of interest e.g. "Blood" or "Skin".
#' @param clrs_specific The specific color palette (should be named).
#' @param clusters_col The column to plot on the x axis (e.g. seurat_clusters).
#' @param fill_col The column to fill by (e.g. annotated_clusters).
#' @param fill_col_name The label for the fill aesthetic.
#' @param x_axis What to plot on the x axis: "Cluster" or "Cell Type".
#' @param details The optional subtitle.
#'
#' @returns A ggplot bar plot of cell counts.
#' @export
plot_counts_cluster <- function(seurat_obj, tissue_type = "", clrs_specific,
                                clusters_col, fill_col, fill_col_name,
                                x_axis = "Cluster", # add_zeroes = FALSE,
                                details = NULL) {
   # if you want to use the default Seurat colors
   if (rlang::is_missing(clrs_specific)) {
      clrs_specific <- hue_pal()(n_distinct(seurat_obj[[clusters_col]]))
   }

   if (rlang::is_missing(fill_col)) {
      fill_col <- NULL
      fill_by <- clusters_col
      fill_col_name <- clusters_col

      plot_title <- paste(tissue_type, x_axis)
   } else {
      fill_by <- fill_col

      # they shouldn't both be missing
      if (rlang::is_missing(fill_col_name)) fill_col_name <- fill_col

      plot_title <- paste(tissue_type, fill_col_name, # paste0(fill_col_name, "s"),
                          "per", x_axis)
   }

   data <- data.frame(table(seurat_obj[[c(clusters_col, fill_col)]])) %>%
      rename(Count = Freq) %>%
      filter(Count != 0) # we don't want the zeroes

   # if (add_zeroes) {
   #   data <- bind_rows(data, data.frame(clusters = add_zeroes,
   #                                      Count = rep(0, length(add_zeroes)))) %>%
   #             arrange(clusters)
   # }

   p <- ggplot(data,
               aes(x = !!sym(clusters_col), y = Count, fill = !!sym(fill_by))) +
      geom_col(color = "black", linewidth = 0.2) +
      labs(title = plot_title, subtitle = details, x = x_axis,
           fill = fill_col_name) +
      scale_fill_manual(values = clrs_specific) +
      theme_bw + labels_standard +
      guides(fill = guide_legend(ncol = 1)) # one column legend

   # doesn't work with multiple fills
   # could just plot the sum on top
   if (is.null(fill_col)) {
      p <- p +
         geom_text(aes(label = Count), vjust = -1, size = 3) +
         theme(legend.position = "none")
   } else {
      p <- p +
         geom_text(aes(label = after_stat(y), group = clusters_col),
                   stat = "summary", fun = sum, vjust = -1)
   }

   if (x_axis == "Cell Type") p <- p + labels_rotate_x

   # return the plot
   return(p)
}


#' Get specific markers from a marker genes database
#'
#' @description
#' This function returns specific markers from the "all" marker genes dataframe
#' based on various filtering criteria including source, cell types, and tissue types.
#'
#' @param markers_df The database of marker genes
#' @param sources Optional vector of sources - who the markers came from if you want specific origins.
#' @param contains Optional string to catch multiple cell types (e.g. "mDCs" and "pDCs").
#' @param tissue_types Optional vector of tissue types (e.g. "blood", "skin", etc.).
#' @param cell_types Optional vector of cell types you want markers for.
#' @param alphabetize_types Whether to return the markers alphabetized for each cell type.
#' @param alphabetize_all Whether to return all of the selected markers alphabetized.
#'
#' @returns A character vector of unique gene features/markers.
#' @export
get_features_from_all <- function(markers_df, sources, contains, tissue_types,
                                  cell_types, alphabetize_types = TRUE,
                                  alphabetize_all = TRUE) {
   # make specific selections if needed
   if (!rlang::is_missing(sources)) {
      markers_df <- filter(markers_df, Source %in% sources)
   }
   if (rlang::is_missing(cell_types)) {
      cell_types <- unique(markers_df$Cell_Type)
   }
   if (!rlang::is_missing(contains)) {
      cell_types <- grep(paste(contains, collapse = "|"),
                         cell_types, value = TRUE)
   }
   if (!rlang::is_missing(tissue_types)) {
      markers_df <- filter(markers_df, Tissue_Type %in% tissue_types)
   }

   features <- c()
   for (cell_type in cell_types) {
      markers <- (filter(markers_df, Cell_Type == cell_type))$Marker
      if (alphabetize_types) markers <- sort(markers)
      features <- append(features, markers)
   }

   if (alphabetize_all) features <- str_sort(features, numeric = TRUE)

   return(unique(features)) # DotPlot doesn't work with duplicated features
}


#' Get the information needed for using `add_info_bar()` on a `DotPlot`
#'
#' @description
#' This function takes the filtered markers dataframe and selects the relevant columns to be used for adding an info bar to a DotPlot.
#' It renames the Marker column to features.plot and formats the Cell_Type_Full column for better readability.
#'
#' @param markers_df The markers data.frame filtered to match your input features
#'
#' @returns A data.frame with Cell_Type_Full and features.plot columns
#' @export
get_cell_types <- function(markers_df) {
   markers_df %>%
      select(Cell_Type_Full, Marker) %>%
      rename(features.plot = Marker) %>%
      mutate(Cell_Type_Full = str_replace_all(Cell_Type_Full, "_", " "))
}


#' Generate a title for a `DotPlot`
#'
#' @description
#' This function generates a title for a DotPlot by combining a provided dataset description with a formatted list of marker sources.
#' The marker sources are sorted, formatted to replace underscores with spaces, and concatenated into a comma-separated string enclosed in parentheses.
#'
#' @param plot_title Dataset description
#' @param marker_sources The list of marker sources
#'
#' @returns A string with sources comma-separated and in parentheses
#' @export
gen_dot_title <- function(plot_title = "", marker_sources) {
   paste0(plot_title, " (",
          str_c(sort(str_replace_all(marker_sources, "_", " ")),
                collapse = ", "),
          ")")
}


#' Display markers from a filtered marker genes database as a table
#'
#' @description
#' This function takes a filtered marker genes dataframe and returns a formatted
#' table showing markers organized by cell type.
#'
#' @param filtered_markers_df A filtered dataframe from the marker genes database
#'   containing at least Cell_Type and Marker columns.
#'
#' @returns A formatted table showing markers grouped by cell type.
#' @export
source_markers <- function(filtered_markers_df) {
   filtered_markers_df %>%
      select(Cell_Type, Marker) %>%
      group_by(Cell_Type) %>%
      distinct() %>%
      summarize(Markers = paste(Marker, collapse = ", "))
}


#' Add a user-specified list of cluster annotations to a Seurat object
#'
#' @description
#' This function adds a user-specified list of cluster annotations to the Seurat object.
#' It can be used for both manual and automated annotation, updating cell identities
#' and adding annotation columns to the metadata.
#'
#' @param seurat_obj The Seurat object to annotate.
#' @param annotations_df Data frame containing cluster-to-cell-type mappings, typically
#'   with "Cluster" and "CellType" columns.
#' @param cell_types_col The name of the column containing the cell type annotations.
#' @param relabel Whether to update the Seurat object's active identities. Sometimes
#'   you just want to add the metadata.
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
   # prepare the annotation information
   annotations <- annotations_df[[cell_types_col]] # you only need the cell type information
   names(annotations) <- levels(seurat_obj[[clusters_col]] %>% pull())

   # relabel the Seurat clusters
   current_idents <- Idents(seurat_obj)
   Idents(seurat_obj) <- clusters_col # reset to original clusters
   seurat_obj <- RenameIdents(seurat_obj, annotations)

   # alphabetize the cell types
   if (alphabetize) {
      Idents(seurat_obj) <- factor(Idents(seurat_obj),
                                   levels = sort(levels(seurat_obj)))
   }

   # useful metadata (e.g. if you want to have multiple annotation sets)
   seurat_obj[[annotations_col]] <- Idents(seurat_obj)

   if (relocate) {
      seurat_obj@meta.data <-
         seurat_obj[[]] %>%
         relocate(!!sym(annotations_col), .after = !!sym(clusters_col))
   }

   # if you just wanted to add the metadata
   if (!relabel) Idents(seurat_obj) <- current_idents

   return(seurat_obj)
}


#' Run automated cell type annotation
#'
#' @description
#' This function runs automated cell type annotation using either `Azimuth` or `CellTypist`.
#'
#' @details
#' Supports Azimuth and CellTypist annotation methods.
#' Assumes that the `Cells()` of `seurat_obj` are properly formatted (i.e. unique).
#' For `Azimuth` with a Seurat v5 object, all of the layers have to be joined.
#' For `CellTypist`, assumes the H5AD file and predictions have already been generated.
#'
#' @param seurat_obj The Seurat object. Must be the path to a H5AD object if using CellTypist.
#' @param annotation_method Which method to use: c("Azimuth", "CellTypist")
#' @param reference Reference or model to use for prediction. Defaults to "pbmcref" (for Azimuth).
#' @param azimuth_assay Assay to use for Azimuth
#' @param azimuth_levels Levels to process for Azimuth e.g. c("l1", "l2", "l3")
#'
#' @returns A data.frame with the annotations for each cell
#' @export
automated_annotation <- function(seurat_obj, annotation_method,
                                 reference = "pbmcref", azimuth_assay = "RNA",
                                 azimuth_levels = c("l1", "l2", "l3")) {
   # validate input
   valid_methods <- c("Azimuth", "CellTypist")
   if (!any(annotation_method %in% valid_methods)) {
      stop("Method must be one of: ", paste(valid_methods, collapse = ", "))
   }

   # run Azimuth annotation
   if (annotation_method == "Azimuth") {
      if (!requireNamespace("Azimuth", quietly = TRUE)) {
         stop("The Azimuth package is required for annotation.")
      }

      cat("Running Azimuth annotation...\n")
      # could list how many of the specified features are not present in the reference
      obj_azimuth <- Azimuth::RunAzimuth(query = seurat_obj,
                                         reference = reference,
                                         assay = azimuth_assay)

      # save annotations by cell
      annotations <- data.frame(cell_id = Cells(obj_azimuth),
                                mapping.score = obj_azimuth[[]]$mapping.score)

      # process each level of annotation
      for (level in azimuth_levels) {
         full_level <- paste0("predicted.celltype.", level)

         if (!full_level %in% colnames(obj_azimuth[[]])) {
            warning(paste("Column", full_level, "not found. Skipping level", level))
            next
         }

         # add in the predicted cell types and the confidence score
         annotations <- bind_cols(annotations,
                                  obj_azimuth[[]] %>%
                                     remove_rownames() %>%
                                     select(all_of(full_level),
                                            paste0(full_level, ".score")))

         cat(paste("Generated Azimuth", level, "annotations\n"))
      }
   }

   # run CellTypist annotation
   if (annotation_method == "CellTypist") {
      cat("Processing CellTypist annotation...\n")

      if (!requireNamespace("reticulate", quietly = TRUE)) {
         stop("reticulate package is required for running CellTypist predictions")
      }

      # load model and data
      model <- celltypist$models$Model$load(model = reference)
      obj_h5ad <- scanpy$read_h5ad(filename = seurat_obj)

      # run predictions (on a per cell level)
      predictions <- celltypist$annotate(filename = obj_h5ad, model = model,
                                         majority_voting = FALSE)

      # get the predicted cell types and the confidence score
      annotations <- predictions$to_adata()
      annotations <- annotations$obs %>%
         select(cell_id, predicted_labels, conf_score) %>%
         remove_rownames()
   }

   return(annotations)
}


#' Map cell types to Seurat clusters
#'
#' @description
#' This function maps cell types to Seurat clusters by counting the number of cells in each cluster that belong to each cell type, and then assigning the most common cell type to each cluster.
#' It returns a data frame with the assigned cell types for each cluster.
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


#' Find the right clustering resolution to obtain the desired number of clusters
#'
#' @description
#' This function iteratively tests clustering resolutions in Seurat to find the resolution
#' that yields the desired number of clusters. It uses the specified graph and returns the
#' Seurat object with clusters if successful, or stops if the desired number is exceeded or not found.
#'
#' @param seurat_obj The Seurat object.
#' @param graph_name The name of the graph to use for clustering.
#' @param desired_k The desired number of clusters.
#'
#' @returns The Seurat object with clusters at the resolution that matches desired_k.
#' @export
find_k_clusters <- function(seurat_obj, graph_name = "RNA_snn", desired_k) {
   for (res in seq(0.1, 2, by = 0.1)) {
      cat(paste0("Checking resolution ", res, ": "))

      seurat_obj <-
         suppressWarnings(FindClusters(seurat_obj, resolution = res,
                                       graph.name = graph_name, algorithm = 1, # 4
                                       verbose = FALSE))
      n_clusters <- n_distinct(seurat_obj$seurat_clusters)
      cat(paste(n_clusters, "clusters\n"))

      if (n_clusters == desired_k) {
         # message(paste("Resolution", res, "gives", desired_k, "clusters"))
         return(seurat_obj)
      } else if (n_clusters > desired_k) {
         stop("The number of desired clusters has been exceeded.")
      } else {
         # don't keep the other resolutions
         seurat_obj[[paste0(graph_name, "_res.", res)]] <- c()
      }
   }
   stop("Could not find resolution to match desired clusters.")
}
