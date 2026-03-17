#' Read in and process output files from nf-core/airrflow
#'
#' @description
#' This function reads airrflow's repertoire files, processes them to add
#' subject and sample information, computes CDR3 amino acid properties, and adds
#' gene family information.
#'
#' @details
#' Only written for BCR data right now.
#' v4.0 has `clone_size_count` and `clone_size_freq` columns.
#' v4.3.1 has `duplicate_count` and `light_only_cell` columns.
#' There can still be some `NA` `c_call`s.
#'
#' @param dataset_path The path to the dataset directory.
#' @param version_airrflow The airrflow version (as a string).
#'
#' @returns A processed AIRR-formatted data.frame with several columns added.
#' @export
process_airrflow <- function(dataset_path, version_airrflow) {
  # deal with possible format e.g. 4.0 or 4.3.1
  version_airrflow_num <-
    str_split_1(version_airrflow, pattern = "")[1:3] %>% str_c(collapse = "")

  # determine where the files are depending on the airrflow version
  path_repertoire <- file.path(dataset_path, "airrflow",
                               "bcr", version_airrflow, "results")
  if (version_airrflow_num < 4.2) { # or 4.3.0?
    path_repertoire <- file.path(path_repertoire,
                                 "repertoire_comparison", "repertoires")
  } else if (version_airrflow_num == "4.3") {
    path_repertoire <- file.path(path_repertoire,
                                 "clonal_analysis", "define_clones",
                                 "all_reps_clone_report", "repertoires")
  } else { # v5
    path_repertoire <- file.path(path_repertoire,
                                 "clonal_analysis", "clonal_assignment",
                                 "all_reps_clone_report", "repertoires")
  }

  # read in the airrflow output files (there's one per subject)
  rep_files <- list.files(path_repertoire, full.names = TRUE)
  if (length(rep_files) == 0) {
    stop("No files found, are you sure you're using the correct directory?")
  }
  combined_bcr <- map(rep_files, read_tsv, show_col_types = FALSE,
                      # sometimes a "sex" column will use M or F which then
                      # incorrectly reads as logical
                      col_types = readr::cols(sex = readr::col_character())) %>%
                  bind_rows()

  # add in useful columns
  combined_bcr <-
    combined_bcr %>%
    # create unique clone IDs (readr reads clone_id as numeric)
    # mutate(clone_id = as.character(clone_id)) %>%
    mutate(clone_id_unique = paste0(subject_id, "_", clone_id), # sample_id
           .after = clone_id)

  # update the cell id column to make it unique (by sample)
  combined_bcr <-
    combined_bcr %>%
    mutate(cell_id_original = cell_id, .after = cell_id) %>%
    separate(cell_id, sep = "-", into = "cell_id", # cell_id_formatted
             remove = FALSE, extra = "drop") %>%
    mutate(cell_id = paste0(sample_id, "_", cell_id))

  # add isotype information based on c_call
  if (!"isotype" %in% colnames(combined_bcr)) {
    combined_bcr <-
      combined_bcr %>%
      mutate(isotype = case_when(str_detect(c_call, "IGHA") ~ "IgA",
                                 str_detect(c_call, "IGHD") ~ "IgD",
                                 str_detect(c_call, "IGHE") ~ "IgE",
                                 str_detect(c_call, "IGHG") ~ "IgG",
                                 str_detect(c_call, "IGHM") ~ "IgM"))
  }

  # filter out heavy chains with unassigned isotypes
  cat(paste("There are",
            nrow(filter(combined_bcr, locus == "IGH", is.na(c_call))),
            "heavy chains with `NA` c_calls.\n"))
  combined_bcr <- combined_bcr %>% filter(!is.na(c_call) | locus != "IGH")

  # filter out cells with unpaired light chains
  # airrflow does this, but we have to do it again since we now filtered out
  # more heavy chains
  unpaired_light_chains <- combined_bcr %>%
    group_by(cell_id) %>%
    filter(!any(locus == "IGH")) %>%
    pull(cell_id)
  cat(paste("There are", length(unpaired_light_chains), "unpaired light",
            "chains (no corresponding heavy chain with the same cell id).\n"))
  combined_bcr <- combined_bcr %>% filter(!cell_id %in% unpaired_light_chains)

  # filter out heavy chains with light chains assigned to the isotype
  # this comes from the "multi" column in the original 10x data
  # see: https://kb.10xgenomics.com/s/article/360001715971-What-is-the-source-of-Multi-chains-in-the-contig-annotations-csv-files
  multi_heavy_chains <- combined_bcr %>%
    filter(locus == "IGH", is.na(isotype)) %>%
    pull(cell_id)
  cat(paste("There are", length(multi_heavy_chains), "heavy chains with light",
            "chain c_calls.\n"))
  combined_bcr <- combined_bcr %>% filter(!cell_id %in% multi_heavy_chains)

  # preserve the original c_call column and create a simplified version
  combined_bcr <-
    combined_bcr %>%
    mutate(c_call_original = c_call, .before = c_call,
           c_call = getFamily(c_call)) # don't return multiple hits

  # add alakazam's amino acid properties for heavy chain sequences
  # excludes non-informative positions
  if (any(combined_bcr$locus == "IGH")) {
    amino_acids <-
      aminoAcidProperties(data = filter(combined_bcr, locus == "IGH"),
                          seq = "junction", nt = TRUE, trim = TRUE,
                          label = "cdr3") # junction if this is left out
    combined_bcr <-
      left_join(combined_bcr,
                select(amino_acids, sequence_id, starts_with("cdr3_")),
                by = join_by(sequence_id))
  }
  # adds the "cdr3_aa_length", "cdr3_aa_gravy", "cdr3_aa_bulk",
  # "cdr3_aa_polarity", "cdr3_aa_aliphatic", "cdr3_aa_charge", "cdr3_aa_acidic",
  # "cdr3_aa_basic" and "cdr3_aa_aromatic" columns

  # add gene family information
  combined_bcr <- add_family_info(combined_airr = combined_bcr)

  cat(paste("There are", nrow(filter(combined_bcr)),
            "total chains in the combined data."))

  return(combined_bcr)
}


#' Plot UMAPs with AIRR overlays
#'
#' @description
#' This function plots the annotated UMAPs alongside the BCR/TCR overlays.
#'
#' @details
#' lightgray is Seurat's default background color and cells_total has to be a list for the labels to work later.
#' This could probably be replaced with `plot_umap()`.
#'
#' @param seurat_obj The Seurat object (must contain Has_BCR/Has_TCR cols)
#' @param tissue_type The tissue type of interest.
#' @param airr_type BCR or TCR
#' @param clrs_specific The specific color palette (should be named).
#' @param barcode_col The barcode column name: barcode, cell_id, or Cell_ID_Unique.
#' @param plot_label Whether or not to include the labels.
#' @param plot_by The grouping method: all samples, by dataset, by sample, by isotype.
#' @param ncol_sample The number of columns for sample-wise plots.
#'
#' @returns A Seurat UMAP plot with AIRR cells highlighted by the specified grouping.
#' @export
plot_immune_overlay <- function(seurat_obj, tissue_type, airr_type,
                                clrs_specific, plot_by = "all",
                                barcode_col = "cell_id",
                                plot_label = FALSE, ncol_sample = 6) {
  # plot options
  pt_size <- 0.2
  label_size <- 4
  sizes_highlight <- 0.2 # should probably be the same size as the points
  plot_title <- paste(tissue_type, airr_type, "Overlay")

  # if you want to use the default Seurat colors
  # if (rlang::is_missing(clrs_specific)) clrs_specific <- hue_pal()(num_clusts)

  # AIRR info
  airr_col <- paste0("Has_", airr_type)
  combined_vdj_gex <- filter(seurat_obj[[]], !!rlang::sym(airr_col))

  if (plot_by == "dataset" || plot_by == "sample") {
    # set up the cells to be highlighted
    cells_total <- c()
    datasets <- as.character(unique(combined_vdj_gex$Dataset))

    for (dataset in datasets) {
      cells_dataset <- filter(combined_vdj_gex,
                              Dataset == dataset)[[barcode_col]]
      cells_total[dataset] <- list(cells_dataset)
    }

    # main plot (coloring the samples by dataset is helpful)
    if (plot_by == "sample") {
      p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                    split.by = "sample_id",
                    repel = TRUE,
                    cells.highlight = cells_total,
                    sizes.highlight = sizes_highlight,
                    ncol = ncol_sample, raster = FALSE) +
        labs(title = paste(plot_title, "by Sample"), color = "Dataset")
    } else { # dataset
      p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                    label = plot_label, label.size = label_size,
                    repel = TRUE,
                    cells.highlight = cells_total,
                    sizes.highlight = sizes_highlight,
                    raster = FALSE) +
        labs(title = paste(plot_title, "by Dataset"), color = "Dataset")
    }

    # the legend
    p <- p + scale_color_manual(labels = c(paste0("non-", airr_type), rev(datasets)),
                                values = c("lightgray", unname(clrs_specific[rev(datasets)])))
  } else if (plot_by == "isotype") {
    # select the cells to be highlighted
    # Idents(seurat_obj) <- "isotype"
    # cells_total <- CellsByIdentities(seurat_obj)
    # cells_total <- lapply(cells_total, function(x) x[!is.na(x)])
    # cells_total["NA"] <- c()
    # isotypes <- rev(names(cells_total)) # need to reverse for plotting
    cells_total <- c()
    isotypes <- sort(unique(combined_vdj_gex$isotype))

    for (isotype in isotypes) {
      cells_isotype <- filter(combined_vdj_gex, isotype == isotype)[[barcode_col]]
      cells_total[isotype] <- list(cells_isotype)
    }

    # main plot
    p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                  label = plot_label, label.size = label_size,
                  repel = TRUE,
                  cells.highlight = cells_total,
                  sizes.highlight = sizes_highlight,
                  raster = FALSE) +
      labs(title = paste(plot_title, "by Isotype"), color = "Isotype")

    # the legend
    p <- p + scale_color_manual(labels = c(paste0("non-", airr_type), rev(isotypes)),
                                values = c("lightgray", unname(clrs_specific[rev(isotypes)])))
  } else { # all
    # set up the cells to be highlighted
    cells_total <- combined_vdj_gex[[barcode_col]]

    # plot all samples
    p <- UMAPPlot(object = seurat_obj, pt.size = pt_size,
                  label = plot_label, label.size = label_size,
                  repel = TRUE,
                  cols.highlight = clrs_specific[airr_type],
                  cells.highlight = cells_total,
                  sizes.highlight = sizes_highlight,
                  raster = FALSE) +
      labs(title = plot_title, color = "Data Type")

    # the legend
    p <- p + scale_color_manual(name = "Data Type",
                                labels = c(paste0("non-", airr_type), airr_type),
                                values = c("lightgray", unname(clrs_specific[airr_type])))
  }

  # standardize the labels
  p <- p & labels_standard & clean_umap

  return(p)
}


#' Add in family and gene information from `alakazam`
#'
#' @description
#' Extracts and adds V, D, and J gene family and gene information to an AIRR-formatted
#' dataframe using the alakazam package functions.
#'
#' @details
#' biomaRt also has a getGene function, so we have to be specific.
#'
#' @param combined_airr An AIRR-formatted data.frame.
#'
#' @returns A data.frame with six new columns containing gene family and gene information.
#' @export
add_family_info <- function(combined_airr) {
  if ("v_call" %in% names(combined_airr)) {
    combined_airr <-
      combined_airr %>%
      mutate(v_call_family = alakazam::getFamily(combined_airr$v_call),
             v_call_gene = alakazam::getGene(combined_airr$v_call))
  }

  if ("d_call" %in% names(combined_airr)) {
    combined_airr <-
      combined_airr %>%
      mutate(d_call_family = alakazam::getFamily(combined_airr$d_call),
             d_call_gene = alakazam::getGene(combined_airr$d_call))
  }

  if ("j_call" %in% names(combined_airr)) {
    combined_airr <-
      combined_airr %>%
      mutate(j_call_family = alakazam::getFamily(combined_airr$j_call),
             j_call_gene = alakazam::getGene(combined_airr$j_call))
  }

  return(combined_airr)
}


#' Convert family and gene information to sorted factors
#'
#' @description
#' Converts the gene family and gene columns added by `add_family_info()` to
#' properly ordered factors using numeric sorting.
#'
#' @details
#' For after `add_family_info()` has been run
#'
#' @param combined_airr An AIRR-formatted data.frame. This function ensures that gene
#' families and genes are ordered correctly (e.g., IGHV1, IGHV2, IGHV10 instead
#' of IGHV1, IGHV10, IGHV2).
#'
#' @returns A data.frame with up to six columns converted to sorted factors.
#' @export
factor_family_info <- function(combined_airr) {
  if ("v_call" %in% names(combined_airr)) {
    combined_airr <-
      combined_airr %>%
      mutate(v_call_family =
               factor(v_call_family,
                      str_sort((unique(v_call_family)), numeric = TRUE)),
             v_call_gene =
               factor(v_call_gene,
                      str_sort((unique(v_call_gene)), numeric = TRUE)),)
  }

  if ("d_call" %in% names(combined_airr)) {
    combined_airr <-
      combined_airr %>%
      mutate(d_call_family =
               factor(d_call_family,
                      str_sort((unique(d_call_family)), numeric = TRUE)),
             d_call_gene =
               factor(d_call_gene,
                      str_sort((unique(d_call_gene)), numeric = TRUE)),)
  }

  if ("j_call" %in% names(combined_airr)) {
    combined_airr <-
      combined_airr %>%
      mutate(j_call_family =
               factor(j_call_family,
                      str_sort((unique(j_call_family)), numeric = TRUE)),
             j_call_gene =
               factor(j_call_gene,
                      str_sort((unique(j_call_gene)), numeric = TRUE)),)
  }

  return(combined_airr)
}


#' Bin the mutation frequency
#'
#' @description
#' Bins the `mu_freq` column in a Seurat object into specified categories (e.g., 0%, 0-1%, 1-5%, etc.) for easier visualization and analysis.
#' The function creates new columns with binned mutation frequencies based on the provided number of bins.
#'
#' @details
#' The bins are (most likely) not going to be equal sizes.
#'
#' @param seurat_obj The Seurat object.
#' @param num_bins The number of bins to split `mu_freq` into. Must be at least one of 2, 3, or 5.
#'
#' @returns The provided Seurat object with a new binned mu_freq column.
#' @export
bin_mu_freq <- function(seurat_obj, num_bins = c(2, 3, 5)) {
  # get the current mutation frequencies
  mu_freqs <- seurat_obj$mu_freq
  mu_freqs[is.na(mu_freqs)] <- -1
  # mu_freqs <- mu_freqs[!is.na(mu_freqs)]
  max_mu_freq <- max(mu_freqs)

  # small epsilon for upper bound of the cut
  eps <- 1e-3

  if (5 %in% num_bins) {
    # bins <- seq(0, round(max_mu_freq, digits = 1), by = 0.05)
    bins <- c(-1, 0, 0.01, 0.04, 0.1, max_mu_freq + eps) # for rounding
    labels <- c("0%", "0% to 1%", "1% to 5%", "5% to 10%",
                paste0("10% to ", round(100 * max_mu_freq), "%"))

    seurat_obj$mu_freq_bins <- cut(mu_freqs, breaks = bins, labels = labels)
  }

  # use 0.01 just in case we have novel alleles
  if (3 %in% num_bins) {
    bins <- c(-1, 0, 0.01, max_mu_freq + eps)
    labels <- c("0%", "0% to 1%", ">1%")

    seurat_obj$mu_freq_bins_fewer <-
      cut(mu_freqs, breaks = bins, labels = labels)
  }

  if (2 %in% num_bins) {
    bins <- c(-1, 0.01, max_mu_freq + eps)
    labels <- c("0% to 1%", ">1%")

    seurat_obj$mu_freq_bins_binary <-
      cut(mu_freqs, breaks = bins, labels = labels)
  }

  return(seurat_obj)
}


#' Process BCR features for integration
#'
#' @description
#' Processes BCR features by renaming columns, converting ordered/categorical
#' variables to numeric representations, and normalizing all features. This
#' prepares BCR data for integration with gene expression data in a Seurat object.
#'
#' @details Uses `recipes` to:
#'   1. Rename columns to distinguish from existing metadata
#'   2. Convert ordered factors to ordinal scores OR one-hot encode categoricals
#'   3. Center and normalize all numeric predictors
#'   4. Transpose the result for compatibility with Seurat assays
#'  You could do something like janitor::clean_names(bcr_features) to remove the underscores from the column names right off the bat, but one-hot encoding will add names with underscores automatically unless you messing with the `naming` argument in `step_dummy()`.
#'
#' @param bcr_features A data frame containing BCR features (e.g., isotype, mutation
#'   frequency). May contain ordered factors, numeric, or categorical variables.
#'
#' @return A transposed matrix of processed features where:
#'   - Ordered variables are converted to numeric scores and suffixed with "-ordered"
#'   - Numeric variables are suffixed with "-scaled"
#'   - Categorical variables are one-hot encoded (all categories kept)
#'   - All numeric predictors are normalized to mean = 0 and sd = 1
#'   - Underscores are removed from feature names (so Seurat doesn't throw a warning)
#' @export
process_bcr_features <- function(bcr_features) {
  # TODO: print out what steps were taken

  # rename (to distinguish from existing metadata) and convert the columns
  if (any(sapply(bcr_features, is.ordered))) {
    bcr_features <- bcr_features %>%
      rename_with(~paste0(., "-ordered"), where(is.ordered)) %>%
      rename_with(~paste0(., "-scaled"), where(is.numeric))

    # convert ordered variables to numeric
    ref_cell <- recipe( ~ ., data = bcr_features) %>%
      step_ordinalscore(all_ordered_predictors())
  } else {
    bcr_features <- bcr_features %>%
      rename_with(~paste0(., "-scaled"), everything())

    # one-hot encoding of categorical variables
    ref_cell <- recipe( ~ ., data = bcr_features) %>%
      # cover missing values just in case
      step_unknown(all_nominal_predictors()) %>%
      step_dummy(all_nominal_predictors(), one_hot = TRUE)
  }

  # TODO: add step_zv()?

  # centers and normalizes to mean of 0 and sd of 1
  ref_cell <- ref_cell %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = bcr_features)
  bcr_features <- bake(ref_cell, new_data = NULL) %>% t()

  # "Feature names cannot have underscores ('_'), replacing with dashes ('-')"
  rownames(bcr_features) <- gsub("_", "-", rownames(bcr_features))

  return(bcr_features)
}


#' Convert the output of an embedding method to a matrix
#'
#' @description
#' This function takes the output from an embedding method (e.g., AntiBERTy, ESM2, immune2vec) and converts it into a `Matrix` matrix format that can be used for downstream analysis.
#' It handles the mapping of cell IDs from the embeddings to the combined AIRR data, ensuring that only cells present in both datasets are retained.
#' The resulting matrix has features as rows and cells as columns, with appropriate column names based on the number of dimensions in the embeddings.
#'
#' @details
#' Assume that all inputs can be provided as tsv files.
#' The AMULETy outputs always have a column named "cell_id".
#' immune2vec does not include the cell ids in the output.
#' Uses `Matrix()` instead of `as.matrix()` since the latter only returns dense matrices.
#' AntiBERTy runs on 512 dimensions, AntiBERTa2 and BALM-paired run on 1024 dimensions, and ESM2 runs on 1280 dimensions (through AMULETy).
#' immune2vec is usually run with 100 dimensions.
#'
#' @param embeddings The data.frame embeddings output.
#' @param combined_airr The combined output from airrflow/Immcantation. Must contain columns called "cell_id" and "cell_id_original".
#' @param combined_airr_input The data.frame provided to immune2vec; contains translated sequences.
#'
#' @returns A `Matrix` of embeddings
#' @export
convert_embeddings <- function(embeddings, combined_airr, combined_airr_input) {
  # check just in case
  if (!all(c("cell_id", "cell_id_original") %in% colnames(combined_airr))) {
    stop("One or both of the required columns (`cell_id`, `cell_id_original`) ",
         "are missing from combined_airr.")
  }

  # immune2vec only
  if (!"cell_id" %in% colnames(embeddings)) {
    # just in case
    if (rlang::is_missing(combined_airr_input)) {
      stop("Please provide the original file used as input for running immune2vec.")
    }

    # we're only interested in the heavy chains
    combined_airr_input <- filter(combined_airr_input, locus == "IGH")

    # it is possible that there are more translated sequences than the combined
    # AIRR data if it was processed through airrflow with AMULETy (since the
    # former splits off the data to be translated before all of the filters are
    # done)
    embeddings$cell_id <- combined_airr_input$cell_id
  }

  # check for cells in the embeddings that are not found in combined_airr
  # which could be due to filters on the latter (e.g. removing NA c_calls)
  embeddings_only_cell_ids <-
    setdiff(embeddings$cell_id, combined_airr$cell_id_original)
  if (length(embeddings_only_cell_ids) > 0) { # or any(is.na(cell_ids))
    embeddings <- filter(embeddings, !cell_id %in% embeddings_only_cell_ids)
  }
  cat(paste(length(embeddings_only_cell_ids), "cells in the embedding with no",
            "corresponding cell ids in the combined AIRR file were removed."))

  # depending on how AMULETy was run and read in, there might be sample names
  # in their own column
  if ("sample_id" %in% colnames(embeddings)) {
    # pull the proper cell names (with sample names appended)
    # use the sample_id column since there is the tiniest chance that barcodes
    # can collide across samples
    cell_ids <- embeddings %>%
      select(cell_id, sample_id) %>%
      rename(cell_id_original = cell_id) %>%
      left_join(combined_airr %>%
                  select(cell_id, cell_id_original, sample_id) %>%
                  distinct(),
                by = join_by(cell_id_original, sample_id)) %>%
      pull(cell_id)

    # remove the sample column
    embeddings <- embeddings %>% select(-sample_id)
  } else {
    # since cell_id_original is unique, a named lookup is more efficient than
    # using dplyr's join
    cell_id_lookup <-
      setNames(combined_airr$cell_id, combined_airr$cell_id_original)
    cell_ids <- as.character(cell_id_lookup[embeddings$cell_id])
  }

  # remove the cell_id column
  embeddings <- embeddings %>% select(-cell_id)

  # set column names (by number of dimensions)
  dims <- ncol(embeddings)
  dim_cols <- paste0("Dim", seq_len(dims)) # don't start at 0

  # convert to a sparse matrix and match the regular format (features x cells)
  embeddings_mat <- as.matrix(embeddings)
  rownames(embeddings_mat) <- cell_ids
  colnames(embeddings_mat) <- dim_cols
  embeddings_mat <- Matrix(t(embeddings_mat), sparse = TRUE)

  return(embeddings_mat)
}
