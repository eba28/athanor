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


#' Build a Seurat object from BCR embeddings
#'
#' @description
#' Creates and processes a BCR-only Seurat object from a matrix of pre-computed
#' embeddings (e.g. AntiBERTa2, AntiBERTy, BALM-paired, ESM2, immune2vec, etc.).
#' Produces `bpca`, `BCR.nn`, `BCR_nn`, `BCR_snn`, and `bcr.umap` reductions/graphs.
#'
#' @details
#' Embeddings are used as-is for `scale.data` (no `ScaleData` call) since they
#' are already on a comparable scale. The `data` layer is populated from `counts`
#' so downstream reads do not fail.
#' If many cells share identical embeddings (e.g. from clonal expansion),
#' `RunUMAP()` may hang on spectral initialization. Use
#' `bcr_embeddings_pipeline_dedup()` to handle this automatically.
#'
#' @param embeddings Matrix of BCR embeddings (features x cells).
#' @param embedding_type Character label for the embedding method (stored in `Misc`).
#' @param combined_airr Optional data frame passed to `gex_add_airr()` to add
#'   AIRR metadata columns. If NULL, the step is skipped.
#' @param new_cols Character vector of columns to add from `combined_airr`.
#'   Only used when `combined_airr` is provided.
#' @param num_pcs Number of principal components to compute.
#' @param num_dims Number of PCA dimensions to use for neighbor finding and UMAP.
#' @param k_param Number of nearest neighbors.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @returns A Seurat object with a BCR assay, a new PCA (`bpca`), new neighbor graphs
#'   (`BCR_nn`, `BCR_snn`, `BCR.nn`), and a new UMAP (`bcr.umap`).
#' @export
bcr_embeddings_pipeline <- function(embeddings, embedding_type,
                                    combined_airr = NULL, new_cols = NULL,
                                    num_pcs = 50, num_dims = 20, k_param = 20,
                                    verbose = TRUE) {
  # if embeddings are not already given as a sparse numeric matrix, convert them
  if (!is(embeddings, "dgCMatrix")) {
    embeddings <- as(embeddings, "dgCMatrix")
    cli::cli_inform(c("v" = "Converted embeddings to sparse matrix format."))
  }

  # set up an empty Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = embeddings, assay = "BCR")
  seurat_obj$cell_id <- Seurat::Cells(seurat_obj)
  # DefaultAssay(seurat_obj) <- "BCR" # not needed since it is the only assay
  Seurat::VariableFeatures(seurat_obj) <- rownames(seurat_obj[["BCR"]])

  # embeddings are already normalized and scaled so those steps can be skipped
  seurat_obj <-
    Seurat::SetAssayData(seurat_obj, layer = "data",
                         new.data = GetAssayData(seurat_obj, assay = "BCR",
                                                 layer = "counts"))
  seurat_obj <-
    Seurat::SetAssayData(seurat_obj, layer = "scale.data",
                         new.data = as.matrix(embeddings)) # required format

  # add BCR info to the object's metadata
  # the number of heavy chains should match how many sequences ran through
  # (except for immune2vec)
  if (!is.null(combined_airr)) {
    seurat_obj <- gex_add_airr(seurat_obj, combined_airr = combined_airr,
                               new_cols = new_cols, verbose = verbose)
  }

  max_dim <- min(nrow(embeddings), ncol(embeddings))
  max_pcs <- max_dim - 1
  if (num_pcs >= max_dim) {
    warning("num_pcs (", num_pcs, ") >= embedding dimensions (", max_dim,
            "); reducing to ", max_pcs, " to avoid SVD convergence failure.")
    num_pcs <- max_pcs
    num_dims <- min(num_dims, max_pcs)
  }

  # irlba throws a warning to "use a standard svd instead" when requesting more
  # than 50% of all singular values, so use exact SVD in that case
  use_approx <- num_pcs < max_dim / 2
  seurat_obj <- Seurat::RunPCA(object = seurat_obj, npcs = num_pcs,
                               reduction.name = "bpca", reduction.key = "bpca_",
                               approx = use_approx, verbose = verbose)
  cli::cli_inform(c("v" = "Computed PCA with {num_pcs} dimensions using {ifelse(use_approx, 'approximate', 'exact')} SVD."))

  # visual check for how many PCs to keep (aka num_dims)
  # ElbowPlot(seurat_obj, reduction = "bpca", ndims = num_pcs)

  seurat_obj <- Seurat::FindNeighbors(object = seurat_obj, reduction = "bpca",
                                      dims = 1:num_dims, k.param = k_param,
                                      return.neighbor = TRUE,
                                      graph.name = "BCR.nn", verbose = verbose)
  seurat_obj <- Seurat::FindNeighbors(object = seurat_obj, reduction = "bpca",
                                      dims = 1:num_dims, k.param = k_param,
                                      compute.SNN = TRUE,
                                      graph.name =
                                        str_c("BCR_", c("", "s"), "nn"),
                                      verbose = verbose)
  seurat_obj <- Seurat::RunUMAP(object = seurat_obj, reduction = "bpca",
                                nn.name = "BCR.nn", n.neighbors = k_param,
                                reduction.name = "bcr.umap",
                                reduction.key = "bcrUMAP_",
                                verbose = verbose)

  # add useful information to the Miscellaneous slot
  Seurat::Misc(seurat_obj, slot = "embedding_type") <- embedding_type
  Seurat::Misc(seurat_obj, slot = "embedding_dims") <- nrow(embeddings)

  # print an overview of what was done
  if (verbose) {
    cli::cli_inform(c("v" = "Created BCR assay with {nrow(embeddings)} features and {ncol(embeddings)} cells."))
    cli::cli_inform(c("v" = "Computed PCA with {num_pcs} dimensions, using {num_dims} for neighbor finding and UMAP."))
    cli::cli_inform(c("v" = "Found neighbors with k = {k_param}."))
  }

  seurat_obj
}


#' Build a Seurat object from BCR embeddings, deduplicating identical embeddings first
#'
#' @description
#' Wrapper around `bcr_embeddings_pipeline()` that detects cells with identical
#' embeddings (e.g. from clonal expansion), runs PCA, neighbor finding, and UMAP
#' on unique embeddings only, then copies coordinates back to all cells. This
#' avoids the spectral initialization hang in `RunUMAP()` caused by degenerate
#' neighbor graphs from zero-distance duplicate points.
#'
#' @inheritParams bcr_embeddings_pipeline
#' @returns A Seurat object as returned by `bcr_embeddings_pipeline()`, with all
#'   original cells present. Cells with identical embeddings receive the same PCA
#'   and UMAP coordinates as their first occurrence.
#' @export
bcr_embeddings_pipeline_dedup <- function(embeddings, embedding_type,
                                          combined_airr = NULL, new_cols = NULL,
                                          num_pcs = 50, num_dims = 20, k_param = 20,
                                          verbose = TRUE) {
  if (!is(embeddings, "dgCMatrix")) {
    embeddings <- as(embeddings, "dgCMatrix")
    cli::cli_inform(c("v" = "Converted embeddings to sparse matrix format."))
  }

  # build per-cell string keys from the sparse structure
  cli::cli_inform(c("i" = "Scanning {ncol(embeddings)} cells for identical embeddings — this may take a moment for large datasets."))
  emb_keys <- vapply(seq_len(ncol(embeddings)), function(j) {
    start <- embeddings@p[j] + 1L
    end   <- embeddings@p[j + 1L]
    if (start > end) return("")
    paste(embeddings@i[start:end], embeddings@x[start:end], sep = "=", collapse = ",")
  }, character(1))
  # canonical_idx[i] = index in embeddings of the first occurrence of cell i's embedding
  canonical_idx <- match(emb_keys, emb_keys)
  unique_mask <- canonical_idx == seq_len(ncol(embeddings))
  n_dups <- sum(!unique_mask)

  if (n_dups == 0) {
    cli::cli_inform(c("v" = "No identical embeddings found; running standard pipeline."))
    return(bcr_embeddings_pipeline(embeddings, embedding_type,
                                   combined_airr = combined_airr, new_cols = new_cols,
                                   num_pcs = num_pcs, num_dims = num_dims,
                                   k_param = k_param, verbose = verbose))
  }

  cli::cli_inform(c("!" = "{n_dups} cell{?s} have identical embeddings out of {ncol(embeddings)} total cells."))
  cli::cli_inform(c("v" = "Running pipeline on {sum(unique_mask)} unique embeddings; duplicates will receive copied coordinates."))

  emb_unique <- embeddings[, unique_mask, drop = FALSE]
  # for each cell, its row index in seurat_unique
  unique_positions <- integer(ncol(embeddings))
  unique_positions[which(unique_mask)] <- seq_len(sum(unique_mask))
  unique_row_idx <- unique_positions[canonical_idx]

  # run core pipeline on unique embeddings only (metadata handled below on full object)
  seurat_unique <- bcr_embeddings_pipeline(emb_unique, embedding_type,
                                           combined_airr = NULL,
                                           num_pcs = num_pcs, num_dims = num_dims,
                                           k_param = k_param, verbose = verbose)

  # build full Seurat object with all cells
  seurat_obj <- Seurat::CreateSeuratObject(counts = embeddings, assay = "BCR")
  seurat_obj$cell_id <- Seurat::Cells(seurat_obj)
  Seurat::VariableFeatures(seurat_obj) <- rownames(seurat_obj[["BCR"]])
  seurat_obj <-
    Seurat::SetAssayData(seurat_obj, layer = "data",
                         new.data = Seurat::GetAssayData(seurat_obj, assay = "BCR",
                                                         layer = "counts"))
  seurat_obj <-
    Seurat::SetAssayData(seurat_obj, layer = "scale.data",
                         new.data = as.matrix(embeddings))

  # expand PCA: duplicate cells get the same coordinates as their canonical cell
  pca_emb <- Seurat::Embeddings(seurat_unique, "bpca")[unique_row_idx, , drop = FALSE]
  rownames(pca_emb) <- colnames(embeddings)
  seurat_obj[["bpca"]] <-
    Seurat::CreateDimReducObject(embeddings = pca_emb,
                                 loadings = Seurat::Loadings(seurat_unique, "bpca"),
                                 stdev = seurat_unique[["bpca"]]@stdev,
                                 key = "bpca_", assay = "BCR")

  # FindNeighbors on full object: identical cells form natural cliques
  seurat_obj <- Seurat::FindNeighbors(object = seurat_obj, reduction = "bpca",
                                      dims = 1:num_dims, k.param = k_param,
                                      return.neighbor = TRUE,
                                      graph.name = "BCR.nn", verbose = verbose)
  seurat_obj <- Seurat::FindNeighbors(object = seurat_obj, reduction = "bpca",
                                      dims = 1:num_dims, k.param = k_param,
                                      compute.SNN = TRUE,
                                      graph.name = str_c("BCR_", c("", "s"), "nn"),
                                      verbose = verbose)

  # expand UMAP: same logic as PCA above
  umap_emb <- Seurat::Embeddings(seurat_unique, "bcr.umap")[unique_row_idx, , drop = FALSE]
  rownames(umap_emb) <- colnames(embeddings)
  seurat_obj[["bcr.umap"]] <-
    Seurat::CreateDimReducObject(embeddings = umap_emb,
                                 key = "bcrUMAP_", assay = "BCR")

  # add BCR info to metadata on the full object
  # the number of heavy chains should match how many sequences ran through
  # (except for immune2vec)
  if (!is.null(combined_airr)) {
    seurat_obj <- gex_add_airr(seurat_obj, combined_airr = combined_airr,
                               new_cols = new_cols, verbose = verbose)
  }

  Seurat::Misc(seurat_obj, slot = "num_dups") <- n_dups
  Seurat::Misc(seurat_obj, slot = "embedding_dims") <- nrow(embeddings)
  Seurat::Misc(seurat_obj, slot = "embedding_type") <- embedding_type

  if (verbose) {
    cli::cli_inform(c("v" = "Created BCR assay with {nrow(embeddings)} features and {ncol(embeddings)} cells ({n_dups} duplicate{?s} assigned copied coordinates)."))
    cli::cli_inform(c("v" = "Computed PCA with {num_pcs} dimensions, using {num_dims} for neighbor finding and UMAP."))
    cli::cli_inform(c("v" = "Found neighbors with k = {k_param}."))
  }

  seurat_obj
}


#'
#' Bin the mutation frequency
#'
#' @description
#' Bins the `mu_freq` column in a Seurat object into specified categories (e.g. 0%, 0-1%, 1-5%, etc.) for easier visualization and analysis.
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


#' Convert the output of an embedding method to a matrix
#'
#' @description
#' This function takes the output from an embedding method (e.g. AntiBERTy, ESM2, immune2vec) and converts it into a `Matrix` matrix format that can be used for downstream analysis.
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
    cli::cli_abort"One or both of the required columns (`cell_id`, `cell_id_original`) ",
         "are missing from combined_airr.")
  }

  # immune2vec only
  if (!"cell_id" %in% colnames(embeddings)) {
    # just in case
    if (rlang::is_missing(combined_airr_input)) {
      cli::cli_abort"Please provide the original file used as input for running immune2vec.")
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
  cli::cli_inform("{length(embeddings_only_cell_ids)} cell{?s} in the embedding with no \\
corresponding cell IDs in the combined AIRR file were removed.")

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

  # convert to a sparse matrix and match the regular format (features by cells)
  embeddings_mat <- as.matrix(embeddings)
  rownames(embeddings_mat) <- cell_ids
  colnames(embeddings_mat) <- dim_cols
  embeddings_mat <- Matrix(t(embeddings_mat), sparse = TRUE)

  return(embeddings_mat)
}


#' Convert family and gene information to sorted factors
#'
#' @description
#' Converts the gene family and gene columns added by [add_family_info()] to
#' properly ordered factors using numeric sorting.
#'
#' @details
#' For after [add_family_info()] has been run
#'
#' @param combined_airr An AIRR-formatted data.frame. This function ensures that gene
#' families and genes are ordered correctly (e.g. IGHV1, IGHV2, IGHV10 instead
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

  # TODO: make this more efficient
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
    cli::cli_abort"No files found, are you sure you're using the correct directory?")
  }
  combined_bcr <- purrr::map(rep_files, read_tsv, show_col_types = FALSE,
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
      mutate(isotype = case_when(stringr::str_detect(c_call, "IGHA") ~ "IgA",
                                 stringr::str_detect(c_call, "IGHD") ~ "IgD",
                                 stringr::str_detect(c_call, "IGHE") ~ "IgE",
                                 stringr::str_detect(c_call, "IGHG") ~ "IgG",
                                 stringr::str_detect(c_call, "IGHM") ~ "IgM"))
  }

  # filter out heavy chains with unassigned isotypes
  cli::cli_inform("{nrow(filter(combined_bcr, locus == 'IGH', is.na(c_call)))} \\
heavy chain{?s} with {.code NA} c_calls.")
  combined_bcr <- combined_bcr %>% filter(!is.na(c_call) | locus != "IGH")

  # filter out cells with unpaired light chains
  # airrflow does this, but we have to do it again since we now filtered out
  # more heavy chains
  unpaired_light_chains <- combined_bcr %>%
    group_by(cell_id) %>%
    filter(!any(locus == "IGH")) %>%
    pull(cell_id)

  cli::cli_inform("{length(unpaired_light_chains)} unpaired light chain{?s} \\
(no corresponding heavy chain with the same cell ID).")
  combined_bcr <- combined_bcr %>% filter(!cell_id %in% unpaired_light_chains)

  # filter out heavy chains with light chains assigned to the isotype
  # this comes from the "multi" column in the original 10x data
  # see: https://kb.10xgenomics.com/s/article/360001715971-What-is-the-source-of-Multi-chains-in-the-contig-annotations-csv-files
  multi_heavy_chains <- combined_bcr %>%
    filter(locus == "IGH", is.na(isotype)) %>%
    pull(cell_id)

  cli::cli_inform("{length(multi_heavy_chains)} heavy chain{?s} with light chain c_calls.")
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

  cli::cli_inform("{nrow(combined_bcr)} total chain{?s} in the combined data.")

  return(combined_bcr)
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
#'   2. Convert ordered factors to ordinal scores OR one-hot encode categorical variables
#'   3. Center and normalize all numeric predictors
#'   4. Transpose the result for compatibility with Seurat assays
#'  You could do something like janitor::clean_names(bcr_features) to remove the underscores from the column names right off the bat, but one-hot encoding will add names with underscores automatically unless you messing with the `naming` argument in `step_dummy()`.
#'
#' @param bcr_features A data frame containing BCR features (e.g. isotype, mutation
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
  has_ordered <- any(sapply(bcr_features, is.ordered))

  # rename numeric columns
  bcr_features <- bcr_features %>%
                    rename_with(~str_c(., "-scaled"), where(is.numeric))

  # rename (to distinguish from existing metadata) and convert the columns
  if (has_ordered) {
    cli::cli_inform(c("i" = "Ordered variables detected: converting with ordinal scoring"))
    bcr_features <- bcr_features %>%
                      rename_with(~str_c(., "-ordered"), where(is.ordered))

    # convert ordered variables to numeric; one-hot encode any remaining nominals
    ref_cell <- recipe( ~ ., data = bcr_features) %>%
                  step_ordinalscore(all_ordered_predictors()) %>%
                  step_unknown(all_nominal_predictors()) %>%
                  step_dummy(all_nominal_predictors(), one_hot = TRUE)
  } else {
    cli::cli_inform(c("i" = "No ordered variables: applying one-hot encoding to nominal predictors"))

    # one-hot encoding of categorical variables
    ref_cell <- recipe( ~ ., data = bcr_features) %>%
                  # cover missing values just in case
                  step_unknown(all_nominal_predictors()) %>%
                  step_dummy(all_nominal_predictors(), one_hot = TRUE)
  }

  # step_zv removes zero-variance columns (e.g. uniform "unknown" levels from step_unknown)
  # before normalization, which would produce NaN/Inf for them
  ref_cell <- ref_cell %>%
                step_zv(all_numeric_predictors()) %>%
                step_normalize(all_numeric_predictors()) %>%
                prep(training = bcr_features)

  n_removed <- length(ref_cell$steps[[which(sapply(ref_cell$steps, \(s) inherits(s, "step_zv")))]]$removals)
  if (n_removed > 0) {
    cli::cli_inform(c("i" = "Removed {n_removed} zero-variance feature{?s}"))
  }

  bcr_features <- bake(ref_cell, new_data = NULL) %>% base::t()

  cli::cli_inform(c("v" = "Processed {nrow(bcr_features)} feature{?s} across {ncol(bcr_features)} cell{?s}"))

  # step_dummy uses underscores for level names; Seurat requires dashes in feature names
  rownames(bcr_features) <- gsub("_", "-", rownames(bcr_features))

  return(bcr_features)
}
