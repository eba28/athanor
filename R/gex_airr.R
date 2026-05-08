#' Add AIRR information to a Seurat object
#'
#' @description
#' This function integrates adaptive immune receptor repertoire (AIRR) data with
#' gene expression (GEX) data in a Seurat object. Currently built to integrate
#' BCR data, including both heavy and light chain information.
#'
#' @details
#' This would typically be used after [seurat_pipeline()] and before [concatenate_gex_bcr()].
#' Right now this is just built to integrate in BCR data and assumes that the BCR data includes light chains.
#' Assumes that `seurat_obj` contains `cell_id` and `annotated_clusters_simpler`.
#' Collapsed light chains are alphabetized.
#' There are columns for IGK and IGL instead of just "light" since a cell can have both.
#' Adds mutation frequency bins and isotype staging information for BCR data.
#'
#' @param seurat_obj The Seurat object containing GEX data.
#' @param airr_type Type of immune receptor data. Currently supports "BCR".
#' @param combined_airr BCR AIRR formatted data frame with heavy and light chains.
#' @param new_cols Vector of column names to select from the AIRR data.
#' @param verbose Logical indicating whether or not to print integration summary information.
#'
#' @returns The Seurat object with AIRR columns added to the metadata, including
#'   Has_BCR, isotype information, mutation frequencies, and pairing status.
#' @export
gex_add_airr <- function(seurat_obj, airr_type = "BCR",
                         combined_airr, new_cols, verbose = TRUE) {
  # TODO: add back a more detailed printout of the light chains e.g.
  # There are 118563 heavy chain BCRs being integrated in, of which 111039 have one paired light chain, 6243 have multiple paired light chains, and 1281 have no paired light chains.

  # TODO: this can also be run to add BCR info to a BCR object, so rename the
  # function and don't have it say GEX (look at the assay name?)

  # make sure there isn't a format mismatch between cell ids
  if (length(intersect(seurat_obj$cell_id, combined_airr$cell_id)) == 0) {
    cli::cli_abort("Please check the cell id format in the Seurat object or AIRR table.")
  }

  if (any(!new_cols %in% colnames(combined_airr))) {
    cli::cli_abort("Make sure that you are only adding columns which exist in the AIRR table.")
  }

  # for comparison later on
  meta_cols_orig <- colnames(seurat_obj[[]])

  # define names and parameters based on AIRR type
  airr_col <- paste0("Has_", airr_type)
  airr_type <- paste0(airr_type, "s")

  if (airr_type == "BCRs") {
    locus_main <- "IGH"
    chain_main <- "heavy"
    chain_other <- "light"
    paired_col <- paste0("paired_", chain_other)
    relocate_v <- c("v_call_family_IGK", "v_call_family_IGL")
    relocate_j <- c("j_call_family_IGK", "j_call_family_IGL")
  } else { # TCRs
    locus_main <- "TRB"
    chain_main <- "beta"
    chain_other <- "alpha"
    paired_col <- paste0("paired_", chain_other)
    relocate_v <- "v_call_family_TRA"
    relocate_j <- "j_call_family_TRA"
  }

  # Seurat only takes unique cell ids, hence the filter to heavy/beta chains
  cells_main <- (filter(combined_airr, locus == locus_main))$cell_id
  cells_other <- (filter(combined_airr, locus != locus_main))$cell_id
  cells_paired <- intersect(cells_main, cells_other) # cells_main[which(cells_main %in% cells_other)]

  # select columns of interest from the AIRR data
  combined_airr_select_main <-
    combined_airr %>%
    filter(locus == locus_main) %>%
    select(all_of(new_cols)) %>%
    mutate(!!airr_col := TRUE, !!paired_col := cell_id %in% cells_paired,
           .before = cell_id)

  # pick out some other chain info (note that there may be multiple chains)
  combined_airr_select_other <-
    combined_airr %>%
    filter(locus != locus_main) %>%
    select(cell_id, locus, v_call_family, j_call_family) %>%
    pivot_wider(names_from = locus,
                values_from = c(v_call_family, j_call_family),
                values_fn = ~ paste(sort(.), collapse = ", "))

  # merge the main and other chain information
  combined_airr_select <-
    left_join(combined_airr_select_main, combined_airr_select_other,
              by = "cell_id") %>%
    relocate(all_of(relocate_v), .before = v_call_gene) %>%
    relocate(all_of(relocate_j), .before = j_call_gene)

  # make select columns into alphabetized factors (useful for plotting later)
  factor_cols <- intersect(c("paired_light", "paired_alpha", "isotype",
                             "c_call", "locus"), names(combined_airr_select))
  if (length(factor_cols) > 0) {
    combined_airr_select <- combined_airr_select %>%
      mutate(across(all_of(factor_cols), ~ factor(., levels = sort(unique(.)))))
  }

  # for columns that have numbers in them
  call_cols <- str_subset(names(combined_airr_select), "call_")
  if (length(call_cols) > 0) {
    combined_airr_select <-
      combined_airr_select %>%
      mutate(across(all_of(call_cols),
                    ~ factor(., levels = str_sort(unique(.), numeric = TRUE))))
  }

  # combine the two sources of barcodes using the unique cell ids
  combined_airr_select_obj <-
    left_join(seurat_obj[["cell_id"]], combined_airr_select,
              by = "cell_id") %>% # multiple = "error"
    # fill in FALSE for any missing values (i.e. those that don't have AIRR data)
    mutate(!!airr_col := replace_na(.data[[airr_col]], FALSE))
  # select(-cell_id)

  # check that cell counts match
  # ncol(seurat_obj@assays$RNA) == nrow(combined_airr_select_obj)

  # integrate AIRR data with the GEX Seurat object
  # you could also just bind_cols with seurat_obj[[]]
  # for (meta_col in colnames(combined_airr_select_obj)) {
  #  seurat_obj[[meta_col]] <- combined_airr_select_obj[[meta_col]]
  # }
  seurat_obj <- AddMetaData(seurat_obj, metadata = combined_airr_select_obj)

  # add more AIRR-specific columns
  if (airr_type == "BCRs") {
    # these columns will be useful later when we are exploring the JC results
    # this could be (and are) within run_wnn(), but it's nice
    # to have these columns when just examining the GEX object
    seurat_obj <- bin_mu_freq(seurat_obj)

    seurat_obj@meta.data <-
      seurat_obj[[]] %>%
      mutate(isotype_stage =
               case_when(isotype %in% c("IgD", "IgM") ~ "Unswitched",
                         isotype %in% c("IgA", "IgE", "IgG") ~ "Switched",
                         .default = NA),
             mu_freq_iso = paste0(mu_freq_bins_fewer, " ", isotype_stage)) %>%
      mutate(locus_light = case_when(!is.na(v_call_family_IGK) &
                                       !is.na(v_call_family_IGL) ~ "IGK, IGL",
                                     !is.na(v_call_family_IGK) ~ "IGK",
                                     !is.na(v_call_family_IGL) ~ "IGL",
                                     .default = NA),
             .after = locus) %>%
      mutate(mu_freq_iso = na_if(mu_freq_iso, "NA NA")) %>%
      # helps with plotting
      mutate(isotype_stage = factor(isotype_stage),
             mu_freq_iso =
               factor(mu_freq_iso,
                      levels = c("0% Unswitched", "0% Switched",
                                 "0% to 1% Unswitched", "0% to 1% Switched",
                                 ">1% Unswitched", ">1% Switched")))

    # let's add a few more approaches to cell typing
    # TODO: decapitalize these?
    seurat_obj@meta.data <-
      seurat_obj[[]] %>%
      mutate(annotated_clusters_bcr =
               ifelse(mu_freq <= 0.01 & isotype %in% c("IgM", "IgD"),
                      "Naive B Cells", "Non-Naive B Cells"))

    if ("annotated_clusters_simpler" %in% colnames(seurat_obj[[]])) {
      # assumes that the plasma cells don't need re-assigning
      seurat_obj@meta.data <-
        seurat_obj[[]] %>%
        mutate(annotated_clusters_gex_bcr =
                 ifelse(annotated_clusters_simpler == "Naive B cells" &
                          (mu_freq > 0.01 | isotype %in% c("IgA", "IgE", "IgG")),
                        "Memory B cells", annotated_clusters_simpler))
    } else {
      cli::cli_inform(c("!" = "{.code annotated_clusters_gex_bcr} not added: \\
        {.code annotated_clusters_simpler} column not found in metadata."))
    }
  } else if (airr_type == "TCRs") {
    # classify cells by TCR chain type (alpha-beta vs gamma-delta)
    has_tra <- "v_call_family_TRA" %in% names(seurat_obj[[]])
    has_trg <- "v_call_family_TRG" %in% names(seurat_obj[[]])

    if (has_tra || has_trg) {
      seurat_obj@meta.data <-
        seurat_obj[[]] %>%
        mutate(tcr_chain_type =
                 case_when(has_tra & !is.na(v_call_family_TRA) ~ "Alpha-Beta",
                           has_trg & !is.na(v_call_family_TRG) ~ "Gamma-Delta",
                           .default = NA),
               .after = locus)
    }
  }

  # make sure that the levels are okay
  seurat_obj@meta.data <- seurat_obj@meta.data %>% droplevels()

  # print a summary of the integration
  if (verbose) {
    count_gex_airr <- nrow(filter(combined_airr_select, get(airr_col)))
    count_paired <- nrow(filter(combined_airr_select, get(paired_col) == "TRUE"))
    count_filtered <- nrow(filter(combined_airr_select_obj, get(airr_col)))
    count_unmatched <- count_gex_airr - count_filtered
    coverage <- mean(Cells(seurat_obj) %in% combined_airr_select$cell_id)
    count_excluded <- 1 - mean(combined_airr_select$cell_id %in% Cells(seurat_obj))

    new_cols <- str_c(setdiff(colnames(seurat_obj[[]]), meta_cols_orig),
                      collapse = ", ")

    cli::cli_inform(c(
      "i" = "There are {count_gex_airr} {chain_main} chain {airr_type} being integrated, \\
{count_paired} of which have at least one paired {chain_other} chain.",
      "i" = "When filtered to matching cell IDs: {count_filtered} cells \\
({label_percent(accuracy = 0.1)(coverage)} of GEX data ({ncol(seurat_obj)} cells)).",
      "i" = "{count_unmatched} {airr_type} \\
({label_percent(accuracy = 0.1)(count_excluded)} of total {chain_main} chains) \\
had no matching GEX cell IDs and were excluded.",
      "i" = "New metadata columns added: {new_cols}."))

    # detailed chain pairing breakdown
    if (airr_type == "BCRs") {
      has_igk <- !is.na(combined_airr_select$v_call_family_IGK)
      has_igl <- !is.na(combined_airr_select$v_call_family_IGL)
      count_one_light  <- sum(xor(has_igk, has_igl))
      count_multi_light <- sum(has_igk & has_igl)
      count_no_light   <- sum(!has_igk & !has_igl)
      cli::cli_inform(c(
        "i" = "Of the {count_gex_airr} heavy chains: {count_one_light} have one paired \\
light chain, {count_multi_light} have multiple paired light chains, \\
and {count_no_light} have no paired light chain."))
    } else if (airr_type == "TCRs" &&
               "v_call_family_TRA" %in% names(combined_airr_select)) {
      has_tra <- !is.na(combined_airr_select$v_call_family_TRA)
      count_one_alpha   <- sum(has_tra & !str_detect(combined_airr_select$v_call_family_TRA,
                                                     ", ", negate = FALSE))
      count_multi_alpha <- sum(has_tra & str_detect(combined_airr_select$v_call_family_TRA,
                                                    ", "))
      count_no_alpha    <- sum(!has_tra)
      cli::cli_inform(c(
        "i" = "Of the {count_gex_airr} beta chains: {count_one_alpha} have one paired \\
alpha chain, {count_multi_alpha} have multiple paired alpha chains, \\
and {count_no_alpha} have no paired alpha chain."))
    }
  }

  # return the updated Seurat object
  seurat_obj
}


#' Merge a GEX Seurat object with a BCR Seurat object
#'
#' @description
#' Combines a gene expression (GEX) Seurat object and a BCR Seurat object by
#' matching on shared `cell_id` values. Handles the common case where the two
#' objects have different cell counts (e.g. not every GEX cell has a BCR sequence).
#'
#' @details
#' Cell matching is done via `cell_id` metadata, not raw Seurat barcodes.
#' If GEX barcodes differ from `cell_id` (e.g. they carry a sample prefix),
#' BCR cell names are renamed via [SeuratObject::RenameCells()] to match before the
#' assay is transferred.
#'
#' Both objects are subset to their shared cells; the result is suitable for concatenation or WNN or other downstream workflows.
#'
#' Due to subsetting, graphs, neighbors and reductions maybe have to be regenerated.
#'
#' @param gex_obj A Seurat object containing GEX (RNA) data with a `cell_id`
#'   metadata column.
#' @param bcr_obj A Seurat object containing a BCR assay with a `cell_id`
#'   metadata column (typically produced by [bcr_embeddings_pipeline()]).
#' @param transfer_reductions Whether to copy BCR reductions (`bpca`,
#'   `bcr.umap`) and graphs (`BCR.nn`, `BCR_nn`, `BCR_snn`) from `bcr_obj`
#'   into the merged object.
#' @param verbose Logical indicating whether or not to print a summary of the merge.
#'
#' @returns A Seurat object with both RNA and BCR assays and BCR metadata columns.
#' @export
merge_gex_bcr <- function(gex_obj, bcr_obj, transfer_reductions = TRUE,
                          verbose = TRUE) {
  # argument checks
  if (!inherits(gex_obj, "Seurat")) cli::cli_abort("gex_obj must be a Seurat object.")
  if (!inherits(bcr_obj, "Seurat")) cli::cli_abort("bcr_obj must be a Seurat object.")
  if (!"cell_id" %in% colnames(gex_obj[[]])) {
    cli::cli_abort("The GEX object must have a 'cell_id' metadata column.")
  }
  if (!"cell_id" %in% colnames(bcr_obj[[]])) {
    cli::cli_abort("The BCR object must have a 'cell_id' metadata column.")
  }
  if (!"BCR" %in% names(bcr_obj@assays)) {
    cli::cli_abort("The BCR object must contain a BCR assay.")
  }

  shared_cell_ids <- intersect(gex_obj$cell_id, bcr_obj$cell_id)
  if (length(shared_cell_ids) == 0) {
    cli::cli_abort("No shared cell ids found; please check 'cell_id' formatting in both objects.")
  }

  # original cell counts and BCR info
  n_gex <- ncol(gex_obj)
  n_bcr <- ncol(bcr_obj)
  n_shared <- length(shared_cell_ids)

  # TODO: use inferring functions
  bcr_dims <- bcr_obj@neighbors$BCR.nn@alg.info$ndim
  bcr_k <- ncol(bcr_obj@neighbors[["BCR.nn"]]@nn.idx)

  if (verbose) {
    cov_gex <- scales::label_percent(accuracy = 0.1)(n_shared / n_gex)
    cov_bcr <- scales::label_percent(accuracy = 0.1)(n_shared / n_bcr)
    cli::cli_inform(c("i" = "{n_shared} total shared cells \\
                      ({cov_gex} of the {n_gex} GEX cells and {cov_bcr} \\
                      of the {n_bcr} BCR cells)."))
  }

  # subset the objects to just the shared cells
  gex_obj <- subset(gex_obj, cell_id %in% shared_cell_ids)
  bcr_obj <- subset(bcr_obj, cell_id %in% shared_cell_ids)

  # add the BCR assay to the GEX object
  seurat_obj <- gex_obj
  seurat_obj[["BCR"]] <- bcr_obj[["BCR"]]

  # transfer over some of the BCR-specific metadata
  # all of the other column should already be present in the GEX object
  seurat_obj@meta.data <-
    left_join(seurat_obj[[]],
              bcr_obj[[]] %>% select(cell_id, nCount_BCR, nFeature_BCR),
              by = join_by(cell_id)) %>%
    {if ("nFeature_ADT" %in% names(.))
      relocate(., nCount_BCR, nFeature_BCR, .after = nFeature_ADT) else .}
  rownames(seurat_obj@meta.data) <- seurat_obj$cell_id

  # filter out cells that have no BCR data
  n_missing <- sum(!seurat_obj$Has_BCR)
  seurat_obj <- subset(seurat_obj, Has_BCR == TRUE)
  bcr_obj <- subset(bcr_obj, Has_BCR == TRUE)

  if (verbose) cli::cli_inform(c("!" = "{n_missing} cells had no BCR metadata and were excluded."))

  # filter out cells without paired light chains
  # have to use quotes because it is a factor
  n_unpaired <- sum(seurat_obj$paired_light == "FALSE")
  seurat_obj <- subset(seurat_obj, paired_light == "TRUE")
  bcr_obj <- subset(bcr_obj, paired_light == "TRUE")

  if (verbose) cli::cli_inform(c("!" = "{n_unpaired} cells had no paired light chain{?s} and were excluded."))

  # just in case
  seurat_obj@meta.data <- seurat_obj@meta.data %>% droplevels()

  # transfer BCR reductions and graphs
  # assumes standard names for the BCR object
  if (transfer_reductions) {
    for (reduc in c("bpca", "bcr.umap")) {
      if (reduc %in% names(bcr_obj@reductions)) {
        seurat_obj@reductions[[reduc]] <- bcr_obj@reductions[[reduc]]
      }
    }
    for (graph in c("BCR_nn", "BCR_snn")) {
      if (graph %in% names(bcr_obj@graphs)) {
        seurat_obj@graphs[[graph]] <- bcr_obj@graphs[[graph]]
      }
    }
    # note that if bcr_obj was subset above, there will be no neighbors to transfer
    for (nn in "BCR.nn") {
      if (nn %in% names(bcr_obj@neighbors)) {
        seurat_obj@neighbors[[nn]] <- bcr_obj@neighbors[[nn]]
      } else {
        # don't need all of regen_reduc()
        # TODO: do this on the bcr_obj?
        seurat_obj@neighbors[[nn]] <-
          FindNeighbors(seurat_obj, reduction = "bpca",
                        dims = 1:bcr_dims, k.param = bcr_k,
                        graph.name = stringr::str_c("BCR", "_", c("", "s"),
                                                    "nn"),
                        verbose = verbose)
      }
    }
  }

  # transfer the other slots
  seurat_obj@commands <- c(seurat_obj@commands, bcr_obj@commands)
  seurat_obj@misc <- c(seurat_obj@misc, bcr_obj@misc)

  # TODO: print out the reductions
  if (verbose) {
    n_gex_only <- n_gex - n_shared
    n_bcr_only <- n_bcr - n_shared

    if (n_gex_only > 0) {
      cli::cli_inform(c("!" = "{n_gex_only} GEX cell{?s} had no BCR matches and were excluded."))
    }
    if (n_bcr_only > 0) {
      cli::cli_inform(c("!" = "{n_bcr_only} BCR cell{?s} had no GEX matches and were excluded."))
    }

    cli::cli_inform(c("v" = "The final merged object has {ncol(seurat_obj)} cells in both RNA and BCR assays."))
  }

  seurat_obj
}


#' Infer k_param from existing RNA neighbors, or return a default.
#'
#' @description
#' Infers the k_param used for neighbor finding from the existing RNA neighbors slot, or
#' returns a default if not found.
#'
#' @param seurat_obj A Seurat object containing a neighbors slot.
#' @param default The default k_param to return if not found in the neighbors slot.
#' @param verbose Logical indicating whether or not to print a message about the source of the
#' inferred k_param.
#'
#' @returns An integer representing the k_param used for neighbor finding.
infer_k_param <- function(seurat_obj, default = 20, verbose = TRUE) {
  # # TODO: only print out the actual values that are filled in
  # if (missing(num_dims) | missing(k_param)) {
  #   # try to use the existing neighbors slot if possible
  #   if (nn_name %in% names(seurat_obj@neighbors)) {
  #     nn <- seurat_obj@neighbors[[nn_name]]
  #
  #     if (missing(num_dims)) num_dims <- nn@alg.info$ndim
  #     if (missing(k_param)) k_param <- ncol(nn@nn.idx)
  #
  #     cli::cli_inform("Using existing neighbor graph for assay {assay} to determine num_dims ({num_dims}) and k_param ({k_param}).")
  #   } else if (nn_command %in% names(seurat_obj@commands)) {
  #     cmd <- seurat_obj@commands[[nn_command]]
  #
  #     if (missing(num_dims)) num_dims <- length(cmd$dims)
  #     if (missing(k_param)) k_param <- cmd$k.param
  #
  #     cli::cli_inform("Using existing command {nn_command} for assay {assay} to determine num_dims ({num_dims}) and k_param ({k_param}).")
  #   }
  #   else {
  #     cli::cli_inform(c("i" = "No existing neighbor graph found for assay {assay}, \\
  #                              so using default values for num_dims (20) and k_param (20)."))
  #     if (missing(num_dims)) num_dims <- 20
  #     if (missing(k_param)) k_param <- 20
  #   }
  # }

  has_nn <- "RNA.nn" %in% names(seurat_obj@neighbors)
  k <- if (has_nn) ncol(seurat_obj@neighbors[["RNA.nn"]]@nn.idx) else as.integer(default)
  if (verbose) {
    source_msg <- if (has_nn) "RNA neighbors" else "default"
    cli::cli_inform(c("i" = "Using k_param = {k} from {source_msg}."))
  }
  k
}

#' Infer num_dims (PCA dims used for neighbor finding) from existing RNA neighbors.
#'
#' @description
#' Infers the number of dimensions used for neighbor finding from the existing RNA neighbors slot, or returns a default if not found.
#'
#' @param seurat_obj A Seurat object containing a neighbors slot.
#' @param default The default number of dimensions to return if not found in the neighbors slot.
#' @param verbose Logical indicating whether or not to print a message about the source of the inferred num_dims.
#'
#' @returns An integer representing the number of dimensions used for neighbor finding.
infer_num_dims <- function(seurat_obj, default = 20, verbose = TRUE) {
  has_nn <- "RNA.nn" %in% names(seurat_obj@neighbors)
  d <- if (has_nn) seurat_obj@neighbors[["RNA.nn"]]@alg.info$ndim else as.integer(default)
  if (verbose) {
    source_msg <- if (has_nn) "existing neighbors" else "default"
    cli::cli_inform(c("i" = "Using num_dims = {d} from {source_msg}."))
  }
  d
}


#' Concatenate GEX and BCR data in a Seurat object
#'
#' @description
#' Creates a combined representation of gene expression (GEX) and B-cell
#' receptor (BCR) data, then runs the standard Seurat pipeline on the result.
#' Supports two input types and two PCA stages — see Details.
#'
#' @details
#' This would typically be used after [seurat_pipeline()] and [gex_add_airr()].
#'
#' The six combinations of `input_type` and `pca_stage` offer different
#' tradeoffs:
#'
#' **`input_type = "embeddings"`, `pca_stage = "raw"`** \cr
#' Same as above but BCR data comes from a pre-computed embedding matrix or the
#' `BCR` assay of a merged object (from [merge_gex_bcr()]). Supply `embeddings`
#' directly or pass a merged object and the BCR assay is detected automatically.
#' Embedding dimensions tend to have more comparable scales to RNA features than
#' raw metadata columns, but scale mismatch still applies.
#'
#' **`input_type = "features"`, `pca_stage = "raw"`** (default) \cr
#' BCR metadata columns (via `cols_to_include`) are processed by
#' [process_bcr_features()] into a numeric features-by-cells matrix and
#' row-bound onto the RNA count/data matrix. A new `RNA_BCR` assay is created
#' and the full Seurat pipeline (normalize, scale, PCA, neighbors, UMAP) is run
#' from scratch on the combined data. The main limitation is scale mismatch:
#' log-normalized RNA values and BCR metadata live in different ranges, so PCA
#' may be dominated by whichever modality has higher total variance even after
#' scaling.
#'
#' **`input_type = "embeddings"`, `pca_stage = "reduced"`** \cr
#' Same as the features variant above but uses BCR embeddings instead of
#' metadata columns. BCR embedding dimensions are appended to transposed GEX PCA
#' embeddings and a joint PCA is run on the combined matrix.
#'
#' **`input_type = "features"`, `pca_stage = "reduced"`** \cr
#' A middle ground between `"raw"` and `"embed"`. Instead of row-binding BCR
#' features onto the raw RNA matrix, they are appended to the transposed GEX PCA
#' embeddings (`rpca`, subset to `num_dims` PCs). The combined matrix
#' (n_gex_pcs + n_bcr_features rows) is stored as a new `RNA_BCR` assay,
#' scaled, and a joint PCA is run. This avoids the scale mismatch of `"raw"`
#' (GEX PCs and BCR metadata are in more comparable ranges) while still
#' performing a new projection that can mix the two modalities. Normalization is
#' skipped since GEX PCs are already processed. `filter_genes` does not apply.
#'
#' **`input_type = "embeddings"`, `pca_stage = "embed"`** \cr
#' Column-binds the existing `rpca` and `bpca` reductions directly. Requires a
#' merged object from [merge_gex_bcr()] with both reductions already computed.
#' This is the most efficient path when BCR embeddings are already available.
#' As with the features "embed" path, `num_dims` controls how many PCs are taken
#' from each reduction before joining.
#'
#' **`input_type = "features"`, `pca_stage = "embed"`** \cr
#' BCR metadata features are first embedded into their own PCA space: a
#' `BCR_meta` assay is created, scaled, and PCA is run to produce `bcr_meta.pca`
#' (capped at `nrow(bcr_features) - 1` PCs). The resulting BCR PCA embeddings
#' are then column-bound with the existing `rpca` embeddings to form a joint PCA
#' space. This is more principled than `"raw"` for metadata features because
#' both modalities are in comparable PCA spaces before being joined, and each
#' modality's internal variance structure is preserved. Use `num_dims` to
#' control how many PCs are taken from each side.
#'
#' For the `"raw"` and `"reduced"` paths, variable features are always set by
#' appending BCR features onto the existing RNA variable features rather than
#' re-running [FindVariableFeatures()]. Re-running would drop BCR features from
#' the selection and trigger Seurat warnings about underscores in feature names
#' and missing count layers.
#'
#' The `normalize`, `num_features`, `num_pcs`, `num_dims`, and `k_param`
#' arguments are passed to [seurat_pipeline()] for the `"raw"` and `"reduced"`
#' paths. For `"reduced"` and `"embed"`, `num_dims` additionally controls how
#' many GEX PCs are taken from `rpca` before concatenation.
#'
#' @param seurat_obj A Seurat object with RNA assay and BCR metadata.
#' @param pca_stage One of `"raw"`, `"reduced"`, or `"embed"`. `"raw"`
#'   concatenates at the count/feature level and runs a joint PCA from scratch.
#'   `"reduced"` appends BCR features onto transposed GEX PCA embeddings and
#'   runs a joint PCA. `"embed"` column-binds existing GEX and BCR PCA spaces
#'   directly. See Details for all combinations.
#' @param gex_reduction Name of the GEX PCA reduction to use as the GEX
#'   component for `pca_stage = "reduced"` and `pca_stage = "embed"`. Defaults
#'   to `"rpca"`. Use `"integrated"` to use a batch-corrected reduction (e.g.
#'   from Harmony).
#' @param input_type `"features"` to use processed BCR metadata columns;
#'   `"embeddings"` to use a pre-computed embedding matrix or the BCR assay
#'   from a merged object.
#' @param cols_to_include Character vector of BCR metadata column names to use
#'   as features (e.g. `c("mu_freq", "isotype")`). Required for
#'   `input_type = "features"`.
#' @param embeddings A features-by-cells matrix of BCR embeddings. Required for
#'   `input_type = "embeddings"` when no merged object is provided.
#' @param filter_genes If specified, filter out genes from this category
#'   (e.g. `"IG"` and/or `"TR"`). Only applies for `pca_stage = "raw"`.
#' @param ensembl_version Ensembl version for gene annotations (e.g.
#'   `"GRCh38.104"`). If `NULL`, auto-detected from
#'   `Misc(seurat_obj, "ensembl_version")` when available.
#' @param cache_file Passed to [get_airr_genes()]. Path to a cached RDS result
#'   to use instead of querying Ensembl live.
#' @param normalize If `TRUE`, normalize the `RNA_BCR` assay before scaling.
#'   Set to `FALSE` if the data has already been normalized. Ignored for
#'   `pca_stage = "reduced"` (GEX PCs do not need normalization).
#' @param num_features Number of variable features for the `"Before"` paths.
#' @param num_pcs Number of principal components to compute.
#' @param num_dims Number of PCA dimensions to use for neighbor finding. For
#'   `"After"` and `"Before_PCs"`, also controls how many GEX PCs are taken
#'   from `rpca` before concatenation.
#' @param k_param Number of nearest neighbors.
#' @param verbose Whether to show output from Seurat functions.
#'
#' @returns A Seurat object with:
#'   - New `RNA_BCR` assay (for `pca_stage = "Before"`)
#'   - PCA reduction (`rna_bcr.pca`)
#'   - UMAP reduction (`rna_bcr.umap`)
#'   - Neighbor graphs computed on the combined data
#' @export
concatenate_gex_bcr <- function(seurat_obj,
                                pca_stage = c("raw", "reduced", "embed"),
                                input_type = c("features", "embeddings"),
                                cols_to_include, embeddings = NULL,
                                gex_reduction = "rpca",
                                filter_genes, ensembl_version = NULL,
                                cache_file = NULL,
                                normalize = TRUE, num_features = 2000,
                                num_pcs = 50, num_dims = 20, k_param = 20,
                                verbose = TRUE) {
  # TODO: add an option to do weighting to influence the effect of the BCRs (post-normalization)

  pca_stage  <- match.arg(pca_stage)
  input_type <- match.arg(input_type)

  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("{.arg seurat_obj} must be a Seurat object.")
  }
  if (!"cell_id" %in% colnames(seurat_obj[[]])) {
    cli::cli_abort("Cell ID column not found in metadata.")
  }
  # auto-detect ensembl version from object misc slot when filtering
  if (!missing(filter_genes) && is.null(ensembl_version)) {
    ensembl_version <- Misc(seurat_obj, "ensembl_version")
    if (!is.null(ensembl_version) && verbose) {
      cli::cli_inform(c("i" = "Using Ensembl version {ensembl_version} from object."))
    }
  }

  # `missing()` only works on formal args of the function it's called from, so
  # validate cols_to_include here rather than inside each branch
  if (input_type == "features" && missing(cols_to_include)) {
    cli::cli_abort(c(
      '{.arg cols_to_include} is required for {.code input_type = "features"}.',
      "i" = 'e.g. {.code cols_to_include = c("mu_freq", "isotype")}.'
    ))
  }

  # detect merged object (has BCR assay + bpca from merge_gex_bcr)
  is_merged <- "BCR" %in% names(seurat_obj@assays) &&
    "bpca" %in% names(seurat_obj@reductions)

  # for non-merged embeddings, subset to shared cells before any branch logic
  if (input_type == "embeddings" && !is_merged && !is.null(embeddings)) {
    shared <- intersect(colnames(embeddings), seurat_obj$cell_id)
    if (length(shared) < ncol(seurat_obj)) {
      if (verbose) {
        cli::cli_inform(c("i" = "Embeddings cover {length(shared)} of {ncol(seurat_obj)} \\
                                 cells; subsetting to shared set."))
      }
      seurat_obj <- subset(seurat_obj, cells = shared)
    }
  }

  if (pca_stage == "raw") {
    # warn loudly if normalizing embeddings — log1p of PLM values is meaningless
    if (input_type == "embeddings" && normalize) {
      cli::cli_warn(c(
        "!" = "{.code normalize = TRUE} with {.code input_type = 'embeddings'} \\
               will log-normalize embedding values alongside RNA counts.",
        "i" = "For pretrained embeddings, prefer {.code normalize = FALSE}."
      ))
    }

    # --- resolve BCR data (features x cells matrix) ---
    if (input_type == "embeddings") {
      if (is_merged) {
        bcr_features <- GetAssayData(seurat_obj, assay = "BCR", layer = "data")
        if (verbose) {
          cli::cli_inform(c("i" = "Using BCR assay ({nrow(bcr_features)} dimensions) \\
                                   from merged object."))
        }
      } else if (!is.null(embeddings)) {
        bcr_features <- embeddings[, Cells(seurat_obj), drop = FALSE]
        if (verbose) {
          cli::cli_inform(c("i" = "Using provided embeddings matrix \\
                                   ({nrow(bcr_features)} dimensions, {ncol(bcr_features)} cells)."))
        }
      } else {
        cli::cli_abort(c(
          'Must provide a merged object or an {.arg embeddings} matrix \\
           for {.code input_type = "embeddings"}.',
          "i" = "Run {.fn merge_gex_bcr} first, or pass an embeddings matrix."
        ))
      }
    } else {
      bcr_features <- seurat_obj[[]] %>% select(all_of(cols_to_include))
      bcr_features <- process_bcr_features(bcr_features, verbose = verbose)
    }

    # Seurat 5 silently rewrites "_" to "-" in feature names at assay creation,
    # de-syncing the assay names from the VariableFeatures list and causing
    # ScaleData to find nothing to scale. Using "." avoids the rewrite.
    rownames(bcr_features) <- gsub("_", ".", rownames(bcr_features))

    # write to counts when normalizing (seurat_pipeline will call NormalizeData
    # which populates data), or directly to data when already normalized
    if (normalize) {
      seurat_obj[["RNA_BCR"]] <- CreateAssay5Object(
        counts = rbind(GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
                       bcr_features))
    } else {
      seurat_obj[["RNA_BCR"]] <- CreateAssay5Object(
        data = rbind(GetAssayData(seurat_obj, assay = "RNA", layer = "data"),
                     bcr_features))
    }
    DefaultAssay(seurat_obj) <- "RNA_BCR"

    # always append BCR features onto existing RNA variable features rather than
    # re-running FindVariableFeatures, which would drop BCR features from the
    # selection and trigger Seurat warnings about underscores in feature names
    # and about the count layer being missing when normalize = FALSE
    collisions <- intersect(VariableFeatures(seurat_obj, assay = "RNA"),
                            rownames(bcr_features))
    if (length(collisions) > 0) {
      cli::cli_warn(c(
        "!" = "{length(collisions)} BCR feature name{?s} collide with existing \\
               RNA variable features: {.val {collisions}}.",
        "i" = "Consider renaming BCR features to avoid ambiguity."
      ))
    }
    n_rna_var <- length(VariableFeatures(seurat_obj, assay = "RNA"))
    VariableFeatures(seurat_obj) <-
      unique(c(VariableFeatures(seurat_obj, assay = "RNA"), rownames(bcr_features)))
    if (verbose) {
      n_total <- length(VariableFeatures(seurat_obj))
      cli::cli_inform(c("i" = "Variable features: {n_rna_var} RNA + \\
                               {nrow(bcr_features)} BCR = {n_total} total."))
    }

    if (!missing(filter_genes)) {
      seurat_obj <- filter_variable_features(seurat_obj, filter_genes,
                                             ensembl_version = ensembl_version,
                                             cache_file = cache_file,
                                             bcr_features = bcr_features)
    }

    seurat_obj <- seurat_pipeline(seurat_obj, assay = "RNA_BCR",
                                  pca_name = "rna_bcr.pca",
                                  normalize = normalize,
                                  num_features = num_features,
                                  num_pcs = num_pcs, num_dims = num_dims,
                                  k_param = k_param,
                                  find_var_features = FALSE,
                                  verbose = verbose)
  } else if (pca_stage == "reduced") {
    if (!gex_reduction %in% names(seurat_obj@reductions)) {
      cli::cli_abort(c(
        "Reduction {.code {gex_reduction}} not found.",
        "i" = "Run {.fn seurat_pipeline} first to compute the GEX PCA."
      ))
    }

    # resolve BCR data — same logic as "raw"
    if (input_type == "embeddings") {
      if (is_merged) {
        bcr_features <- GetAssayData(seurat_obj, assay = "BCR", layer = "data")
        if (verbose) {
          cli::cli_inform(c("i" = "Using BCR assay ({nrow(bcr_features)} dimensions) \\
                                   from merged object."))
        }
      } else if (!is.null(embeddings)) {
        bcr_features <- embeddings[, Cells(seurat_obj), drop = FALSE]
        if (verbose) {
          cli::cli_inform(c("i" = "Using provided embeddings matrix \\
                                   ({nrow(bcr_features)} dimensions, {ncol(bcr_features)} cells)."))
        }
      } else {
        cli::cli_abort(c(
          'Must provide a merged object or an {.arg embeddings} matrix \\
           for {.code input_type = "embeddings"}.',
          "i" = "Run {.fn merge_gex_bcr} first, or pass an embeddings matrix."
        ))
      }
    } else {
      bcr_features <- seurat_obj[[]] %>% select(all_of(cols_to_include))
      bcr_features <- process_bcr_features(bcr_features, verbose = verbose)
    }

    # subset GEX PCA to num_dims PCs and transpose to features-by-cells
    n_gex_dims <- min(num_dims, ncol(Embeddings(seurat_obj, gex_reduction)))
    gex_pca_mat <- t(Embeddings(seurat_obj, gex_reduction)[, seq_len(n_gex_dims), drop = FALSE])

    # same "_" → "." rename as "raw" branch — avoids Seurat 5's silent rewrite
    rownames(gex_pca_mat)  <- gsub("_", ".", rownames(gex_pca_mat))
    rownames(bcr_features) <- gsub("_", ".", rownames(bcr_features))

    combined_mat <- rbind(gex_pca_mat, bcr_features)
    if (verbose) {
      cli::cli_inform(c("i" = "Combining {n_gex_dims} GEX PCs + \\
                               {nrow(bcr_features)} BCR features = \\
                               {nrow(combined_mat)} rows."))
    }

    # write to counts then copy to data — ScaleData reads from data, not counts,
    # but CreateAssay5Object requires at least one of counts/data to be provided
    seurat_obj[["RNA_BCR"]] <- CreateAssay5Object(counts = combined_mat)
    LayerData(seurat_obj, assay = "RNA_BCR", layer = "data") <-
      LayerData(seurat_obj, assay = "RNA_BCR", layer = "counts")
    DefaultAssay(seurat_obj) <- "RNA_BCR"
    VariableFeatures(seurat_obj) <- rownames(combined_mat)

    seurat_obj <- seurat_pipeline(seurat_obj, assay = "RNA_BCR",
                                  pca_name = "rna_bcr.pca",
                                  normalize = FALSE,
                                  num_features = num_features,
                                  num_pcs = num_pcs, num_dims = num_dims,
                                  k_param = k_param,
                                  find_var_features = FALSE,
                                  verbose = verbose)
  } else {
    # pca_stage == "embed": column-bind GEX and BCR PCA spaces
    if (!gex_reduction %in% names(seurat_obj@reductions)) {
      cli::cli_abort(c(
        "Reduction {.code {gex_reduction}} not found.",
        "i" = 'Run {.fn RunPCA} with {.code reduction.name = "{gex_reduction}"} first.'
      ))
    }

    if (input_type == "features") {
      # embed BCR metadata features into their own PCA space first, then cbind
      bcr_features <- seurat_obj[[]] %>% select(all_of(cols_to_include))
      bcr_features <- process_bcr_features(bcr_features, verbose = verbose)

      # cap PCs at number of BCR features (can't exceed rank of the matrix)
      n_bcr_pcs <- min(nrow(bcr_features) - 1L, num_dims)
      if (verbose) {
        cli::cli_inform(c("i" = "Computing BCR PCA from {nrow(bcr_features)} metadata \\
                                 features ({n_bcr_pcs} PCs)."))
      }
      # BCR metadata is already on a sane scale; write to data so ScaleData
      # doesn't need a counts layer
      # TODO: call this RNA_BCR for consistency??
      colnames(bcr_features) <- Cells(seurat_obj)
      seurat_obj[["BCR_meta"]] <- CreateAssay5Object(counts = bcr_features)
      LayerData(seurat_obj, assay = "BCR_meta", layer = "data") <-
        LayerData(seurat_obj, assay = "BCR_meta", layer = "counts")
      VariableFeatures(seurat_obj, assay = "BCR_meta") <- rownames(bcr_features)
      seurat_obj <- ScaleData(seurat_obj, assay = "BCR_meta", verbose = verbose)
      scale_bcr <- GetAssayData(seurat_obj, assay = "BCR_meta", layer = "scale.data")
      use_approx <- n_bcr_pcs < min(nrow(scale_bcr), ncol(scale_bcr)) / 2
      seurat_obj <- RunPCA(seurat_obj, assay = "BCR_meta", npcs = n_bcr_pcs,
                           reduction.name = "bcr_meta.pca",
                           reduction.key = "bcrmetapca_",
                           approx = use_approx, verbose = verbose)
      bcr_reduc <- "bcr_meta.pca"
    } else {
      # embeddings: use existing bpca from merge_gex_bcr
      if (!is_merged) {
        cli::cli_abort(c(
          '{.code input_type = "embeddings"} with {.code pca_stage = "embed"} \\
           requires a merged object with a BCR assay and {.code bpca} reduction.',
          "i" = "Run {.fn merge_gex_bcr} first."
        ))
      }
      if (!"bpca" %in% names(seurat_obj@reductions)) {
        cli::cli_abort(c(
          "Reduction {.code bpca} not found.",
          "i" = 'Run {.fn RunPCA} with {.code reduction.name = "bpca"} first.'
        ))
      }
      bcr_reduc <- "bpca"
    }

    n_gex_dims <- min(num_dims, ncol(Embeddings(seurat_obj, gex_reduction)))
    n_bcr_dims <- min(num_dims, ncol(Embeddings(seurat_obj, bcr_reduc)))
    total_dims  <- n_gex_dims + n_bcr_dims

    combined_pca <-
      cbind(Embeddings(seurat_obj, gex_reduction)[, seq_len(n_gex_dims), drop = FALSE],
            Embeddings(seurat_obj, bcr_reduc)[, seq_len(n_bcr_dims), drop = FALSE])
    colnames(combined_pca) <- str_c("rnabcrpca_", seq_len(total_dims))

    seurat_obj[["rna_bcr.pca"]] <-
      CreateDimReducObject(embeddings = combined_pca, key = "rnabcrpca_",
                           assay = DefaultAssay(seurat_obj))
    if (verbose) {
      cli::cli_inform(c("i" = "Combined PCA: {n_gex_dims} GEX + \\
                               {n_bcr_dims} BCR = {total_dims} dims."))
    }

    seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rna_bcr.pca",
                                k.param = k_param, dims = 1:total_dims,
                                verbose = verbose,
                                graph.name = c("RNA_BCR_nn", "RNA_BCR_snn"))
    seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "rna_bcr.pca",
                                k.param = k_param, dims = 1:total_dims,
                                return.neighbor = TRUE,
                                verbose = verbose, graph.name = "RNA_BCR.nn")
    seurat_obj <- RunUMAP(object = seurat_obj, nn.name = "RNA_BCR.nn",
                          n.neighbors = k_param, reduction.name = "rna_bcr.umap",
                          reduction.key = "rnabcrUMAP_", verbose = verbose)
  }

  Misc(seurat_obj, slot = "concat_pca_stage")  <- pca_stage
  Misc(seurat_obj, slot = "concat_input_type") <- input_type

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Concatenated GEX & BCR pipeline complete.",
      "i" = "Use {.fn Seurat::DimPlot} with {.arg reduction = 'rna_bcr.umap'} or \\
             {.fn athanor::plot_dimplot} with {.arg reduc = 'rna_bcr.umap'} to visualize."
    ))
  }

  return(seurat_obj)
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
#' @param seurat_obj A Seurat object containing GEX data (at the least). If this
#'   is already a merged object produced by [merge_gex_bcr()] (i.e., it contains
#'   a BCR assay and a `bpca` reduction), `embeddings` and `embedding_type` may
#'   be omitted and the existing BCR infrastructure will be reused.
#' @param embeddings Matrix of BCR embeddings (genes by cells). Optional if
#'   `seurat_obj` already contains a BCR assay and `bpca` reduction.
#' @param embedding_type The embeddings method. Optional for merged objects;
#'   falls back to `seurat_obj@misc$embedding_type` if not provided.
#' @param pc_gex The number of PCs for the GEX assay.
#' @param pc_bcr The number of PCs for the BCR assay.
#' @param k The number of neighbors to use for each modality.
#' @param cluster Whether or not to perform clustering.
#' @param cluster_res Named list of clustering resolutions for GEX, BCR, and WNN.
#' @param modality_weights Named vector of modality weights. If NULL, Seurat will calculate automatically.
#' @param verbose Logical indicating whether or not to print messages.
#'
#' @returns A Seurat object with WNN run.
#' @export
run_wnn <- function(seurat_obj, embeddings, embedding_type, pc_gex = 20,
                    pc_bcr = 20, k_param = 20,
                    gex_reduction = "rpca", bcr_reduction = "bpca",
                    cluster = FALSE,
                    cluster_res = list("GEX" = 1, "BCR" = 1, "WNN" = 1),
                    modality_weights = NULL, verbose = TRUE) {
  # TODO: update this to be able to run on other omics e.g. GEX & ADT
  # TODO: give the option to use an integrated GEX assay

  # input validation
  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("seurat_obj must be a Seurat object.")
  }
  if (!"cell_id" %in% colnames(seurat_obj[[]])) {
    cli::cli_abort("Cell ID column not found in metadata.")
  }
  if (missing(pc_gex)) {
    pc_gex <- if (gex_reduction %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[[gex_reduction]])
    } else if ("pca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["pca"]])
    } else {
      20
    }
    cli::cli_inform(c("i" = "Using pc_gex = {pc_gex} from existing reductions."))
  }
  if (missing(pc_bcr)) {
    pc_bcr <- if (bcr_reduction %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[[bcr_reduction]])
    } else {
      20
    }
    cli::cli_inform(c("i" = "Using pc_bcr = {pc_bcr} from existing reductions."))
  }
  if (missing(k_param)) {
    k_param <- infer_k_param(seurat_obj, verbose = verbose)
  }

  # detect if the object is already merged (has BCR assay + bcr_reduction from merge_gex_bcr)
  is_merged <- "BCR" %in% names(seurat_obj@assays) &&
    bcr_reduction %in% names(seurat_obj@reductions)

  if (is_merged) {
    cli::cli_inform(c("i" = "Merged object detected: using existing GEX and BCR assays and reductions."))
    if (missing(embedding_type)) {
      embedding_type <- seurat_obj@misc$embedding_type # should exist
    }
  } else {
    if (missing(embeddings)) {
      cli::cli_abort(c(
        "Must provide {.arg embeddings} or pass a merged object.",
        "i" = "Run {.fn merge_gex_bcr} first, or provide an embeddings matrix."
      ))
    }
    if (ncol(embeddings) == 0) {
      cli::cli_abort("No cells found in the provided embeddings matrix.")
    }
  }

  # remove any clustering that might exist and reset the factor levels
  # assumes a UMAP might have been run for the RNA data (instead of tSNE)
  cols_to_remove <- c()
  if (any(grepl("^snn_res", names(seurat_obj[[]])))) {
    cols_to_remove <-
      c(cols_to_remove,
        grep("^snn_res", names(seurat_obj[[]]), value = TRUE))
  }
  if ("seurat_clusters" %in% names(seurat_obj[[]])) {
    cols_to_remove <- c(cols_to_remove, "seurat_clusters")
  }
  seurat_obj@meta.data <-
    seurat_obj[[]] %>%
    {if (length(cols_to_remove) > 0)
      select(., -all_of(cols_to_remove)) else .} %>%
    droplevels()

  if (!is_merged) {
    # TODO: double check this part
    # TODO: use merge_gex_bcr??

    # build BCR object (subset to cells present in the GEX object)
    bcr_obj <- bcr_embeddings_pipeline(
      embeddings[, intersect(colnames(embeddings), seurat_obj$cell_id)],
      embedding_type = embedding_type,
      num_pcs = pc_bcr, num_dims = pc_bcr, k_param = k_param, verbose = verbose)

    # add BCR assay; Seurat v5 does not propagate nCount/nFeature automatically
    suppressWarnings(seurat_obj[["BCR"]] <- bcr_obj[["BCR"]])
    seurat_obj$nCount_BCR <- bcr_obj[[]][Cells(seurat_obj), "nCount_BCR"]
    seurat_obj$nFeature_BCR <- bcr_obj[[]][Cells(seurat_obj), "nFeature_BCR"]
    if ("nFeature_ADT" %in% names(seurat_obj[[]])) {
      seurat_obj@meta.data <-
        seurat_obj[[]] %>%
        relocate(nCount_BCR, nFeature_BCR, .after = nFeature_ADT)
    }

    if (ncol(seurat_obj@assays$RNA) != ncol(seurat_obj@assays$BCR)) {
      cli::cli_abort("The number of cells in the RNA and BCR assays do not match.")
    }

    # transfer reductions and graphs
    seurat_obj[[bcr_reduction]] <- bcr_obj[["bpca"]]
    seurat_obj[["bcr.umap"]] <- bcr_obj[["bcr.umap"]]
    seurat_obj@neighbors[["BCR.nn"]] <- bcr_obj@neighbors[["BCR.nn"]]
    seurat_obj@graphs[["BCR_nn"]] <- bcr_obj@graphs[["BCR_nn"]]
    seurat_obj@graphs[["BCR_snn"]] <- bcr_obj@graphs[["BCR_snn"]]

    # add run info
    Misc(seurat_obj, slot = "embedding_type") <- embedding_type
    Misc(seurat_obj, slot = "embedding_dims") <- nrow(seurat_obj@assays[["BCR"]])
  }

  # find multimodal neighbors
  if (is.null(modality_weights)) {
    seurat_obj <-
      Seurat::FindMultiModalNeighbors(object = seurat_obj,
                                      reduction.list = list(gex_reduction,
                                                            bcr_reduction),
                                      dims.list = list(1:pc_gex, 1:pc_bcr),
                                      k.nn = k_param,
                                      knn.graph.name = "w_nn",
                                      snn.graph.name = "w_snn",
                                      weighted.nn.name = "w.nn",
                                      modality.weight.name =
                                        str_c(c("RNA", "BCR"), ".weight"),
                                      return.intermediate = TRUE,
                                      modality.weight = NULL,
                                      verbose = verbose)

    seurat_obj@meta.data <-
      seurat_obj[[]] %>%
      mutate(weight_assay =
               case_when(seurat_obj[[]][["RNA.weight"]] > 0.5 ~ "RNA",
                         seurat_obj[[]][["RNA.weight"]] < 0.5 ~ "BCR",
                         is.na(seurat_obj[[]][["RNA.weight"]]) ~ NA,
                         .default = "Tie"))

    wnn_nn_name    <- "w.nn"
    rna_weight_col <- "RNA.weight"
    wnn_umap_name  <- "wnn.umap"
    wnn_umap_key   <- "wnnUMAP_"
  } else {
    # custom weights: "equal" or a length-2 numeric vector c(gex_wt, bcr_wt)
    if (identical(modality_weights, "equal")) {
      gex_wt    <- 0.5
      bcr_wt    <- 0.5
      mw_suffix <- "equal"
    } else if (is.numeric(modality_weights) && length(modality_weights) == 2) {
      if (abs(sum(modality_weights) - 1) > 1e-6) {
        cli::cli_abort("{.arg modality_weights} values must sum to 1.")
      }
      gex_wt    <- modality_weights[[1]]
      bcr_wt    <- modality_weights[[2]]
      mw_suffix <- paste0(gex_wt, "_", bcr_wt)
    } else {
      cli::cli_abort(c(
        "Invalid {.arg modality_weights}.",
        "i" = "Provide {.code NULL} for automatic weighting, \\
               {.code \"equal\"} for 50/50, or a length-2 numeric vector \\
               summing to 1 (e.g. {.code c(0.7, 0.3)})."
      ))
    }

    if (is.null(seurat_obj@misc$modality.weight)) {
      cli::cli_abort(c(
        "No existing WNN run found in {.code @misc$modality.weight}.",
        "i" = "Run {.fn run_wnn} with {.code modality_weights = NULL} first \\
               to generate the initial WNN, then re-run with custom weights."
      ))
    }

    mw_mod      <- seurat_obj@misc$modality.weight
    mw_mod_names <- Cells(seurat_obj)
    mw_mod@modality.weight.list[[gex_reduction]] <-
      setNames(rep(gex_wt, ncol(seurat_obj)), mw_mod_names)
    mw_mod@modality.weight.list[[bcr_reduction]] <-
      setNames(rep(bcr_wt, ncol(seurat_obj)), mw_mod_names)

    wnn_nn_name    <- paste0("w_", mw_suffix, ".nn")
    rna_weight_col <- paste0("RNA_", mw_suffix, ".weight")
    wnn_umap_name  <- paste0("wnn_", mw_suffix, ".umap")
    wnn_umap_key   <- paste0("wnn", mw_suffix, "UMAP_")

    seurat_obj <-
      Seurat::FindMultiModalNeighbors(object = seurat_obj,
                                      reduction.list = list(gex_reduction,
                                                            bcr_reduction),
                                      dims.list = list(1:pc_gex, 1:pc_bcr),
                                      k.nn = k_param,
                                      knn.graph.name = paste0("w_", mw_suffix, "_nn"),
                                      snn.graph.name = paste0("w_", mw_suffix, "_snn"),
                                      weighted.nn.name = wnn_nn_name,
                                      modality.weight.name =
                                        str_c(c("RNA", "BCR"), "_", mw_suffix, ".weight"),
                                      return.intermediate = TRUE,
                                      modality.weight = mw_mod,
                                      verbose = verbose)

    seurat_obj@meta.data <-
      seurat_obj[[]] %>%
      mutate(!!paste0("weight_assay_", mw_suffix) :=
               case_when(seurat_obj[[]][[rna_weight_col]] > 0.5 ~ "RNA",
                         seurat_obj[[]][[rna_weight_col]] < 0.5 ~ "BCR",
                         is.na(seurat_obj[[]][[rna_weight_col]]) ~ NA,
                         .default = "Tie"))
  }

  # check for NA values in the neighbor index
  na_nn <- is.na(seurat_obj[[rna_weight_col]])
  if (any(na_nn)) {
    bad_cells <- seurat_obj[[]]$cell_id[na_nn]
    bcr_mat <- GetAssayData(seurat_obj, assay = "BCR", layer = "data")
    bcr_mat <- as.matrix(bcr_mat[, bad_cells])
    # add 1 because the first column is not included in the duplicate count
    dup_count <- sum(duplicated(t(bcr_mat))) + 1

    Misc(seurat_obj, slot = "embedding_dups") <- bad_cells

    cli::cli_inform(c(
      "x" = "{sum(na_nn)} cell{?s} (out of {ncol(seurat_obj)}) have no valid \\
             neighbors in the WNN space.",
      "i" = "This is most likely due to identical {bcr_reduction} embeddings. \\
             Out of the {sum(na_nn)} affected cell{?s}, {dup_count} \\
             {?has/have} duplicate embeddings.",
      ">" = "UMAP will not be run for this WNN reduction.",
      ">" = "Consider filtering duplicate {bcr_reduction} embeddings before \\
             calling {.fn run_wnn}: \\
             {.code seurat_obj <- subset(seurat_obj, !cell_id %in% embedding_dups)}"
    ))
  } else {
    seurat_obj <-
      Seurat::RunUMAP(object = seurat_obj,
                      nn.name = wnn_nn_name,
                      n.neighbors = k_param,
                      reduction.name = wnn_umap_name,
                      reduction.key = wnn_umap_key,
                      verbose = verbose)
  }

  # the Leiden algorithm (4) has been shown to be better than Louvain (1), but
  # in order to use it you have to install the 'leidenbase' package
  if (cluster) {
    algo <- 1

    # cluster the BCR assay
    seurat_obj <- FindClusters(object = seurat_obj, graph.name = "BCR_snn",
                               resolution = cluster_res[["BCR"]],
                               algorithm = algo, verbose = verbose)
    # cluster the GEX assay
    # TODO: do this in seurat_pipeline??
    seurat_obj <- FindClusters(object = seurat_obj, graph.name = "RNA_snn",
                               resolution = cluster_res[["GEX"]],
                               algorithm = algo, verbose = verbose)
    # cluster the WNN assay
    seurat_obj <- FindClusters(object = seurat_obj, graph.name = "w_snn",
                               resolution = cluster_res[["WNN"]],
                               algorithm = algo, verbose = verbose)

    # set the cluster identities and fix the order (RNA and BCR are fine)
    meta_res_wnn <- paste0("w_snn_res.", cluster_res[["WNN"]])
    seurat_obj[[]][[meta_res_wnn]] <- fct_inseq(seurat_obj[[]][[meta_res_wnn]])
    Idents(seurat_obj) <- meta_res_wnn
  }

  if (verbose) {
    if (!any(na_nn)) {
      cli::cli_inform(c(
        "v" = "WNN neighbors calculated and UMAP run.",
        "i" = "Use {.fn Seurat::DimPlot} with {.code reduction = '{wnn_umap_name}'} \\
               or {.fn athanor::plot_dimplot} with {.code reduc = '{wnn_umap_name}'}."))
    }

    if (is.null(modality_weights)) {
      cli::cli_inform(c(
        "v" = "Modality weights were automatically calculated.",
        "i" = "See {.fn Seurat::FindMultiModalNeighbors} for details."))
    } else {
      cli::cli_inform(c(
        "v" = "Custom modality weights applied ({gex_reduction}: {gex_wt}, \\
               {bcr_reduction}: {bcr_wt})."))
    }
  }

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
#' @param other_type The second assay.
#'
#' @returns A text message.
#' @export
extract_wnn_vars <- function(seurat_obj, gex_pca = "rpca",
                             other_pca = "bpca", other_type = "BCR") {
  # TODO: update this to work on a concatenated object too
  # TODO: print info about reduction

  # assay details
  assay_details <- character(0)
  bcr_embeddings_name <- grep("^BCR", names(seurat_obj@assays), value = TRUE)

  if ("RNA" %in% names(seurat_obj@assays)) {
    assay_details <- c(assay_details,
                       paste(nrow(seurat_obj@assays$RNA), "genes"))
  }
  if ("ADT" %in% names(seurat_obj@assays)) {
    assay_details <- c(assay_details,
                       paste(nrow(seurat_obj@assays$ADT),
                             "cell surface protein markers"))
  }
  if (length(bcr_embeddings_name) > 0) {
    assay_details <- c(assay_details,
                       paste(nrow(seurat_obj@assays[[bcr_embeddings_name]]),
                             "embedding dimensions"))
  }

  cli::cli_inform(c(
    "i" = "This object has {ncol(seurat_obj)} cells with \\
           {str_c(assay_details, collapse = ', ')}."
  ))

  # WNN reduction info
  # TODO: fix the count
  if ("w.nn" %in% names(seurat_obj@neighbors) &&
      all(c(gex_pca, other_pca) %in% names(seurat_obj@reductions))) {

    # objs_wnn[["immune2vec"]]@commands[["Seurat..FindMultiModalNeighbors"]]@params[["dims.list"]]
    cli::cli_inform(c(
      "i" = "WNN was run with {ncol(seurat_obj@reductions[[gex_pca]])} GEX PCs \\
             and {ncol(seurat_obj@reductions[[other_pca]])} {other_type} PCs."
    ))
  }

  # clustering info (largest resolution only)
  if ("seurat_clusters" %in% colnames(seurat_obj[[]])) {
    nclusters_rna <- n_distinct(
      select(seurat_obj[[]], starts_with("RNA_snn")) %>% select(last_col()))
    nclusters_bcr <- n_distinct(
      select(seurat_obj[[]], starts_with("BCR_snn")) %>% select(last_col()))
    nclusters_wnn <- n_distinct(
      select(seurat_obj[[]], starts_with("w_snn")) %>% select(last_col()))

    cli::cli_inform(c(
      "i" = "Clusters (largest resolution): {nclusters_rna} RNA, \\
             {nclusters_bcr} BCR, {nclusters_wnn} WNN."
    ))
  }

  invisible(seurat_obj)
}
