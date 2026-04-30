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
#' @param verbose Whether to print integration summary information.
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
    stop("Please check the cell id format in the Seurat object or AIRR table.")
  }

  if (any(!new_cols %in% colnames(combined_airr))) {
    stop("Make sure that you are only adding columns which exist in the AIRR table.")
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
  cells_main <- (filter(combined_airr, locus == "IGH"))$cell_id
  cells_other <- (filter(combined_airr, locus != "IGH"))$cell_id
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

    # TODO: print out that this wasn't added if applicable
    if ("annotated_clusters_simpler" %in% colnames(seurat_obj[[]])) {
             # assumes that the plasma cells don't need re-assigning
      seurat_obj@meta.data <-
        seurat_obj[[]] %>%
        mutate(annotated_clusters_gex_bcr =
                 ifelse(annotated_clusters_simpler == "Naive B cells" &
                          (mu_freq > 0.01 | isotype %in% c("IgA", "IgE", "IgG")),
                        "Memory B cells", annotated_clusters_simpler))
    }
  }

  # make sure that the levels are okay
  seurat_obj@meta.data <- seurat_obj@meta.data %>% droplevels()

  # print an verbose of the integration
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
#' BCR cell names are renamed via [Seurat::RenameCells()] to match before the
#' assay is transferred.
#'
#' With `join = "inner"` (default), both objects are subset to their shared
#' cells; the result is suitable for WNN or concatenation workflows.
#'
#' With `join = "left"`, all GEX cells are retained. Cells without BCR data
#' receive a zero-filled embedding row in the BCR assay and `Has_BCR = FALSE`.
#' Zero embeddings are meaningless and downstream analyses should filter on
#' `Has_BCR` before using the BCR assay.
#'
#' BCR metadata columns not already present in `gex_obj` are transferred via
#' [Seurat::AddMetaData()], with `NA` (or `FALSE` for `Has_BCR`) filled in for
#' non-BCR cells.
#'
#' @param gex_obj A Seurat object containing GEX (RNA) data with a `cell_id`
#'   metadata column.
#' @param bcr_obj A Seurat object containing a BCR assay with a `cell_id`
#'   metadata column (typically produced by [bcr_embeddings_pipeline()]).
#' @param join How to handle cells not shared between the two objects.
#'   `"inner"` (default) keeps only shared cells; `"left"` keeps all GEX cells
#'   and zero-fills the BCR assay for unmatched cells.
#' @param transfer_reductions Whether to copy BCR reductions (`bpca`,
#'   `bcr.umap`) and graphs (`BCR.nn`, `BCR_nn`, `BCR_snn`) from `bcr_obj`
#'   into the merged object. Only available when `join = "inner"`.
#' @param verbose Whether to print a summary of the merge.
#'
#' @returns A Seurat object with both RNA and BCR assays and BCR metadata columns.
#' @export
merge_gex_bcr <- function(gex_obj, bcr_obj, join = c("inner", "left"),
                           transfer_reductions = TRUE, verbose = TRUE) {
  join <- match.arg(join)

  if (!inherits(gex_obj, "Seurat")) stop("gex_obj must be a Seurat object.")
  if (!inherits(bcr_obj, "Seurat")) stop("bcr_obj must be a Seurat object.")
  if (!"cell_id" %in% colnames(gex_obj[[]])) {
    stop("gex_obj must have a cell_id metadata column.")
  }
  if (!"cell_id" %in% colnames(bcr_obj[[]])) {
    stop("bcr_obj must have a cell_id metadata column.")
  }
  if (!"BCR" %in% names(bcr_obj@assays)) {
    stop("bcr_obj must contain a BCR assay.")
  }

  shared_cell_ids <- intersect(gex_obj$cell_id, bcr_obj$cell_id)
  if (length(shared_cell_ids) == 0) {
    stop("No shared cell IDs found. Check cell_id formatting in both objects.")
  }

  n_gex <- ncol(gex_obj)
  n_bcr <- ncol(bcr_obj)

  # subset BCR object to shared cells
  bcr_shared_barcodes <- Cells(bcr_obj)[bcr_obj$cell_id %in% shared_cell_ids]
  bcr_obj_sub <- subset(bcr_obj, cells = bcr_shared_barcodes)

  # remap BCR cell names to GEX barcodes if they differ from cell_id
  gex_cellid_to_barcode <- setNames(rownames(gex_obj[[]]), gex_obj$cell_id)
  gex_barcodes_for_bcr <- unname(gex_cellid_to_barcode[bcr_obj_sub$cell_id])
  if (!identical(Cells(bcr_obj_sub), gex_barcodes_for_bcr)) {
    bcr_obj_sub <- RenameCells(bcr_obj_sub, new.names = gex_barcodes_for_bcr)
  }

  if (join == "inner") {
    gex_obj <- subset(gex_obj, cells = gex_barcodes_for_bcr)
    suppressWarnings(gex_obj[["BCR"]] <- bcr_obj_sub[["BCR"]])
  } else {
    # left join: build a full BCR count matrix, zero-filling non-BCR cells
    bcr_feature_names <- rownames(bcr_obj_sub[["BCR"]])
    n_features <- length(bcr_feature_names)
    all_gex_barcodes <- Cells(gex_obj)
    non_bcr_barcodes <- setdiff(all_gex_barcodes, gex_barcodes_for_bcr)

    bcr_counts_shared <- GetAssayData(bcr_obj_sub, assay = "BCR", layer = "counts")

    zero_mat <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0),
      dims = c(n_features, length(non_bcr_barcodes)),
      dimnames = list(bcr_feature_names, non_bcr_barcodes)
    )

    # reorder shared cells to match gex_obj order, then append non-BCR cells
    bcr_counts_full <- cbind(
      bcr_counts_shared[, intersect(all_gex_barcodes, colnames(bcr_counts_shared)),
                        drop = FALSE],
      zero_mat
    )[, all_gex_barcodes, drop = FALSE]

    suppressWarnings(gex_obj[["BCR"]] <- CreateAssayObject(counts = bcr_counts_full))

    if (transfer_reductions && length(non_bcr_barcodes) > 0) {
      cli::cli_warn(c(
        "!" = "transfer_reductions is not supported for join = \"left\" when \\
GEX cells without BCR data are present.",
        "i" = "Reductions computed on BCR-only cells cannot be applied to \\
{length(non_bcr_barcodes)} non-BCR cells."
      ))
      transfer_reductions <- FALSE
    }
  }

  # transfer nCount_BCR and nFeature_BCR, aligned to gex_obj cells
  bcr_counts_meta <-
    bcr_obj_sub[[]] %>%
    select(cell_id, nCount_BCR, nFeature_BCR)

  gex_counts_meta <- left_join(gex_obj[["cell_id"]], bcr_counts_meta, by = "cell_id")
  gex_obj$nCount_BCR <- gex_counts_meta$nCount_BCR
  gex_obj$nFeature_BCR <- gex_counts_meta$nFeature_BCR

  # relocate count columns after the last RNA/ADT count column
  last_count_col <- tail(grep("^nFeature_(RNA|ADT)$", names(gex_obj[[]]),
                               value = TRUE), 1)
  if (length(last_count_col) > 0) {
    gex_obj@meta.data <- gex_obj[[]] %>%
      relocate(nCount_BCR, nFeature_BCR, .after = all_of(last_count_col))
  }

  # transfer any BCR metadata columns not already in gex_obj;
  # always exclude Has_BCR here — we recompute it below from actual BCR presence
  skip_cols <- c("orig.ident", "nCount_BCR", "nFeature_BCR", "Has_BCR",
                 names(gex_obj[[]]))
  bcr_meta_extra <-
    bcr_obj_sub[[]] %>%
    select(-any_of(skip_cols))

  if (ncol(bcr_meta_extra) > 0) {
    gex_meta_extra <- left_join(gex_obj[["cell_id"]], bcr_meta_extra, by = "cell_id")
    gex_obj <- AddMetaData(gex_obj, metadata = gex_meta_extra)
  }

  # overwrite Has_BCR to reflect actual BCR assay presence in the merged object,
  # regardless of whether it existed before (e.g. from a prior gex_add_airr() call)
  gex_obj$Has_BCR <- Cells(gex_obj) %in% gex_barcodes_for_bcr

  gex_obj@meta.data <- gex_obj@meta.data %>% droplevels()

  # transfer BCR reductions and graphs (inner join only, or left join with no gaps)
  if (transfer_reductions) {
    for (rd in c("bpca", "bcr.umap")) {
      if (rd %in% names(bcr_obj_sub@reductions)) {
        gex_obj[[rd]] <- bcr_obj_sub[[rd]]
      }
    }
    for (nb in "BCR.nn") {
      if (nb %in% names(bcr_obj_sub@neighbors)) {
        gex_obj@neighbors[[nb]] <- bcr_obj_sub@neighbors[[nb]]
      }
    }
    for (gr in c("BCR_nn", "BCR_snn")) {
      if (gr %in% names(bcr_obj_sub@graphs)) {
        gex_obj@graphs[[gr]] <- bcr_obj_sub@graphs[[gr]]
      }
    }
  }

  if (verbose) {
    n_shared <- length(shared_cell_ids)
    n_gex_only <- n_gex - n_shared
    n_bcr_only <- n_bcr - n_shared
    coverage <- label_percent(accuracy = 0.1)(n_shared / n_gex)
    cli::cli_inform(c(
      "i" = "{n_shared} shared cells ({coverage} of {n_gex} GEX cells).",
      if (n_gex_only > 0)
        "i" = "{n_gex_only} GEX cell{?s} had no BCR match \\
 and were {if (join == 'inner') 'excluded' else 'retained with zero BCR embeddings'}.",
      if (n_bcr_only > 0)
        "i" = "{n_bcr_only} BCR cell{?s} had no GEX match and were excluded.",
      "v" = "Merged object has {ncol(gex_obj)} cells with both RNA and BCR assays."
    ))
  }

  gex_obj
}
