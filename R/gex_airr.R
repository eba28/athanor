#' Add AIRR information to a Seurat object
#'
#' @description
#' This function integrates adaptive immune receptor repertoire (AIRR) data with
#' gene expression (GEX) data in a Seurat object. Currently built to integrate
#' BCR data, including both heavy and light chain information.
#'
#' @details
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
#' @param overview Whether to print integration summary information.
#'
#' @returns The Seurat object with AIRR columns added to the metadata, including
#'   Has_BCR, isotype information, mutation frequencies, and pairing status.
#' @export
gex_add_airr <- function(seurat_obj, airr_type = "BCR",
                         combined_airr, new_cols, overview = TRUE) {
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
    seurat_obj@meta.data <-
      seurat_obj[[]] %>%
      mutate(annotated_clusters_bcr =
               ifelse(mu_freq <= 0.01 & isotype %in% c("IgM", "IgD"),
                      "Naive B Cells", "Non-Naive B Cells"),
             # assumes that the plasma cells don't need re-assigning
             annotated_clusters_gex_bcr =
               ifelse(annotated_clusters_simpler == "Naive B cells" &
                        (mu_freq > 0.01 | isotype %in% c("IgA", "IgE", "IgG")),
                      "Memory B cells", annotated_clusters_simpler))
  }

  # make sure that the levels are okay
  seurat_obj@meta.data <- seurat_obj@meta.data %>% droplevels()

  # print an overview of the integration
  if (overview) {
    count_gex_airr <- nrow(filter(combined_airr_select, get(airr_col)))
    count_paired <- nrow(filter(combined_airr_select, get(paired_col) == "TRUE"))
    count_filtered <- nrow(filter(combined_airr_select_obj, get(airr_col)))
    count_unmatched <- count_gex_airr - count_filtered
    coverage <- mean(Cells(seurat_obj) %in% combined_airr_select$cell_id)
    count_excluded <- 1 - mean(combined_airr_select$cell_id %in% Cells(seurat_obj))

    cat(sprintf("There are %d %s chain %s being integrated, %d of which have at least one paired %s chain. When filtered to just the ones that match by cell_id, they cover %d cells (%s of the GEX data (%d cells)). %d of these %s (%s of the total %s chains) did not have a match with the GEX cell ids and were excluded.\n",
                count_gex_airr, chain_main, airr_type, count_paired, chain_other,
                count_filtered, label_percent(accuracy = 0.1)(coverage),
                ncol(seurat_obj), count_unmatched, airr_type,
                label_percent(accuracy = 0.1)(count_excluded), count_gex_airr))

    # double check the columns that were added
    cat(paste0("The following new columns were added to the Seurat object's metadata: ",
               str_c(setdiff(colnames(seurat_obj[[]]), meta_cols_orig),
                     collapse = ", "), "."))
  }


  # return the updated Seurat object
  seurat_obj
}
