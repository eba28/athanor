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

  # original cell counts
  n_gex <- ncol(gex_obj)
  n_bcr <- ncol(bcr_obj)
  n_shared <- length(shared_cell_ids)

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
    relocate(nCount_BCR, nFeature_BCR, .after = nFeature_ADT)
  rownames(seurat_obj@meta.data) <- seurat_obj$cell_id

  # filter out cells that have no BCR data
  n_missing <- sum(!obj$Has_BCR)
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
    for (graph in c("BCR_nn", "BCR_snn")) {
      if (graph %in% names(bcr_obj@graphs)) {
        seurat_obj@graphs[[graph]] <- bcr_obj@graphs[[graph]]
      }
    }
    # note that if bcr_obj was subset above, there will be no neighbors to transfer
    for (nn in "BCR.nn") {
      if (nn %in% names(bcr_obj@neighbors)) {
        seurat_obj@neighbors[[nn]] <- bcr_obj@neighbors[[nn]]
      }
    }
    for (reduc in c("bpca", "bcr.umap")) {
      if (reduc %in% names(bcr_obj@reductions)) {
        seurat_obj@reductions[[reduc]] <- bcr_obj@reductions[[reduc]]
      }
    }
  }

  # transfer the other slots
  seurat_obj@commands <- c(seurat_obj@commands, bcr_obj@commands)
  seurat_obj@misc <- c(seurat_obj@misc, bcr_obj@misc)

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
#' @param seurat_obj A Seurat object containing RNA assay and BCR metadata. If
#'   this is a merged object produced by [merge_gex_bcr()] (i.e., it has a BCR
#'   assay), `cols_to_include` may be omitted and the BCR assay features will be
#'   used directly.
#' @param pca_stage Add BCR information before PCA or after PCA.
#' @param cols_to_include Character vector of BCR metadata column names to include
#'   in the concatenated assay e.g. c("mu_freq", "isotype") or embedded dimensions.
#'   Optional when a BCR assay is already present in `seurat_obj`.
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
                                ensembl_version = NULL, verbose = TRUE) {
  pca_stage <- match.arg(pca_stage)
  # TODO: double check the formatting of filter genes
  # TODO: give an option to use a pre-defined list of filter_genes

  # input validation
  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("seurat_obj must be a Seurat object.")
  }
  if (missing(k_param)) {
    # use the k from the GEX neighbors if it exists, otherwise default to 20
    if ("RNA.nn" %in% names(seurat_obj@neighbors)) {
      nn <- seurat_obj@neighbors[["RNA.nn"]]
      k_param <- ncol(nn@nn.idx)
      cli::cli_inform(c("i" = "Using k = {k_param}", " from RNA neighbors."))
    } else {
      k_param <- 20
      cli::cli_inform(c("i" = "Using default k = {k_param}."))
    }
  }

  if (pca_stage == "Before") {
    # TODO: check if the RNA assay contains the names of the BCR assay
    # intersect(cols_to_include, rownames(GetAssayData(seurat_obj, assay = "RNA")))

    if (rlang::is_missing(cols_to_include)) {
      if (!"BCR" %in% names(seurat_obj@assays)) {
        cli::cli_abort(c(
          "Must provide {.arg cols_to_include} or a merged object \\
           with a BCR assay.",
          "i" = "Run {.fn merge_gex_bcr} first, or specify \\
                 {.arg cols_to_include}."
        ))
      }
      # BCR assay data is already features x cells and processed
      bcr_features <- GetAssayData(seurat_obj, assay = "BCR", layer = "data")
      cli::cli_inform(c("i" = "Using BCR assay features \\
        ({nrow(bcr_features)} dimensions) from merged object."))
    } else {
      # select and process BCR metadata columns
      bcr_features <- seurat_obj[[]] %>% select(all_of(cols_to_include))
      # TODO: add an option to do weighting to influence the effect of the BCRs (post-normalization)
      # TODO: check if the embeddings output has already been normalized by amulety
      bcr_features <- process_bcr_features(bcr_features)
    }
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

    # irlba throws a warning to "use a standard svd instead" when requesting more
    # than 50% of all singular value, so let's use exact SVD if that happens
    # (which is also faster when the embedding dimension is small)
    scale_data <- Seurat::GetAssayData(seurat_obj, layer = "scale.data")
    max_dim <- min(nrow(scale_data), ncol(scale_data))
    use_approx <- num_pcs < max_dim / 2
    seurat_obj <- RunPCA(object = seurat_obj, npcs = num_pcs,
                         reduction.name = "rna_bcr.pca",
                         reduction.key = "rnabcrpca_",
                         approx = use_approx)
    cli::cli_inform(c("v" = "Computed PCA with {num_pcs} dimensions using {ifelse(use_approx, 'approximate', 'exact')} SVD."))
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

  seurat_obj <- regen_reduc(seurat_obj = seurat_obj, pca_name = "rna_bcr.pca",
                            assay = "RNA_BCR", num_dims = num_dims,
                            k_param = k_param, verbose = verbose)

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
                    pc_bcr = 20, k_param = 20, cluster = FALSE,
                    cluster_res = list("GEX" = 1, "BCR" = 1, "WNN" = 1),
                    modality_weights = NULL, verbose = TRUE) {
  # TODO: update this to be able to run on other omics e.g. GEX & ADT
  # TODO: add the option to filter out genes again???
  # TODO: give the option to use an integrated GEX assay

  # input validation
  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("seurat_obj must be a Seurat object.")
  }
  if (!"cell_id" %in% colnames(seurat_obj[[]])) {
    cli::cli_abort("Cell ID column not found in metadata.")
  }
  if (missing(pc_gex)) {
    pc_gex <- if ("rpca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["rpca"]])
    } else if ("pca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["pca"]])
    } else {
      20
    }
    cli::cli_inform(c("i" = "Using pc_gex = {pc_gex}", " from existing reductions."))
  }
  if (missing(pc_bcr)) {
    pc_bcr <- if ("bpca" %in% names(seurat_obj@reductions)) {
      ncol(seurat_obj@reductions[["bpca"]])
    } else {
      20
    }
    cli::cli_inform(c("i" = "Using pc_bcr = {pc_bcr}", " from existing reductions."))
  }
  if (missing(k_param)) {
    # use the k from the GEX neighbors if it exists, otherwise default to 20
    if ("RNA.nn" %in% names(seurat_obj@neighbors)) {
      nn <- seurat_obj@neighbors[["RNA.nn"]]
      k_param <- ncol(nn@nn.idx)
      cli::cli_inform(c("i" = "Using k = {k_param}", " from RNA neighbors."))
    } else {
      k_param <- 20
      cli::cli_inform(c("i" = "Using default k = {k_param}."))
    }
  }

  # detect if the object is already merged (has BCR assay + bpca from merge_gex_bcr)
  is_merged <- "BCR" %in% names(seurat_obj@assays) &&
    "bpca" %in% names(seurat_obj@reductions)

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
    seurat_obj[["bpca"]] <- bcr_obj[["bpca"]]
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
                                      reduction.list = list("rpca", "bpca"),
                                      dims.list = list(1:pc_gex, 1:pc_bcr),
                                      k.nn = k_param,
                                      # match the RNA and BCR style
                                      knn.graph.name = "w_nn",
                                      snn.graph.name = "w_snn",
                                      weighted.nn.name = "w.nn",
                                      modality.weight.name =
                                        str_c(c("RNA", "BCR"), ".weight"),
                                      return.intermediate = TRUE,
                                      modality.weight = modality_weights,
                                      verbose = verbose)

    # add a metadata column listing which assay was chosen for each cell
    seurat_obj@meta.data <-
      seurat_obj[[]] %>%
      mutate(weight_assay =
               case_when(seurat_obj[[]][["RNA.weight"]] > 0.5 ~ "RNA",
                         seurat_obj[[]][["RNA.weight"]] < 0.5 ~ "BCR",
                         is.na(seurat_obj[[]][["RNA.weight"]]) ~ NA,
                         .default = "Tie"))
  } else {
    # must provide a Seurat object with WNN already run so we can rebuild
    # the ModalityWeights object
    # TODO: see if there's a way to make it from scratch
    if (modality_weights == "equal") {
      mw_mod <- seurat_obj@misc$modality.weight
      mw_mod_names <- Cells(seurat_obj)

      mw_mod@modality.weight.list$rpca <-
        setNames(rep(0.5, ncol(seurat_obj)), mw_mod_names)
      mw_mod@modality.weight.list$bpca <-
        setNames(rep(0.5, ncol(seurat_obj)), mw_mod_names)

      seurat_obj <-
        Seurat::FindMultiModalNeighbors(object = seurat_obj,
                                        reduction.list = list("rpca", "bpca"),
                                        dims.list = list(1:pc_gex, 1:pc_bcr),
                                        k.nn = k_param,
                                        # match the RNA and BCR style
                                        knn.graph.name = "w_equal_nn",
                                        snn.graph.name = "w_equal_snn",
                                        weighted.nn.name = "w_equal.nn",
                                        modality.weight.name =
                                          str_c(c("RNA", "BCR"), "_equal.weight"),
                                        return.intermediate = TRUE,
                                        modality.weight = mw_mod,
                                        verbose = verbose)

      # add a metadata column listing which assay was chosen for each cell
      seurat_obj@meta.data <-
        seurat_obj[[]] %>%
        mutate(weight_assay_equal =
                 case_when(seurat_obj[[]][["RNA_equal.weight"]] > 0.5 ~ "RNA",
                           seurat_obj[[]][["RNA_equal.weight"]] < 0.5 ~ "BCR",
                           is.na(seurat_obj[[]][["RNA_equal.weight"]]) ~ NA,
                           .default = "Tie"))
    } else {
      cli::cli_abort("Please provide a valid customization approach. Only 'equal' is currently supported.")
    }
  }

  # check for NA values
  # TODO: add this check to other UMAPs
  na_nn <- is.na(seurat_obj$RNA.weight)
  if (any(na_nn)) {
    bad_cells <- seurat_obj[[]]$cell_id[na_nn]
    bcr_mat <- GetAssayData(seurat_obj, assay = "BCR", layer = "data")
    bcr_mat <- as.matrix(bcr_mat[, bad_cells])
    # add 1 because the first column is not included in the duplicate count
    dup_count <- sum(duplicated(t(bcr_mat))) + 1

    # TODO: record all duplicates, not just the ones that failed??
    Misc(seurat_obj, slot = "embedding_dups") <- bad_cells

    cli::cli_inform(c("x" = "{sum(na_nn)} NA rows found in w.nn neighbor \\
                            indices (out of {ncol(seurat_obj)} total cells). \\
                            This means that these cells do not have any valid \\
                            neighbors in the WNN space.",
                      "i" = "This is most likely due to these cells having \\
                             identical BCR embeddings. Out of the {sum(na_nn)} \\
                             cells with NA neighbors, {dup_count} have \\
                             identical BCR embeddings.",
                      ">" = "UMAP will not be run for the WNN reduction because
                      of the NA values."))
  } else {
    mw_name <- modality_weights
    if (!is.null(modality_weights)) mw_name <- paste0("_", mw_name)

    seurat_obj <-
      Seurat::RunUMAP(object = seurat_obj,
                      nn.name = paste0("w", mw_name, ".nn"),
                      n.neighbors = k_param, # might not be needed
                      reduction.name = paste0("wnn", mw_name, ".umap"),
                      # note that if modality_weights has any underscores in it,
                      # Seurat will remove them when making the key
                      reduction.key = paste0("wnn", modality_weights, "UMAP_"),
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
      cli::cli_inform(c("v" = "WNN neighbors calculated and UMAP run.",
                        "i" = "Use {.fn Seurat::DimPlot} with \\
                               {.arg reduction = 'wnn.umap'} or \\
                               {.fn athanor::plot_dimplot} with \\
                               {.arg reduc = 'wnn.umap'} to visualize the results."))
    }

    # explain the weighting approach
    if (is.null(modality_weights)) {
      cli::cli_inform(c("v" = "Modality weights were automatically calculated based on the provided assays.",
                        "i" = "Check {.fn Seurat::FindMultiModalNeighbors} for more details."))
      # summarize the resulting weight distribution with a total count
      # weight_summary <- seurat_obj[[]] %>%
      #   group_by(weight_assay) %>%
      #   summarise(count = n()) %>%
      #   ungroup()
      # cli::cli_inform(c("v" = "Summary of modality weights across cells:",
      #                   "i" = "{.code print(weight_summary)}"))
    } else {
      cli::cli_inform(c("v" = "Custom modality weights were used:",
                        "i" = "{modality_weights}"))
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
  # base message (after AIRR integration)
  message <- paste("This object has", ncol(seurat_obj), "cells, ")

  # assay details
  assay_details <- c()
  bcr_embeddings_name <- grep("^BCR", names(seurat_obj@assays), value = TRUE)

  # GEX info
  if ("RNA" %in% names(seurat_obj@assays)) {
    assay_details <-
      c(assay_details, paste(nrow(seurat_obj@assays$RNA), "genes"))
  }

  # ADT info
  if ("ADT" %in% names(seurat_obj@assays)) {
    assay_details <-
      c(assay_details,
        paste(nrow(seurat_obj@assays$ADT), "cell surface protein markers"))
  }

  # BCR info
  if (length(bcr_embeddings_name) > 0) {
    assay_details <-
      c(assay_details,
        paste(nrow(seurat_obj@assays[[bcr_embeddings_name]]),
              "embedding dimensions"))
  }

  message <- paste0(message, str_c(assay_details, collapse = ", "), ".")

  # reductions info
  if ("weighted.nn" %in% names(seurat_obj@neighbors)) {
    message <- paste(message, "WNN was run with",
                     ncol(seurat_obj@reductions[[gex_pca]]), "GEX PCs and",
                     ncol(seurat_obj@reductions[[other_pca]]), other_type, "PCs.")
  }

  # TODO: remove this?
  # clustering info
  # does not account for multiple clustering resolutions
  if ("seurat_clusters" %in% colnames(seurat_obj[[]])) {
    nclusters_rna <-
      select(seurat_obj[[]], starts_with("RNA_snn")) %>% select(last_col())
    nclusters_bcr <-
      select(seurat_obj[[]], starts_with("BCR_snn")) %>% select(last_col())
    nclusters_wnn <-
      select(seurat_obj[[]], starts_with("w_snn")) %>% select(last_col())

    message <- paste(message,
                     "There were", n_distinct(nclusters_rna), "RNA clusters,",
                     n_distinct(nclusters_bcr), "BCR clusters, and",
                     n_distinct(nclusters_wnn),
                     "WNN clusters identified (in the largest resolutions).")
    # cat(str_glue("There were {nclusters_rna} RNA clusters, {nclusters_bcr} BCR clusters, and {nclusters_wnn} WNN clusters identified."))
  }

  # output the message
  cli::cli_inform(message)
}
