load_manual_object <- function() {
  readRDS(testthat::test_path("fixtures", "manual.rds"))
}

# Load fixture and add minimal BCR metadata columns for concatenation tests
make_bcr_fixture <- function() {
  obj <- load_manual_object()
  set.seed(42)
  obj$mu_freq <- runif(ncol(obj), 0, 0.05)
  obj$isotype <- factor(sample(c("IgM", "IgG", "IgA"), ncol(obj), replace = TRUE))
  obj
}


# concatenate_gex_bcr -- argument validation ####

test_that("concatenate_gex_bcr errors on non-Seurat input", {
  expect_error(concatenate_gex_bcr(list()), "must be a Seurat object")
})

test_that("concatenate_gex_bcr errors when cols_to_include is missing for features mode", {
  obj <- make_bcr_fixture()
  expect_error(
    concatenate_gex_bcr(obj, stage = "raw", input_type = "features"),
    "cols_to_include"
  )
})


# concatenate_gex_bcr -- embed path ####

# For the embed path the GEX and BCR blocks are each centered and then divided
# by their Frobenius norm before being column-bound.  The test below verifies
# that both blocks end up with sum-of-squares == 1, i.e. neither modality
# dominates the combined space.
test_that("concatenate_gex_bcr reduced_both path normalizes GEX and BCR blocks to equal Frobenius norm", {
  obj <- make_bcr_fixture()
  n_gex_dims <- 5

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "reduced_both", input_type = "features",
                        cols_to_include = c("mu_freq", "isotype"),
                        gex_reduction = "rpca", num_dims = n_gex_dims,
                        k_param = 5, verbose = FALSE)
  ))

  emb <- Embeddings(result, "rna_bcr.pca")
  gex_frob_sq <- sum(emb[, seq_len(n_gex_dims)]^2)
  bcr_frob_sq <- sum(emb[, seq(n_gex_dims + 1, ncol(emb))]^2)

  expect_equal(gex_frob_sq, 1, tolerance = 1e-6)
  expect_equal(bcr_frob_sq, 1, tolerance = 1e-6)
})

test_that("concatenate_gex_bcr reduced_both path appends BCR dims to the combined PCA", {
  obj <- make_bcr_fixture()
  n_gex_dims <- 5

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "reduced_both", input_type = "features",
                        cols_to_include = c("mu_freq", "isotype"),
                        gex_reduction = "rpca", num_dims = n_gex_dims,
                        k_param = 5, verbose = FALSE)
  ))

  emb <- Embeddings(result, "rna_bcr.pca")

  # combined PCA has more dims than GEX alone → BCR was appended
  expect_gt(ncol(emb), n_gex_dims)
  # UMAP was run on the combined space
  expect_true("rna_bcr.umap" %in% names(result@reductions))
})

test_that("concatenate_gex_bcr reduced_both path BCR columns have non-zero per-column variance", {
  obj <- make_bcr_fixture()
  n_gex_dims <- 5

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "reduced_both", input_type = "features",
                        cols_to_include = c("mu_freq", "isotype"),
                        gex_reduction = "rpca", num_dims = n_gex_dims,
                        k_param = 5, verbose = FALSE)
  ))

  emb <- Embeddings(result, "rna_bcr.pca")
  bcr_block <- emb[, seq(n_gex_dims + 1, ncol(emb)), drop = FALSE]

  expect_true(all(apply(bcr_block, 2, var) > 0))
})


# concatenate_gex_bcr -- raw path ####

test_that("concatenate_gex_bcr raw path creates RNA_BCR assay larger than RNA alone", {
  obj <- make_bcr_fixture()
  n_rna_features <- nrow(obj[["RNA"]])

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "raw", input_type = "features",
                        cols_to_include = "mu_freq",
                        normalize = FALSE,
                        num_pcs = 10, num_dims = 5, k_param = 5,
                        verbose = FALSE)
  ))

  expect_true("RNA_BCR" %in% names(result@assays))
  expect_gt(nrow(result[["RNA_BCR"]]), n_rna_features)
})

test_that("concatenate_gex_bcr raw path includes BCR features in VariableFeatures", {
  obj <- make_bcr_fixture()

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "raw", input_type = "features",
                        cols_to_include = "mu_freq",
                        normalize = FALSE,
                        num_pcs = 10, num_dims = 5, k_param = 5,
                        verbose = FALSE)
  ))

  # "mu_freq" becomes "mu-freq-scaled" after process_bcr_features; "." matches "-"
  expect_true(any(grepl("mu.freq", VariableFeatures(result))))
})

test_that("concatenate_gex_bcr raw path RNA variable features are preserved", {
  obj <- make_bcr_fixture()
  orig_var_feats <- VariableFeatures(obj, assay = "RNA")

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "raw", input_type = "features",
                        cols_to_include = "mu_freq",
                        normalize = FALSE,
                        num_pcs = 10, num_dims = 5, k_param = 5,
                        verbose = FALSE)
  ))

  # all original RNA variable features still present
  expect_true(all(orig_var_feats %in% VariableFeatures(result)))
})


# concatenate_gex_bcr -- reduced path ####

test_that("concatenate_gex_bcr reduced_gex path creates RNA_BCR assay with GEX PCs and BCR features", {
  obj <- make_bcr_fixture()
  n_gex_dims <- 5
  # process_bcr_features on a single numeric column yields one feature
  n_bcr_features <- 1

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "reduced_gex", input_type = "features",
                        cols_to_include = "mu_freq",
                        gex_reduction = "rpca", num_dims = n_gex_dims,
                        num_pcs = 5, k_param = 5, verbose = FALSE)
  ))

  expect_true("RNA_BCR" %in% names(result@assays))
  expect_equal(nrow(result[["RNA_BCR"]]), n_gex_dims + n_bcr_features)
  expect_true("rna_bcr.pca" %in% names(result@reductions))
})


# concatenate_gex_bcr -- misc metadata ####

test_that("concatenate_gex_bcr records stage and input_type in Misc slots", {
  obj <- make_bcr_fixture()

  result <- suppressMessages(suppressWarnings(
    concatenate_gex_bcr(obj, stage = "raw", input_type = "features",
                        cols_to_include = "mu_freq",
                        normalize = FALSE,
                        num_pcs = 10, num_dims = 5, k_param = 5,
                        verbose = FALSE)
  ))

  expect_equal(Misc(result, "concat_stage"), "raw")
  expect_equal(Misc(result, "concat_input_type"), "features")
})
