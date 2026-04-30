# helper functions ####
make_airr_df <- function(barcodes, num_dims = 5, seed = 42) {
  set.seed(seed)

  df <- as.data.frame(matrix(rnorm(length(barcodes) * num_dims),
                             nrow = length(barcodes)))
  df$cell_id <- barcodes

  return(df)
}

make_airrflow_dir <- function(base_dir, version) {
  version_num <- str_split_1(version, pattern = "")[1:3] %>% str_c(collapse = "")
  results_path <- file.path(base_dir, "airrflow", "bcr", version, "results")

  rep_path <- if (version_num < 4.2) {
    file.path(results_path, "repertoire_comparison", "repertoires")
  } else if (version_num == "4.3") {
    file.path(results_path, "clonal_analysis", "define_clones",
              "all_reps_clone_report", "repertoires")
  } else {
    file.path(results_path, "clonal_analysis", "clonal_assignment",
              "all_reps_clone_report", "repertoires")
  }

  dir.create(rep_path, recursive = TRUE)
  rep_path
}

make_bcr_row <- function(cell_id, locus, c_call,
                         sample_id = "S1", subject_id = "Subj1", sex = "Female",
                         clone_id = "1", v_call = "IGHV1-2*01",
                         d_call = "IGHD3-10*01", j_call = "IGHJ4*02") {
  # 60-nt junction: divisible by 3, starts with conserved Cys (TGT), valid for
  # `aminoAcidProperties(nt = TRUE, trim = TRUE)`
  junction <- "TGTGCGAGAGATTACGGTATGGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA"

  data.frame(sequence_id = paste0(cell_id, "_", locus),
             cell_id = cell_id, locus = locus, c_call  = c_call,
             sample_id = sample_id, subject_id  = subject_id, sex = sex,
             clone_id = clone_id, junction = junction,
             v_call = v_call, d_call = d_call, j_call = j_call)
}

make_combined_airr <- function(barcodes, sample = "sample") {
  data.frame(cell_id = paste0(sample, "_", sub("-.*", "", barcodes)),
             cell_id_original = barcodes)
}

make_embeddings <- function(num_dims = 100, num_cells = 100, seed = 42) {
  set.seed(seed)

  m <- matrix(rnorm(num_dims * num_cells), nrow = num_dims, ncol = num_cells)
  rownames(m) <- paste0("Dim-", seq_len(num_dims))
  colnames(m) <- paste0("Cell-", seq_len(num_cells))

  return(m)
}

make_seurat_mu_freq <- function(mu_freqs, seed = 42) {
  set.seed(seed)

  counts <- matrix(rpois(length(mu_freqs) * 5, 1), nrow = 5,
                   dimnames = list(paste0("G", seq_len(5)),
                                   paste0("C", seq_along(mu_freqs))))
  obj <- CreateSeuratObject(counts = as(counts, "dgCMatrix"))
  obj$mu_freq <- mu_freqs

  return(obj)
}

write_airrflow_tsv <- function(rep_path, data, filename = "reps.tsv") {
  readr::write_tsv(data, file.path(rep_path, filename))
}

# test data ####
embeddings <- make_embeddings(num_dims = 100, num_cells = 100, seed = 42)
bcr_obj <- bcr_embeddings_pipeline(embeddings, embedding_type = "Simulated",
                                   num_dims = 10, k_param = 10, verbose = FALSE)


# add_family_info ####
test_that("add_family_info adds v_call columns when v_call is present", {
  df <- data.frame(v_call = c("IGHV1-2*01", "IGHV3-30*01"))
  result <- add_family_info(df)

  expect_true("v_call_family" %in% names(result))
  expect_true("v_call_gene" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("add_family_info adds all six columns when v, d, j calls are present", {
  df <- data.frame(v_call = c("IGHV1-2*01", "IGHV3-30*01"),
                   d_call = c("IGHD3-10*01", "IGHD6-19*01"),
                   j_call = c("IGHJ4*02", "IGHJ6*02"))
  result <- add_family_info(df)
  expected_cols <- c("v_call_family", "v_call_gene",
                     "d_call_family", "d_call_gene",
                     "j_call_family", "j_call_gene")

  expect_true(all(expected_cols %in% names(result)))
})

test_that("add_family_info skips absent call columns", {
  df <- data.frame(v_call = c("IGHV1-2*01"))
  result <- add_family_info(df)

  expect_false("d_call_family" %in% names(result))
  expect_false("j_call_gene" %in% names(result))
})


# bcr_embeddings_pipeline ####
test_that("bcr_embeddings_pipeline returns a Seurat object with a BCR assay", {
  expect_s4_class(bcr_obj, "Seurat")
  expect_true("BCR" %in% names(bcr_obj@assays))
  expect_equal(ncol(bcr_obj), ncol(embeddings))
})

test_that("bcr_embeddings_pipeline produces expected reductions and graphs", {
  expect_true("bpca" %in% names(bcr_obj@reductions))
  expect_true("bcr.umap" %in% names(bcr_obj@reductions))
  expect_true("BCR.nn" %in% names(bcr_obj@neighbors))
  expect_true("BCR_nn" %in% names(bcr_obj@graphs))
  expect_true("BCR_snn" %in% names(bcr_obj@graphs))
})

test_that("bcr_embeddings_pipeline stores Misc slots correctly", {
  expect_equal(Seurat::Misc(bcr_obj, slot = "embedding_type"), "Simulated")
  expect_equal(Seurat::Misc(bcr_obj, slot = "embedding_dims"), nrow(embeddings))
})

test_that("bcr_embeddings_pipeline sets cell_id to Cells()", {
  expect_true("cell_id" %in% colnames(bcr_obj[[]]))
  expect_equal(as.character(bcr_obj$cell_id), Cells(bcr_obj))
})

test_that("bcr_embeddings_pipeline sets scale.data from embeddings directly", {
  scale_data <- Seurat::GetAssayData(bcr_obj, assay = "BCR",
                                     layer = "scale.data")

  expect_equal(as.matrix(scale_data), embeddings)
})


# bin_mu_freq ####
test_that("bin_mu_freq adds mu_freq_bins with 5 bins", {
  obj <- make_seurat_mu_freq(c(0, 0.005, 0.025, 0.07, 0.15))
  result <- bin_mu_freq(obj, num_bins = 5)

  expect_true("mu_freq_bins" %in% colnames(result[[]]))
  expect_s3_class(result$mu_freq_bins, "factor")
  expect_equal(nlevels(result$mu_freq_bins), 5)
})

test_that("bin_mu_freq adds mu_freq_bins_fewer with 3 bins", {
  obj <- make_seurat_mu_freq(c(0, 0.005, 0.02))
  result <- bin_mu_freq(obj, num_bins = 3)

  expect_true("mu_freq_bins_fewer" %in% colnames(result[[]]))
  expect_equal(nlevels(result$mu_freq_bins_fewer), 3)
})

test_that("bin_mu_freq adds mu_freq_bins_binary with 2 bins", {
  obj <- make_seurat_mu_freq(c(0.005, 0.02, 0.1))
  result <- bin_mu_freq(obj, num_bins = 2)

  expect_true("mu_freq_bins_binary" %in% colnames(result[[]]))
  expect_equal(nlevels(result$mu_freq_bins_binary), 2)
})

test_that("bin_mu_freq adds all three column sets when all bins requested", {
  obj <- make_seurat_mu_freq(c(0, 0.005, 0.025, 0.07, 0.15))
  result <- bin_mu_freq(obj, num_bins = c(2, 3, 5))

  expect_true(all(c("mu_freq_bins", "mu_freq_bins_fewer", "mu_freq_bins_binary") %in%
                    colnames(result[[]])))
})


# convert_embeddings ####
test_that("convert_embeddings returns a sparse Matrix with features as rows", {
  barcodes <- paste0("AAAC", seq_len(5), "-1")
  combined_airr <- make_combined_airr(barcodes)
  embeddings_df <- make_airr_df(barcodes, num_dims = 8)

  result <- suppressMessages(convert_embeddings(embeddings_df, combined_airr))

  expect_true(inherits(result, "Matrix"))
  expect_equal(ncol(result), 5)
  expect_equal(nrow(result), 8)
})

test_that("convert_embeddings names dimensions Dim1..DimN", {
  barcodes <- paste0("AAAC", seq_len(3), "-1")
  combined_airr <- make_combined_airr(barcodes)
  embeddings_df <- make_airr_df(barcodes, num_dims = 4)

  result <- suppressMessages(convert_embeddings(embeddings_df, combined_airr))

  expect_equal(rownames(result), paste0("Dim", seq_len(4)))
})

test_that("convert_embeddings removes cells absent from combined_airr", {
  barcodes <- paste0("AAAC", seq_len(5), "-1")
  combined_airr <- make_combined_airr(barcodes[1:3])
  embeddings_df <- make_airr_df(barcodes, num_dims = 4)

  result <- suppressMessages(convert_embeddings(embeddings_df, combined_airr))

  expect_equal(ncol(result), 3)
})

test_that("convert_embeddings errors when required columns missing from combined_airr", {
  bad_airr <- data.frame(cell_id = "x", wrong_col = "y")
  embeddings_df <- data.frame(cell_id = "x", V1 = 1.0)

  expect_error(convert_embeddings(embeddings_df, bad_airr))
})


# factor_family_info ####
test_that("factor_family_info converts v_call_family to a numeric-sorted factor", {
  df <- data.frame(v_call = c("x", "x", "x"),
                   v_call_family = c("IGHV1", "IGHV10", "IGHV2"),
                   v_call_gene = c("IGHV1-2", "IGHV10-3", "IGHV2-5"))
  result <- factor_family_info(df)

  expect_s3_class(result$v_call_family, "factor")
  expect_equal(levels(result$v_call_family), c("IGHV1", "IGHV2", "IGHV10"))
})

test_that("factor_family_info handles all three gene loci", {
  df <- data.frame(v_call = "x", v_call_family = "IGHV3",
                   v_call_gene = "IGHV3-1",
                   d_call = "x", d_call_family = "IGHD2",
                   d_call_gene = "IGHD2-2",
                   j_call = "x", j_call_family = "IGHJ1", j_call_gene = "IGHJ1")
  result <- factor_family_info(df)

  expect_s3_class(result$v_call_family, "factor")
  expect_s3_class(result$d_call_family, "factor")
  expect_s3_class(result$j_call_family, "factor")
})

test_that("factor_family_info skips absent call columns", {
  df <- data.frame(v_call = "x", v_call_family = "IGHV1",
                   v_call_gene = "IGHV1-2")
  result <- factor_family_info(df)

  expect_false("d_call_family" %in% names(result))
})


# process_airrflow ####
test_that("process_airrflow errors when no repertoire files are found", {
  dir <- withr::local_tempdir()
  make_airrflow_dir(dir, "4.0")

  expect_error(process_airrflow(dir, "4.0"), "No files found")
})

test_that("process_airrflow creates clone_id_unique as subject_id + clone_id", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01", subject_id = "S1",
                             clone_id = "7"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC", subject_id = "S1",
                             clone_id = "7"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_true("clone_id_unique" %in% colnames(result))
  expect_equal(result$clone_id_unique[result$locus == "IGH"], "S1_7")
})

test_that("process_airrflow stores original cell_id and adds sample_id prefix", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01", sample_id = "S1"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC", sample_id = "S1"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_true("cell_id_original" %in% colnames(result))
  expect_equal(result$cell_id_original[result$locus == "IGH"], "AAAC1-1")
  expect_equal(result$cell_id[result$locus == "IGH"], "S1_AAAC1")
})

test_that("process_airrflow derives isotype from c_call when column is absent", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHG1*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_true("isotype" %in% colnames(result))
  expect_equal(result$isotype[result$locus == "IGH"], "IgG")
})

test_that("process_airrflow converts sex column to character", {
  # TODO: write this
})

test_that("process_airrflow preserves a pre-existing isotype column unchanged", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  data$isotype <- c("IgM", NA)
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_equal(result$isotype[result$locus == "IGH"], "IgM")
})

test_that("process_airrflow filters out IGH rows with NA c_call", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = NA),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"),
                make_bcr_row(cell_id = "AAAC2-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC2-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_false("S1_AAAC1" %in% result$cell_id)
  expect_true("S1_AAAC2" %in% result$cell_id)
})

test_that("process_airrflow filters out light chains with no paired heavy chain", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"),
                make_bcr_row(cell_id = "AAAC2-1", locus = "IGK",
                             v_call = "IGKC"))  # no IGH
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_false("S1_AAAC2" %in% result$cell_id)
  expect_true("S1_AAAC1" %in% result$cell_id)
})

test_that("process_airrflow filters out IGH rows with a non-IGH c_call (multi chains)", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"),
                make_bcr_row(cell_id = "AAAC2-1", locus = "IGH",
                             v_call = "IGLC1"),  # light c_call on IGH
                make_bcr_row(cell_id = "AAAC2-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_false("S1_AAAC2" %in% result$cell_id)
  expect_true("S1_AAAC1" %in% result$cell_id)
})

test_that("process_airrflow stores original c_call in c_call_original", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_true("c_call_original" %in% colnames(result))
  expect_equal(result$c_call_original[result$locus == "IGH"], "IGHM*01")
})

test_that("process_airrflow simplifies c_call to gene family level", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_equal(result$c_call[result$locus == "IGH"], "IGHM")
})

test_that("process_airrflow adds v_call_family and v_call_gene columns", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_true("v_call_family" %in% colnames(result))
  expect_true("v_call_gene" %in% colnames(result))
})

test_that("process_airrflow combines rows from multiple TSV files", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.0")
  data1 <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                              v_call = "IGHM*01", sample_id = "S1"),
                 make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                              v_call = "IGKC", sample_id = "S1"))
  data2 <- rbind(make_bcr_row(cell_id = "AAAC2-1", locus = "IGH",
                              v_call = "IGHA1", sample_id = "S2"),
                 make_bcr_row(cell_id = "AAAC2-1", locus = "IGK",
                              v_call = "IGKC", sample_id = "S2"))
  write_airrflow_tsv(rep_path, data1, "sub1.tsv")
  write_airrflow_tsv(rep_path, data2, "sub2.tsv")

  result <- suppressMessages(process_airrflow(dir, "4.0"))

  expect_true("S1_AAAC1" %in% result$cell_id)
  expect_true("S2_AAAC2" %in% result$cell_id)
})

test_that("process_airrflow uses the v4.3 path for version 4.3.1", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "4.3.1")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "4.3.1"))

  expect_true("S1_AAAC1" %in% result$cell_id)
})

test_that("process_airrflow uses the v5 path for version 5.0.0", {
  dir <- withr::local_tempdir()
  rep_path <- make_airrflow_dir(dir, "5.0.0")
  data <- rbind(make_bcr_row(cell_id = "AAAC1-1", locus = "IGH",
                             v_call = "IGHM*01"),
                make_bcr_row(cell_id = "AAAC1-1", locus = "IGK",
                             v_call = "IGKC"))
  write_airrflow_tsv(rep_path, data)

  result <- suppressMessages(process_airrflow(dir, "5.0.0"))

  expect_true("S1_AAAC1" %in% result$cell_id)
})


# process_bcr_features ####
test_that("process_bcr_features returns a matrix normalized to mean 0", {
  df <- data.frame(mu_freq = c(0.01, 0.05, 0.10, 0.02, 0.08))
  result <- suppressMessages(process_bcr_features(df))

  expect_true(is.matrix(result))
  expect_lt(abs(mean(result)), 1e-10)
})

test_that("process_bcr_features removes underscores from row names", {
  df <- data.frame(mu_freq = c(0.01, 0.05, 0.10))
  result <- suppressMessages(process_bcr_features(df))

  expect_false(any(grepl("_", rownames(result))))
})

test_that("process_bcr_features one-hot encodes unordered categorical variables", {
  df <- data.frame(isotype = c("IgG", "IgM", "IgA", "IgG", "IgM"))
  result <- suppressMessages(process_bcr_features(df))

  expect_true(nrow(result) >= 3)
})

test_that("process_bcr_features converts ordered factors with -ordered suffix", {
  df <- data.frame(mu_bins = factor(c("low", "med", "high", "low", "high"),
                   levels = c("low", "med", "high"), ordered = TRUE))
  result <- suppressMessages(process_bcr_features(df))

  expect_true(is.matrix(result))
  expect_true(any(grepl("ordered", rownames(result))))
})

test_that("process_bcr_features drops zero-variance features", {
  # constant column has zero variance and should be removed by step_zv
  df <- data.frame(mu_freq = c(0.01, 0.05, 0.10, 0.02, 0.08),
                   constant = rep(1, 5))
  result <- suppressMessages(process_bcr_features(df))

  expect_false(any(grepl("constant", rownames(result))))
})
