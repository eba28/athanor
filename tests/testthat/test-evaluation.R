seurat_obj <- readRDS(file.path(testthat::test_path(), "fixtures", "manual.rds"))
n_cells <- ncol(seurat_obj)
cell_names <- colnames(seurat_obj)
adt_features <- rownames(seurat_obj[["ADT"]])

# fill in the neighbors slot
seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "pca",
                            assay = "RNA", k.param = 20,
                            return.neighbor = TRUE, graph.name = "RNA.nn")


# .map_assay_name ####
test_that(".map_assay_name maps GEX and WNN correctly", {
  expect_equal(.map_assay_name("GEX"), "RNA")
  expect_equal(.map_assay_name("WNN"), "w")
  expect_equal(.map_assay_name("BCR"), "BCR")
  expect_equal(.map_assay_name("RNA"), "RNA")
})


# .resolve_neighbors ####
test_that(".resolve_neighbors returns nn.idx matrix for valid slot", {
  nn_idx <- .resolve_neighbors(seurat_obj, "RNA")
  expect_true(is.matrix(nn_idx))
  expect_equal(nrow(nn_idx), n_cells)
})

test_that(".resolve_neighbors errors on missing slot", {
  expect_snapshot(
    .resolve_neighbors(seurat_obj, "BCR"),
    error = TRUE
  )
})

test_that(".resolve_neighbors maps GEX to RNA.nn slot", {
  nn_gex <- .resolve_neighbors(seurat_obj, "GEX")
  nn_rna <- .resolve_neighbors(seurat_obj, "RNA")
  expect_equal(nn_gex, nn_rna)
})


# .get_adt_matrix ####
test_that(".get_adt_matrix returns cells x features matrix", {
  mat <- .get_adt_matrix(seurat_obj)
  expect_equal(nrow(mat), n_cells)
  expect_equal(ncol(mat), length(adt_features))
  expect_equal(rownames(mat), cell_names)
})

test_that(".get_adt_matrix subsets to requested features", {
  mat <- .get_adt_matrix(seurat_obj, features_adt = c("ADT-1", "ADT-2"))
  expect_equal(ncol(mat), 2)
  expect_equal(colnames(mat), c("ADT-1", "ADT-2"))
})

test_that(".get_adt_matrix errors when no features match", {
  expect_snapshot(
    .get_adt_matrix(seurat_obj, features_adt = "CD999"),
    error = TRUE
  )
})


# calc_adt_dists ####
test_that("calc_adt_dists returns named per-cell vector", {
  result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                           features_adt = "ADT-1", k = 20,
                           return_mean = FALSE)
  expect_length(result, n_cells)
  expect_named(result, cell_names)
  expect_true(all(result >= 0, na.rm = TRUE))
})

test_that("calc_adt_dists return_mean = TRUE gives single numeric", {
  result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                           features_adt = "ADT-1", k = 20,
                           return_mean = TRUE)
  expect_length(result, 1)
  expect_true(is.numeric(result))
  expect_true(result >= 0)
})

test_that("calc_adt_dists works across all distance metrics", {
  for (metric in c("mean_abs", "manhattan", "euclidean")) {
    result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                             features_adt = "ADT-1", k = 20,
                             distance_metric = metric, return_mean = TRUE)
    expect_length(result, 1)
    expect_true(result >= 0)
  }
})

test_that("calc_adt_dists works with multiple features", {
  result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                           features_adt = c("ADT-1", "ADT-2"), k = 20,
                           return_mean = FALSE)
  expect_length(result, n_cells)
})

test_that("calc_adt_dists errors on missing neighbor slot", {
  expect_snapshot(
    calc_adt_dists(seurat_obj, base_assay = "BCR",
                   features_adt = "ADT-1", k = 20),
    error = TRUE
  )
})

test_that("calc_adt_dists errors when no features match", {
  expect_snapshot(
    calc_adt_dists(seurat_obj, base_assay = "RNA",
                   features_adt = "CD999", k = 20),
    error = TRUE
  )
})

test_that("calc_adt_dists GEX alias gives same result as RNA", {
  result_rna <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                               features_adt = "ADT-1", k = 20,
                               return_mean = TRUE)
  result_gex <- calc_adt_dists(seurat_obj, base_assay = "GEX",
                               features_adt = "ADT-1", k = 20,
                               return_mean = TRUE)
  expect_equal(result_rna, result_gex)
})


# calc_adt_nn_within_range ####
test_that("calc_adt_nn_within_range returns named proportions in [0, 1]", {
  result <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                     base_assay = "RNA", k = 20)
  expect_length(result, n_cells)
  expect_named(result, cell_names)
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("calc_adt_nn_within_range return_counts gives non-negative integers", {
  result <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                     base_assay = "RNA", k = 20,
                                     return_counts = TRUE)
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 20, na.rm = TRUE))
})

test_that("calc_adt_nn_within_range stricter range gives fewer matches", {
  loose <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                    base_assay = "RNA", k = 20, range = 0.5)
  strict <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                     base_assay = "RNA", k = 20, range = 0.05)
  expect_true(mean(loose, na.rm = TRUE) >= mean(strict, na.rm = TRUE))
})


# calc_adt_quantile ####
test_that("calc_adt_quantile method quantile returns named proportions in [0, 1]", {
  result <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                              base_assay = "RNA", k = 20,
                              method = "quantile")
  expect_length(result, n_cells)
  expect_named(result, cell_names)
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("calc_adt_quantile method percentile_diff returns values in [0, 1]", {
  result <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                              base_assay = "RNA", k = 20,
                              method = "percentile_diff")
  expect_length(result, n_cells)
  expect_named(result, cell_names)
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("calc_adt_quantile warns when return_counts used with percentile_diff", {
  expect_snapshot(
    calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                      base_assay = "RNA", k = 20,
                      method = "percentile_diff", return_counts = TRUE)
  )
})

test_that("calc_adt_quantile return_counts gives integer-valued counts", {
  result <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                              base_assay = "RNA", k = 20,
                              method = "quantile", return_counts = TRUE)
  expect_true(all(result == floor(result), na.rm = TRUE))
  expect_true(all(result >= 0 & result <= 20, na.rm = TRUE))
})


# calc_correlation ####
test_that("Correlation doesn't work on non-existent features", {
  expect_error(calc_correlation(seurat_obj, features_adt = "CD19"),
               regexp = "The requested ADT feature is not present in assay")
})

test_that("Need the neighbors slot to not be empty", {
  expect_error(calc_correlation(seurat_obj, features_adt = "ADT-1"),
               regexp = "No neighbor graphs found in object.")
})


# calc_moran ####
test_that("calc_moran returns a single numeric in [-1, 1]", {
  result <- calc_moran(seurat_obj, features_adt = "ADT-1",
                       graph_name = "RNA_nn")
  expect_length(result, 1)
  expect_true(is.numeric(result))
  expect_true(result >= -1 && result <= 1)
})

test_that("calc_moran errors on missing feature", {
  expect_snapshot(
    calc_moran(seurat_obj, features_adt = "CD999", graph_name = "RNA_nn"),
    error = TRUE
  )
})

test_that("calc_moran errors on missing graph", {
  expect_snapshot(
    calc_moran(seurat_obj, features_adt = "ADT-1", graph_name = "bad_graph"),
    error = TRUE
  )
})
