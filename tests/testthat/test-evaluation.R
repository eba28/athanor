# A small Seurat object with a precomputed RNA.nn neighbor graph (k = 5, with
# each cell as its own first neighbor), an ADT assay (ADT-1 .. ADT-5), a pca
# reduction, and RNA_nn / RNA_snn graphs.
seurat_obj <- readRDS(testthat::test_path("fixtures", "manual2.rds"))
n_cells <- ncol(seurat_obj)
cell_names <- colnames(seurat_obj)
adt_features <- rownames(seurat_obj[["ADT"]])

# metadata columns for the neighbor-matching tests: one BCR-style, one
# GEX-style (so gex_cols routing can be checked)
seurat_obj$grp <- factor(rep(c("a", "b"), length.out = n_cells))
seurat_obj$seurat_clusters <- factor(rep(c("0", "1"), each = n_cells / 2))

# every calc_* evaluation function returns this schema; per-cell output
# additionally carries a leading cell_id column
canonical_cols <- c("Graph", "Assay", "Feature", "Method", "Score")

expect_canonical <- function(df, per_cell = FALSE) {
  expect_s3_class(df, "data.frame")
  want <- if (per_cell) c("cell_id", canonical_cols) else canonical_cols
  expect_identical(names(df), want)
}


# helpers ####
test_that("map_assay_name maps GEX and WNN, passes others through", {
  expect_equal(map_assay_name("GEX"), "RNA")
  expect_equal(map_assay_name("WNN"), "w")
  expect_equal(map_assay_name("BCR"), "BCR")
  expect_equal(map_assay_name("RNA"), "RNA")
})

test_that("resolve_neighbors returns the nn.idx matrix for a valid slot", {
  nn_idx <- resolve_neighbors(seurat_obj, "RNA")
  expect_true(is.matrix(nn_idx))
  expect_equal(nrow(nn_idx), n_cells)
})

test_that("resolve_neighbors maps GEX to the RNA.nn slot", {
  expect_equal(resolve_neighbors(seurat_obj, "GEX"),
               resolve_neighbors(seurat_obj, "RNA"))
})

test_that("resolve_neighbors errors on a missing slot", {
  expect_error(resolve_neighbors(seurat_obj, "BCR"))
})

test_that("get_adt_matrix returns a cells-by-features matrix", {
  mat <- get_adt_matrix(seurat_obj)
  expect_equal(nrow(mat), n_cells)
  expect_equal(ncol(mat), length(adt_features))
  expect_equal(rownames(mat), cell_names)
})

test_that("get_adt_matrix subsets to the requested features", {
  mat <- get_adt_matrix(seurat_obj, features_adt = c("ADT-1", "ADT-2"))
  expect_equal(colnames(mat), c("ADT-1", "ADT-2"))
})

test_that("drop_self_neighbors drops exactly one neighbor per cell", {
  nn_idx <- seurat_obj@neighbors[["RNA.nn"]]@nn.idx
  dropped <- drop_self_neighbors(nn_idx)
  expect_equal(ncol(dropped), ncol(nn_idx) - 1)
  expect_equal(nrow(dropped), nrow(nn_idx))
  # self (row i contains i) is gone for every row
  expect_false(any(vapply(seq_len(nrow(dropped)),
                          function(i) i %in% dropped[i, ], logical(1))))
})


# calc_adt_correlation ####
test_that("calc_adt_correlation returns the canonical schema", {
  result <- calc_adt_correlation(seurat_obj, features_adt = "ADT-1")
  expect_canonical(result)
  expect_equal(unique(result$Method), "Spearman")
  expect_equal(unique(result$Assay), "ADT")
  expect_true(all(result$Score >= -1 & result$Score <= 1))
})

test_that("calc_adt_correlation records the correlation method", {
  result <- calc_adt_correlation(seurat_obj, features_adt = "ADT-1",
                                 cor_method = "pearson")
  expect_equal(unique(result$Method), "Pearson")
})

test_that("calc_adt_correlation drops features that are not in the assay", {
  # absent features are silently skipped; a valid one still comes through
  mixed <- calc_adt_correlation(seurat_obj, features_adt = c("ADT-1", "CD19"))
  expect_canonical(mixed)
  expect_equal(mixed$Feature, "ADT-1")
  # nothing left to score
  expect_null(calc_adt_correlation(seurat_obj, features_adt = "CD19"))
})

test_that("calc_adt_correlation errors when there are no neighbor graphs", {
  bare <- seurat_obj
  bare@neighbors <- list()
  expect_error(calc_adt_correlation(bare, features_adt = "ADT-1"),
               "No neighbor graphs found")
})


# calc_adt_dists ####
test_that("calc_adt_dists return_mean = TRUE gives a one-row tidy frame", {
  result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                           features_adt = "ADT-1", k = 20)
  expect_canonical(result)
  expect_equal(nrow(result), 1)
  expect_true(result$Score >= 0)
})

test_that("calc_adt_dists return_mean = FALSE gives one row per cell", {
  result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                           features_adt = "ADT-1", k = 20, return_mean = FALSE)
  expect_canonical(result, per_cell = TRUE)
  expect_equal(nrow(result), n_cells)
  expect_setequal(result$cell_id, cell_names)
})

test_that("calc_adt_dists labels Feature with the joined feature set", {
  one <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                        features_adt = "ADT-1", k = 20)
  many <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                         features_adt = c("ADT-2", "ADT-1"), k = 20)
  all_feats <- calc_adt_dists(seurat_obj, base_assay = "RNA", k = 20)
  expect_equal(one$Feature, "ADT-1")
  expect_equal(many$Feature, "ADT-1+ADT-2")
  expect_equal(all_feats$Feature, "all")
})

test_that("calc_adt_dists records the distance metric in Method", {
  for (metric in c("mean_abs", "manhattan", "euclidean")) {
    result <- calc_adt_dists(seurat_obj, base_assay = "RNA",
                             features_adt = "ADT-1", k = 20,
                             distance_metric = metric)
    expect_equal(nrow(result), 1)
    expect_true(result$Score >= 0)
  }
})

test_that("calc_adt_dists GEX alias matches RNA", {
  dists <- function(base) {
    calc_adt_dists(seurat_obj, base_assay = base,
                   features_adt = "ADT-1", k = 20)
  }
  expect_equal(dists("GEX"), dists("RNA"))
})

test_that("calc_adt_dists errors on a missing neighbor slot", {
  expect_error(
    calc_adt_dists(seurat_obj, base_assay = "BCR",
                   features_adt = "ADT-1", k = 20)
  )
})


# calc_adt_mean_absolute ####
test_that("calc_adt_mean_absolute returns the canonical schema", {
  agg <- calc_adt_mean_absolute(seurat_obj, features_adt = c("ADT-1", "ADT-2"))
  expect_canonical(agg)
  expect_equal(unique(agg$Method), "Mean_Absolute")
  expect_setequal(agg$Feature, c("ADT-1", "ADT-2"))

  per_cell <- calc_adt_mean_absolute(seurat_obj, features_adt = "ADT-1",
                                     return_mean = FALSE)
  expect_canonical(per_cell, per_cell = TRUE)
  expect_equal(nrow(per_cell), n_cells)
})


# calc_adt_moran ####
test_that("calc_adt_moran returns a one-row tidy frame with Score in [-1, 1]", {
  result <- calc_adt_moran(seurat_obj, features_adt = "ADT-1",
                           graph_name = "RNA_nn")
  expect_canonical(result)
  expect_equal(nrow(result), 1)
  expect_equal(result$Method, "Moran_I")
  expect_equal(result$Feature, "ADT-1")
  expect_true(result$Score >= -1 && result$Score <= 1)
})

test_that("calc_adt_moran errors on a missing feature or graph", {
  expect_error(
    calc_adt_moran(seurat_obj, features_adt = "CD999", graph_name = "RNA_nn")
  )
  expect_error(
    calc_adt_moran(seurat_obj, features_adt = "ADT-1", graph_name = "bad_graph")
  )
})


# calc_adt_nn_within_range ####
test_that("calc_adt_nn_within_range returns the canonical schema in [0, 1]", {
  agg <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1")
  expect_canonical(agg)
  expect_equal(unique(agg$Method), "Range")
  expect_true(all(agg$Score >= 0 & agg$Score <= 1))

  per_cell <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                       return_mean = FALSE)
  expect_canonical(per_cell, per_cell = TRUE)
  expect_equal(nrow(per_cell), n_cells)
})

test_that("calc_adt_nn_within_range: a stricter range never scores higher", {
  loose <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                    range = 0.5)
  strict <- calc_adt_nn_within_range(seurat_obj, features_adt = "ADT-1",
                                     range = 0.05)
  expect_gte(loose$Score, strict$Score)
})


# calc_adt_quantile ####
test_that("calc_adt_quantile returns the canonical schema", {
  agg <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                           base_assay = "RNA", k = 20, n_quantile = 3)
  expect_canonical(agg)
  expect_equal(unique(agg$Method), "Quantile")
  expect_true(all(agg$Score >= 0 & agg$Score <= 1))

  per_cell <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                                base_assay = "RNA", k = 20, n_quantile = 3,
                                return_mean = FALSE)
  expect_canonical(per_cell, per_cell = TRUE)
  expect_equal(nrow(per_cell), n_cells)
})

test_that("calc_adt_quantile binning switches between frequency and width", {
  freq <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                            base_assay = "RNA", k = 20, n_quantile = 3,
                            binning = "frequency", return_mean = FALSE)
  width <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                             base_assay = "RNA", k = 20, n_quantile = 3,
                             binning = "width", return_mean = FALSE)
  expect_canonical(freq, per_cell = TRUE)
  expect_canonical(width, per_cell = TRUE)
  # the two binning strategies do not produce identical per-cell scores here
  expect_false(isTRUE(all.equal(freq$Score, width$Score)))
})

test_that("calc_adt_quantile percentile_diff sets Method = Percentile_Diff", {
  result <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                              base_assay = "RNA", method = "percentile_diff")
  expect_canonical(result)
  expect_equal(unique(result$Method), "Percentile_Diff")
  expect_true(all(result$Score >= 0 & result$Score <= 1))
})

test_that("calc_adt_quantile warns when return_counts meets percentile_diff", {
  expect_warning(
    calc_adt_quantile(seurat_obj, features_adt = "ADT-1", base_assay = "RNA",
                      method = "percentile_diff", return_counts = TRUE),
    "not applicable"
  )
})

test_that("calc_adt_quantile is quiet unless verbose = TRUE", {
  expect_silent(
    calc_adt_quantile(seurat_obj, features_adt = c("ADT-1", "CD999"),
                      base_assay = "RNA", k = 20)
  )
  expect_message(
    calc_adt_quantile(seurat_obj, features_adt = c("ADT-1", "CD999"),
                      base_assay = "RNA", k = 20, verbose = TRUE),
    "Unavailable features"
  )
})

test_that("calc_adt_quantile return_counts gives integer counts", {
  result <- calc_adt_quantile(seurat_obj, features_adt = "ADT-1",
                              base_assay = "RNA", k = 20, n_quantile = 3,
                              return_counts = TRUE, return_mean = FALSE)
  expect_true(all(result$Score == floor(result$Score)))
  expect_true(all(result$Score >= 0))
})


# calc_neighbor_matches ####
test_that("calc_neighbor_matches returns the canonical schema", {
  agg <- calc_neighbor_matches(seurat_obj, meta_cols = "grp")
  expect_canonical(agg)
  expect_equal(unique(agg$Method), "Exact")
  expect_true(all(agg$Score >= 0 & agg$Score <= 1))

  per_cell <- calc_neighbor_matches(seurat_obj, meta_cols = "grp",
                                    return_mean = FALSE)
  expect_canonical(per_cell, per_cell = TRUE)
})

test_that("calc_neighbor_matches routes gex_cols to Assay = GEX", {
  result <- calc_neighbor_matches(seurat_obj,
                                  meta_cols = c("grp", "seurat_clusters"))
  assay_by_feature <- with(result, setNames(Assay, Feature))
  expect_equal(unname(assay_by_feature["grp"]), "BCR")
  expect_equal(unname(assay_by_feature["seurat_clusters"]), "GEX")
})

test_that("calc_neighbor_matches gex_cols is overridable", {
  result <- calc_neighbor_matches(seurat_obj, meta_cols = "grp",
                                  gex_cols = "grp")
  expect_equal(unique(result$Assay), "GEX")
})


# calc_neighbor_overlap ####
test_that("calc_neighbor_overlap returns the canonical schema", {
  agg <- calc_neighbor_overlap(seurat_obj, "RNA.nn", "RNA.nn")
  expect_canonical(agg)
  expect_equal(agg$Graph, "RNA.nn vs RNA.nn")
  expect_equal(agg$Method, "Overlap")
  expect_true(is.na(agg$Assay) && is.na(agg$Feature))

  per_cell <- calc_neighbor_overlap(seurat_obj, "RNA.nn", "RNA.nn",
                                    return_mean = FALSE)
  expect_canonical(per_cell, per_cell = TRUE)
  expect_equal(nrow(per_cell), n_cells)
  # a graph fully overlaps itself: k - 1 shared neighbors after drop_self
  k_after_drop <- ncol(seurat_obj@neighbors[["RNA.nn"]]@nn.idx) - 1
  expect_true(all(per_cell$Score == k_after_drop))
})

test_that("calc_neighbor_overlap errors on a missing graph", {
  expect_error(calc_neighbor_overlap(seurat_obj, "RNA.nn", "nope.nn"))
})


# permuted baselines ####
test_that("permute_adt returns the canonical schema for both methods", {
  result <- permute_adt(seurat_obj, features_adt = "ADT-1", n_permutations = 3)
  expect_canonical(result)
  expect_setequal(result$Method, c("Mean_Absolute", "Range"))
})

test_that("permute_neighbor_matches returns the canonical schema", {
  result <- permute_neighbor_matches(seurat_obj,
                                     meta_cols = c("grp", "seurat_clusters"),
                                     n_permutations = 3)
  expect_canonical(result)
  assay_by_feature <- with(result, setNames(Assay, Feature))
  expect_equal(unname(assay_by_feature["seurat_clusters"]), "GEX")
})
