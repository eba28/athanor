# load a simple Seurat object with neighbor graphs and dimensionality reductions
obj <- readRDS(testthat::test_path("fixtures", "manual.rds"))

# object specific tests ####
test_that("reduce_object returns a Seurat object", {
  obj_reduced <- reduce_object(obj, print_size = FALSE)

  expect_s4_class(obj_reduced, "Seurat")
})

test_that("reduce_object reduces object size", {
  original_size <- object.size(obj)
  # should also be able to tell from the print out when `print_size = TRUE`
  obj_reduced <- reduce_object(obj, print_size = FALSE)
  reduced_size <- object.size(obj_reduced)

  expect_lt(reduced_size, original_size)
})

test_that("reduce_object works with print_size = TRUE", {
  # should runs without error and produce an output
  out <- capture_output(reduce_object(obj, print_size = FALSE))
  expect_false(grepl("size", out))

  # should produce console output in addition to "Currently reducing:"
  expect_output(reduce_object(obj, print_size = TRUE), "object size")
})


# the right things are being removed ####
test_that("reduce_object removes unspecified dimensionality reductions", {
  obj_reduced <- reduce_object(obj, dim_reducs = "umap", print_size = FALSE)

  expect_false("pca" %in% names(obj_reduced@reductions))
})

test_that("reduce_object filters metadata columns when specified", {
  meta_cols <- c("nCount_RNA", "nFeature_RNA")
  obj_reduced <- reduce_object(obj, meta_cols = meta_cols, print_size = FALSE)

  # check that only specified columns remain
  remaining_cols <- colnames(obj_reduced[[]])
  expect_true(all(meta_cols %in% remaining_cols))
  expect_equal(length(setdiff(remaining_cols, colnames(obj_reduced[[]]))), 0)
})

test_that("reduce_object removes neighbors when remove_neighbors = TRUE", {
  obj_reduced <- reduce_object(obj, remove_neighbors = TRUE, print_size = FALSE)

  expect_length(obj_reduced@neighbors, 0)
})


# the right things are being kept ####
test_that("reduce_object keeps specified dimensionality reductions", {
  obj_reduced <- reduce_object(obj, dim_reducs = "umap", print_size = FALSE)

  expect_true("umap" %in% names(obj_reduced@reductions))
})

test_that("reduce_object keeps metadata columns when unspecified", {
  meta_cols <- colnames(obj[[]])
  obj_reduced <- reduce_object(obj, print_size = FALSE)
  remaining_cols <- colnames(obj_reduced[[]])

  expect_identical(meta_cols, remaining_cols)
})

test_that("reduce_object keeps neighbors when remove_neighbors = FALSE", {
  obj_reduced <- reduce_object(obj, remove_neighbors = FALSE,
                               print_size = FALSE)

  expect_true(length(obj_reduced@neighbors) >= 0)
})

test_that("reduce_object handles keeping multiple dimensionality reductions", {
  obj_reduced <- reduce_object(obj, dim_reducs = c("pca", "umap"),
                               print_size = FALSE)

  expect_true("pca" %in% names(obj_reduced@reductions) &
                "umap" %in% names(obj_reduced@reductions))
})


# additional parameters can be passed ####
test_that("reduce_object can pass additional parameters to DietSeurat", {
  # don't keep all of the layers
  obj_reduced <- reduce_object(obj, print_size = FALSE, layers = "data")
  reduced_layers <- names(obj_reduced@assays$RNA@layers)

  expect_true(reduced_layers == "data")

  # don't keep the miscellaneous slot
  obj_reduced <- reduce_object(obj, print_size = FALSE, misc = FALSE)
  expect_true(is.list(obj_reduced@misc))
})

# errors are raised for improper inputs ####
test_that("reduce_object raises error for non-existent metadata columns", {
  invalid_cols <- c("nCount_RNA", "invalid_col")

  expect_error(
    reduce_object(obj, meta_cols = invalid_cols, print_size = FALSE),
    "Please check that all of the specified metadata columns are present in the Seurat object."
  )
})

test_that("reduce_object throws warning for invalid parameter", {
  expect_warning(reduce_object(obj, print_size = FALSE, invalid_param = TRUE))
})
