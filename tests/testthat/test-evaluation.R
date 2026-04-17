seurat_obj <- readRDS(file.path(testthat::test_path(), "fixtures", "manual.rds"))

# correlation ####
test_that("Correlation doesn't work on non-existent features", {
  expect_error(calc_correlation(seurat_obj, features_adt = "CD19"),
               regexp = "The requested ADT feature is not present in assay")
})

test_that("Need the neighbors slot to not be empty", {
  expect_error(calc_correlation(seurat_obj, features_adt = "ADT-1"),
               regexp = "No neighbor graphs found in object.")
})

# fill in the neighbors slot
seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "pca",
                            assay = "RNA", k.param = 20,
                            return.neighbor = TRUE, graph.name = "RNA.nn")
