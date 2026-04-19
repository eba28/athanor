make_embeddings <- function(num_dims = 10, num_cells = 20, seed = 42) {
  set.seed(seed)
  m <- matrix(rnorm(num_dims * num_cells), nrow = num_dims, ncol = num_cells)
  rownames(m) <- paste0("Dim-", seq_len(num_dims))
  colnames(m) <- paste0("Cell-", seq_len(num_cells))
  m
}


# bcr_embeddings_pipeline ####
test_that("bcr_embeddings_pipeline returns a Seurat object with BCR assay", {
  embeddings <- make_embeddings()

  result <- bcr_embeddings_pipeline(embeddings, embedding_type = "test",
                                    num_pcs = 5, num_dims = 5, k_param = 5)

  expect_s4_class(result, "Seurat")
  expect_true("BCR" %in% names(result@assays))
  expect_equal(ncol(result), ncol(embeddings))
})

test_that("bcr_embeddings_pipeline produces expected reductions and graphs", {
  embeddings <- make_embeddings()

  result <- bcr_embeddings_pipeline(embeddings, embedding_type = "test",
                                    num_pcs = 5, num_dims = 5, k_param = 5)

  expect_true("bpca" %in% names(result@reductions))
  expect_true("bcr.umap" %in% names(result@reductions))
  expect_true("BCR.nn" %in% names(result@neighbors))
  expect_true("BCR_nn" %in% names(result@graphs))
  expect_true("BCR_snn" %in% names(result@graphs))
})

test_that("bcr_embeddings_pipeline stores Misc slots correctly", {
  embeddings <- make_embeddings(num_dims = 10)

  result <- bcr_embeddings_pipeline(embeddings, embedding_type = "immune2vec",
                                    num_pcs = 5, num_dims = 5, k_param = 5)

  expect_equal(Misc(result, slot = "embedding_type"), "immune2vec")
  expect_equal(Misc(result, slot = "embeddings_dims"), 10)
})

test_that("bcr_embeddings_pipeline sets cell_id to Cells()", {
  embeddings <- make_embeddings()

  result <- bcr_embeddings_pipeline(embeddings, embedding_type = "test",
                                    num_pcs = 5, num_dims = 5, k_param = 5)

  expect_true("cell_id" %in% colnames(result[[]]))
  expect_equal(result$cell_id, Cells(result))
})

test_that("bcr_embeddings_pipeline sets scale.data from embeddings directly", {
  embeddings <- make_embeddings()

  result <- bcr_embeddings_pipeline(embeddings, embedding_type = "test",
                                    num_pcs = 5, num_dims = 5, k_param = 5)

  scale_data <- GetAssayData(result, assay = "BCR", layer = "scale.data")
  expect_equal(as.matrix(scale_data), embeddings)
})
