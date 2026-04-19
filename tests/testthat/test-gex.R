load_manual_object <- function() {
	readRDS(testthat::test_path("fixtures", "manual.rds"))
}


# add_annotations ####
test_that("add_annotations relabels clusters and adds metadata", {
  seurat_obj <- load_manual_object()
  annotations_df <- data.frame(CellType = "B Cells")

  annotated_obj <- add_annotations(seurat_obj, annotations_df,
                                   cell_types_col = "CellType",
                                   relabel = TRUE, relocate = TRUE,
                                   alphabetize = TRUE,
                                   clusters_col = "seurat_clusters",
                                   annotations_col = "annotated_clusters")

  expect_identical(as.character(Seurat::Idents(annotated_obj)[1]), "B Cells")
  expect_true("annotated_clusters" %in% colnames(annotated_obj[[]]))

  expect_identical(as.character(unique(annotated_obj$annotated_clusters)),
                   "B Cells")
  expect_equal(levels(annotated_obj$annotated_clusters), "B Cells")

  expect_gt(which(colnames(annotated_obj[[]]) == "annotated_clusters"),
            which(colnames(annotated_obj[[]]) == "seurat_clusters"))
})


# automated_annotation ####
test_that("automated_annotation rejects unsupported methods", {
  expect_error(
    automated_annotation(seurat_obj = NA, annotation_method = "Other"),
    "Method must be one of"
  )
})


# cell_type_clusters ####
test_that("cell_type_clusters maps cluster IDs to the dominant annotation", {
  seurat_obj <- load_manual_object()
  seurat_obj$seurat_clusters <- factor("0")
  annotations_df <- data.frame(CellType = "B Cells")

  annotated_obj <- add_annotations(seurat_obj, annotations_df,
                                   cell_types_col = "CellType",
                                   relabel = TRUE, relocate = TRUE,
                                   alphabetize = TRUE,
                                   clusters_col = "seurat_clusters",
                                   annotations_col = "annotated_clusters")

  cluster_map <- cell_type_clusters(annotated_obj,
                                    clusters_col = "seurat_clusters",
                                    annotations_col = "annotated_clusters")

  # TODO: check that seurat_clusters has levels
  expect_equal(nrow(cluster_map), 1)
  expect_equal(cluster_map$annotated_clusters, factor("B Cells"))
  expect_equal(cluster_map$Clusters, "0")
})


test_that("cell_type_clusters summarizes mappings accurately", {
  df <- data.frame(seurat_clusters = c(0, 0, 0, 1, 1),
                   predicted_labels = c(rep("T Cells", 2), rep("B Cells", 3)))
  mock_obj <- list(df)
  class(mock_obj) <- "Seurat"

  # since the function uses seurat_obj[[]], we need to copy that behavior
  # or use a real object:
  # obj <- load_manual_object()
  obj <- CreateSeuratObject(counts = matrix(0, nrow = 1, ncol = 5))
  obj@meta.data <- df

  sum_df <- cell_type_clusters(obj, annotations_col = "predicted_labels")

  expect_s3_class(sum_df, "data.frame")
  expect_equal(sum_df$Clusters[sum_df$predicted_labels == "T Cells"], "0")
  expect_equal(sum_df$Clusters[sum_df$predicted_labels == "B Cells"], "1")
})


# get_airr_genes ####
test_that("get_airr_genes returns the Ensembl version and unique genes", {
	skip_if_not_installed("biomaRt")

	biomaRt_ns <- asNamespace("biomaRt")
	mocked_names <- c("listEnsemblArchives", "useEnsembl", "searchDatasets",
	                  "listFilters", "listAttributes", "getBM")

	original_values <- lapply(mocked_names, function(name) {
		get(name, envir = biomaRt_ns, inherits = FALSE)
	})

	for (name in mocked_names) {
		unlockBinding(name, biomaRt_ns)
	}

	on.exit({
		for (i in seq_along(mocked_names)) {
			name <- mocked_names[[i]]
			unlockBinding(name, biomaRt_ns)
			assign(name, original_values[[i]], envir = biomaRt_ns)
			lockBinding(name, biomaRt_ns)
		}
	}, add = TRUE)

	assign("listEnsemblArchives",
			 function() {
				 data.frame(name = "Ensembl 114", current_release = "*",
								stringsAsFactors = FALSE)
			 },
			 envir = biomaRt_ns)
	assign("useEnsembl",
			 function(...) {
				 structure(list(), class = "Mart")
			 },
			 envir = biomaRt_ns)
	assign("searchDatasets",
			 function(...) {
				 data.frame(dataset = "hsapiens_gene_ensembl",
								stringsAsFactors = FALSE)
			 },
			 envir = biomaRt_ns)
	assign("listFilters",
			 function(...) {
				 data.frame(name = "biotype")
			 },
			 envir = biomaRt_ns)
	assign("listAttributes",
			 function(...) {
				 data.frame(name = c("ensembl_gene_id", "external_gene_name",
											"gene_biotype", "hgnc_symbol",
											"description"),
								page = rep("feature_page", 5),
								stringsAsFactors = FALSE)
			 },
			 envir = biomaRt_ns)
	assign("getBM",
			 function(...) {
				 data.frame(
					 ensembl_gene_id = c("id1", "id2", "id3", "id4"),
					 external_gene_name = c("IGHV1-1", "", "TRAV1-1", "IGHV1-1"),
					 gene_biotype = c("IG_V_gene", "IG_V_gene", "TR_V_gene",
											"IG_V_gene"),
					 hgnc_symbol = c("IGHV1-1", "", "TRAV1-1", "IGHV1-1"),
					 description = c("immunoglobulin heavy variable 1-1",
										  "", "T cell receptor alpha variable 1-1",
										  "immunoglobulin heavy variable 1-1"),
					 stringsAsFactors = FALSE
				 )
			 },
			 envir = biomaRt_ns)

	result <- get_airr_genes(genome = "hsapiens", category = c("IG", "TR"))

	expect_equal(result$ensembl_version, "114")
	expect_equal(result$remove_genes, c("IGHV1-1", "TRAV1-1"))
})

test_that("get_airr_genes returns correct structure", {
  # Mocking biomaRt calls to avoid network dependency in CI/CD
  stub(get_airr_genes, "biomaRt::listEnsemblArchives",
       data.frame(name = "Ensembl 114", current_release = "*"))
  stub(get_airr_genes, "biomaRt::getBM",
       data.frame(external_gene_name = c("TRBV1", "IGHV1"),
                  gene_biotype = "IG_V_gene"))

  # Run function
  result <- get_airr_genes(genome = "hsapiens", category = "IG")

  expect_type(result, "list")
  expect_named(result, c("ensembl_version", "remove_genes"))
  expect_equal(result$ensembl_version, "114")
  expect_true("TRBV1" %in% result$remove_genes)
})


# find_k_clusters ####
test_that("find_k_clusters finds the target cluster count", {
	seurat_obj <- load_manual_object()

	clustered_obj <- find_k_clusters(seurat_obj, graph_name = "RNA_snn",
												desired_k = 1)

	expect_s4_class(clustered_obj, "Seurat")
	expect_equal(dplyr::n_distinct(clustered_obj$seurat_clusters), 1)
})

test_that("find_k_clusters stops when desired_k is exceeded", {
  obj <- CreateSeuratObject(counts = matrix(rpois(500, 1), ncol = 50))
  # Standard preprocessing needed for FindClusters to work
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, verbose = FALSE)

  # Try to find 100 clusters in a 50-cell dataset (should fail)
  expect_error(find_k_clusters(obj, desired_k = 100),
               "Could not find resolution to match desired clusters.")

  # Test exceeding (mocking a low k)
  expect_error(find_k_clusters(obj, desired_k = 1),
               "The number of desired clusters has been exceeded.")
})


# seurat_pipeline ####
test_that("seurat_pipeline runs end to end on the Seurat fixture", {
  seurat_obj <- load_manual_object()

  processed_obj <- NULL
  expect_warning(
    processed_obj <- seurat_pipeline(seurat_obj, nfeatures_RNA = 0,
                                     perc_mt = 100, num_features = 100,
                                     num_pcs = 10, num_dims = 10,
                                     cluster_res = 0.4, verbose = FALSE),
    "No filtration was performed upon this object."
  )

  expect_s4_class(processed_obj, "Seurat")
  expect_true("rpca" %in% names(processed_obj@reductions))
  expect_true("rna.umap" %in% names(processed_obj@reductions))
  expect_true("RNA.nn" %in% names(processed_obj@neighbors))
  expect_true("RNA_snn_res.0.4" %in% colnames(processed_obj[[]]))
})

test_that("seurat_pipeline skips cell filtering when nfeatures_RNA and perc_mt are omitted", {
  seurat_obj <- load_manual_object()
  n_before <- ncol(seurat_obj)

  processed_obj <- seurat_pipeline(seurat_obj, num_features = 100,
                                   num_pcs = 10, num_dims = 10, verbose = FALSE)

  expect_equal(ncol(processed_obj), n_before)
})

test_that("seurat_pipeline processes object and filters genes", {
  counts <- matrix(rpois(1000, 5), ncol = 10)
  rownames(counts) <- paste0("GENE-", 1:100)
  colnames(counts) <- paste0("CELL-", 1:10)
  obj <- CreateSeuratObject(counts = counts)

  stub(seurat_pipeline, "get_airr_genes",
       list(ensembl_version = "114", remove_genes = c("GENE-1", "GENE-2")))

  processed_obj <- seurat_pipeline(obj, num_features = 20, num_pcs = 5,
                                   num_dims = 5, filter_genes = "IG",
                                   verbose = FALSE)

  expect_s4_class(processed_obj, "Seurat")
  expect_true("rpca" %in% names(processed_obj@reductions))
  expect_true("rna.umap" %in% names(processed_obj@reductions))
  expect_false(any(c("GENE-1", "GENE-2") %in% VariableFeatures(processed_obj)))
})


# filter_variable_features ####
test_that("filter_variable_features removes genes and saves Ensembl version", {
  seurat_obj <- load_manual_object()
  ig_gene <- "IGHV1-1"
  VariableFeatures(seurat_obj) <- c(VariableFeatures(seurat_obj)[1:50], ig_gene)

  stub(filter_variable_features, "get_airr_genes",
       list(ensembl_version = "114", remove_genes = ig_gene))

  result <- filter_variable_features(seurat_obj, filter_genes = "IG")

  expect_false(ig_gene %in% VariableFeatures(result))
  expect_equal(Misc(result, slot = "ensembl_version"), "114")
})

test_that("filter_variable_features reports GEX-only count when bcr_features provided", {
  seurat_obj <- load_manual_object()
  ig_gene <- "IGHV1-1"
  VariableFeatures(seurat_obj) <- c(VariableFeatures(seurat_obj)[1:10], ig_gene)

  bcr_features <- matrix(0, nrow = 3, ncol = ncol(seurat_obj),
                         dimnames = list(c("bcr-dim-1", "bcr-dim-2", "bcr-dim-3"),
                                         Cells(seurat_obj)))

  stub(filter_variable_features, "get_airr_genes",
       list(ensembl_version = "114", remove_genes = ig_gene))

  expect_output(
    filter_variable_features(seurat_obj, filter_genes = "IG",
                             bcr_features = bcr_features),
    "GEX-only"
  )
})
