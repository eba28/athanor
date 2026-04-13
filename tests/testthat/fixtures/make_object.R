# make a Seurat object for testing [similar to sim_gex_manual()]

library(Matrix)
library(Seurat)

# set up the counts matrices
num_genes <- 1000
num_cells <- 100
num_proteins <- 5
raw_counts <- as.integer(rexp(num_genes * num_cells, rate = 0.5))
raw_counts <- Matrix(data = raw_counts,
                     nrow = num_genes, ncol = num_cells,
                     sparse = TRUE)
adt_counts <- as.integer(rexp(num_proteins * num_cells, rate = 0.5))
adt_counts <- Matrix(data = adt_counts,
                     nrow = num_proteins, ncol = num_cells,
                     sparse = TRUE)

# format the genes and cells
rownames(raw_counts) <- paste("Gene", 1:num_genes, sep = "-")
colnames(raw_counts) <- paste("Cell", 1:num_cells, sep = "-")
rownames(adt_counts) <- paste("ADT", 1:num_proteins, sep = "-")
colnames(adt_counts) <- paste("Cell", 1:num_cells, sep = "-")

# create the object
manual <- CreateSeuratObject(counts = raw_counts, project = "Manual")
manual[["ADT"]] <- CreateAssayObject(counts = adt_counts, assay = "ADT", project = "Manual")

# add metadata
manual$cell_id <- Cells(manual)

# run the Seurat pipeline
manual <- seurat_pipeline(manual, nfeatures_RNA = 0, perc_mt = 100,
                          num_features = 1000, num_pcs = 30, num_dims = 20)

# save the object
saveRDS(manual, file.path(testthat::test_path(), "fixtures", "manual.rds"))
