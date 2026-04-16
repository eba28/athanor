# make a Seurat object for testing [similar to sim_gex_seurat_obj()]

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
seurat_obj <-
  CreateSeuratObject(counts = raw_counts, project = "Manual")
seurat_obj[["ADT"]] <-
  CreateAssayObject(counts = adt_counts, assay = "ADT", project = "Manual")

# add metadata
seurat_obj$cell_id <- Cells(seurat_obj)

# run the Seurat pipeline
seurat_obj <- seurat_pipeline(seurat_obj, nfeatures_RNA = 0, perc_mt = 100,
                              num_features = 1000, num_pcs = 30, num_dims = 20)

# fill in the misc slot
Misc(seurat_obj, slot = "category") <- "Manual"

# save the object
saveRDS(seurat_obj, file.path(testthat::test_path(), "fixtures", "manual.rds"))



# make an object where the ADT features are not random

# set up the counts matrices
num_genes <- 100
num_cells <- 10
num_proteins <- 5
raw_counts <- as.integer(rexp(num_genes * num_cells, rate = 0.5))
raw_counts <- Matrix(data = raw_counts,
                     nrow = num_genes, ncol = num_cells,
                     sparse = TRUE)
adt_counts <- sample(c(0,1), num_proteins * num_cells, replace = TRUE)
adt_counts <- Matrix(data = adt_counts,
                     nrow = num_proteins, ncol = num_cells,
                     sparse = TRUE)

# format the genes and cells
rownames(raw_counts) <- paste("Gene", 1:num_genes, sep = "-")
colnames(raw_counts) <- paste("Cell", 1:num_cells, sep = "-")
rownames(adt_counts) <- paste("ADT", 1:num_proteins, sep = "-")
colnames(adt_counts) <- paste("Cell", 1:num_cells, sep = "-")

# check the ADT expression
data.frame(adt_counts)
## Cell.1 Cell.2 Cell.3 Cell.4 Cell.5 Cell.6 Cell.7 Cell.8 Cell.9 Cell.10
## ADT-1      0      1      1      1      1      1      0      1      0       0
## ADT-2      0      1      0      1      0      0      0      0      0       1
## ADT-3      0      0      1      1      0      0      0      0      1       0
## ADT-4      1      1      0      0      0      1      0      0      0       0
## ADT-5      0      1      1      0      0      0      1      0      0       0

# create the object
seurat_obj <-
  CreateSeuratObject(counts = raw_counts, project = "Manual")
seurat_obj[["ADT"]] <-
  CreateAssayObject(counts = adt_counts, assay = "ADT", project = "Manual")

# add metadata
seurat_obj$cell_id <- Cells(seurat_obj)

# run the Seurat pipeline
# seurat_obj <- seurat_pipeline(seurat_obj, nfeatures_RNA = 0, perc_mt = 100,
#                               num_features = 100, num_pcs = 5, num_dims = 5)
# Error in svd.function(A = t(x = object), nv = npcs, ...) :
# max(nu, nv) must be strictly less than min(nrow(A), ncol(A))

# save the object
# saveRDS(seurat_obj, file.path(testthat::test_path(), "fixtures", "manual2.rds"))
