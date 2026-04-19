# make Seurat objects for testing [similar to sim_gex_seurat_obj()]

library(Matrix)
library(Seurat)

# for consistency
set.seed(42)

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
saveRDS(seurat_obj,
        file.path(testthat::test_path(), "fixtures", "manual.rds"))



# make a small object with GEX and ADT assays, 10 cells, neighbors calculated
num_genes <- 100
num_cells <- 10
num_proteins <- 5

raw_counts <- as.integer(rexp(num_genes * num_cells, rate = 0.5))
raw_counts <- Matrix(data = raw_counts,
                     nrow = num_genes, ncol = num_cells,
                     sparse = TRUE)

adt_counts <- sample(1:5, num_proteins * num_cells, replace = TRUE)
adt_counts <- Matrix(data = adt_counts,
                     nrow = num_proteins, ncol = num_cells,
                     sparse = TRUE)

rownames(raw_counts) <- paste("Gene", 1:num_genes, sep = "-")
colnames(raw_counts) <- paste("Cell", 1:num_cells, sep = "-")
rownames(adt_counts) <- paste("ADT", 1:num_proteins, sep = "-")
colnames(adt_counts) <- paste("Cell", 1:num_cells, sep = "-")

seurat_obj2 <- CreateSeuratObject(counts = raw_counts, project = "Manual2")
seurat_obj2[["ADT"]] <- CreateAssayObject(counts = adt_counts)
seurat_obj2$cell_id <- Cells(seurat_obj2)

seurat_obj2 <- NormalizeData(seurat_obj2, verbose = FALSE)
seurat_obj2 <- NormalizeData(seurat_obj2, normalization.method = "CLR",
                             margin = 2, assay = "ADT", verbose = FALSE)
seurat_obj2 <- FindVariableFeatures(seurat_obj2, nfeatures = 50,
                                    verbose = FALSE)
seurat_obj2 <- ScaleData(seurat_obj2, verbose = FALSE)
# with only 10 cells, npcs must be < 10
seurat_obj2 <- RunPCA(seurat_obj2, npcs = 5, verbose = FALSE)
## Warning in svd.function(A = t(x = object), nv = npcs, ...) :
##   You're computing too large a percentage of total singular values, use a standard svd instead.
seurat_obj2 <- FindNeighbors(seurat_obj2, reduction = "pca", dims = 1:5,
                             k.param = 5,
                             graph.name = str_c("RNA_", c("", "s"), "nn"),
                             verbose = FALSE)
seurat_obj2 <- FindNeighbors(seurat_obj2, reduction = "pca", dims = 1:5,
                             k.param = 5, return.neighbor = TRUE,
                             graph.name = "RNA.nn", verbose = FALSE)

seurat_obj2@neighbors$RNA.nn@nn.idx
## [,1] [,2] [,3] [,4] [,5]
## [1,]    1    2    8    5    6
## [2,]    2    1    3    4    7
## [3,]    3    9    2    6    4
## [4,]    4    8    5   10    6
## [5,]    5   10    9    1    8
## [6,]    6    8    3    1    4
## [7,]    7    1    8    4    9
## [8,]    8    4    1    5    6
## [9,]    9    5   10    3    6
#3 [10,]   10    5    9    4    3

seurat_obj2 <- RunUMAP(seurat_obj2, nn.name = "RNA.nn")

saveRDS(seurat_obj2,
        file.path(testthat::test_path(), "fixtures", "manual2.rds"))
