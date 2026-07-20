# Builds the small real-data fixture used by vignettes/real_data_pipeline.Rmd.
#
# Source: GSE270142 (Nicholas et al., diabetes), subject 108_ND (a healthy,
# non-diabetic control), available locally at
# /media/edelaron/Midas/diabetes/nicholas/. This script subsamples ~1500 GEX
# cells (preserving the natural fraction that have paired BCR data) so the
# fixture is small enough to bundle in the package, then writes the result to
# inst/extdata/diabetes_nicholas_108_ND/. Not run as part of R CMD check;
# rerun by hand if the fixture needs to be regenerated.


# I'm storing the data in vignettes/ for now instead of making it part of the package itself

library(dplyr)
library(readr)
library(stringr)
library(Matrix)

set.seed(42)

path_raw <- "/media/edelaron/Midas/diabetes/nicholas"
# path_out <- "inst/extdata/diabetes_nicholas_108_ND"
path_out <- "vignettes/data/diabetes_nicholas_108_ND"
sample_id <- "108_ND"
n_cells <- 1500 # TODO: try different cell counts

# raw 10x counts (RNA + ADT) for this subject (non-standard GEO file names,
# so read the trio manually instead of via Seurat::Read10X())
path_sample <- file.path(path_raw, "default", "108_ND")
mat <- Matrix::readMM(file.path(path_sample, "GSE270142_108_matrix.mtx.gz"))
features <- read_tsv(file.path(path_sample, "GSE270142_108_features.tsv.gz"),
                     col_names = c("id", "symbol", "feature_type"),
                     show_col_types = FALSE)
barcodes <- read_tsv(file.path(path_sample, "GSE270142_108_barcodes.tsv.gz"),
                     col_names = "barcode", show_col_types = FALSE)$barcode
rownames(mat) <- make.unique(features$symbol)
colnames(mat) <- barcodes
mat <- as(mat, "CsparseMatrix")

is_rna <- features$feature_type == "Gene Expression"
is_adt <- features$feature_type == "Antibody Capture"
counts <- list(`Gene Expression` = mat[is_rna, ],
              `Antibody Capture` = mat[is_adt, ])

# subsample barcodes (this naturally preserves the real fraction with BCR data)
barcodes_sub <- sort(sample(barcodes, size = n_cells))
cell_id_sub <- paste0(sample_id, "_", str_remove(barcodes_sub, "-1$"))

rna_counts <- counts$`Gene Expression`[, barcodes_sub]
adt_counts <- counts$`Antibody Capture`[, barcodes_sub]
colnames(rna_counts) <- cell_id_sub
colnames(adt_counts) <- cell_id_sub

saveRDS(rna_counts, file.path(path_out, "rna_counts.rds"))
saveRDS(adt_counts, file.path(path_out, "adt_counts.rds"))

# real AntiBERTa2 BCR embeddings for this subject (heavy + light concatenated)
embeddings_raw <-
  read_tsv(file.path(path_raw, "airrflow", "bcr", "5.0.0", "results",
                     "embeddings", "antiberta2", "108_ND.tsv"),
           show_col_types = FALSE) %>%
  mutate(cell_id = paste0(sample_id, "_", str_extract(cell_id, "^[ACGT]{16}"))) %>%
  filter(cell_id %in% cell_id_sub)

embeddings <- t(as.matrix(select(embeddings_raw, -cell_id)))
colnames(embeddings) <- embeddings_raw$cell_id
# Seurat rewrites "_" to "-" in feature names at assay creation, so avoid
# underscores in the dimension names up front (see seurat_pipeline()'s note
# on this in R/general.R)
rownames(embeddings) <- paste0("dim.", seq_len(nrow(embeddings)))

saveRDS(embeddings, file.path(path_out, "bcr_embeddings.rds"))

# real AIRR metadata (mu_freq, isotype, cdr3 features, v/j calls, ...) for
# this subject, subset from the full cohort table
combined_airr <-
  read_tsv(file.path(path_raw, "embeddings", "2.1.2",
                     "combined_airrflow_bcr.tsv"),
           show_col_types = FALSE) %>%
  filter(subject_id == 108, cell_id %in% cell_id_sub)

saveRDS(combined_airr, file.path(path_out, "combined_airr.rds"))
