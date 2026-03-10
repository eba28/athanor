#' Display a nicely formatted table
#'
#' @description
#' This function creates a formatted table with scrolling capability for use in
#' R Markdown documents. It applies striped styling and makes the table scrollable
#' within specified dimensions.
#'
#' @param table The input data frame to be printed
#' @param kable_height The height of the output table (you can set it to NULL to
#'   display the full table without scrolling). Default is "500px".
#' @param kable_width The width of the output table. Default is "100%".
#'
#' @returns A formatted table.
#' @export
print_kable <- function(table, kable_height = "500px", kable_width = "100%") {
  if (nrow(table) > 0) {
    kable(table) %>%
      kable_styling("striped") %>%
      scroll_box(height = kable_height, width = kable_width)
  }
}


#' Display a sortable, scrollable `DataTable` table
#'
#' @param table The input data frame to be printed.
#' @param dt_width The width of the output table.
#'
#' @returns A DataTable table.
#' @export
print_dt <- function(table, dt_width = "800px") {
  DT::datatable(data = table, options = list(scrollX = TRUE),
                rownames = FALSE, filter = "top", width = dt_width)
}


#' Reduce a Seurat object's size
#'
#' @description
#' This function reduces a Seurat object by removing count matrices and keeping only
#' specified dimensionality reductions. This is especially useful for creating
#' lightweight objects for Shiny apps or sharing.
#'
#' @details
#' Modify this as needed if your object is built differently (e.g. tSNE instead).
#' Change the annotation column name if needed.
#' Uses DietSeurat to remove counts while preserving reductions and metadata.
#' This is especially useful if you are making a Shiny app.
#'
#' @param seurat_obj Processed Seurat object to reduce.
#' @param dim_reducs Vector of dimensionality reductions to keep.
#' @param print_size Whether to print info about how much the object was reduced.
#' @param load_annotations Whether to load and add cell type annotations.
#' @param annotations_file File path to CSV file containing cluster and cell type annotations.
#'
#' @returns A reduced Seurat object with specified reductions kept.
#' @export
reduce_object <- function(seurat_obj, dim_reducs = "umap", print_size = TRUE,
                          load_annotations = FALSE, annotations_file) {
  cat(paste("Currently reducing:", deparse(substitute(seurat_obj)), "\n"))

  # modify this as desired
  obj_reduced <- DietSeurat(seurat_obj, layers = NULL, dimreducs = dim_reducs)

  # print the before and after sizes
  if (print_size) {
    cat(paste("Original object size:",
              format(object.size(seurat_obj), units = "auto"), "\n"))
    cat(paste("Reduced object size:",
              format(object.size(obj_reduced), units = "auto"), "\n"))
  }

  # add in the annotations (which should be in cluster order already) if desired
  if (load_annotations) {
    # read in the data - you only need the cell type information
    annotations <- read_csv(annotations_file,
                            col_types = cols(.default = "c"))$CellType

    # add the annotations to the object
    # assumes the cluster and annotations columns are seurat_clusters and annotated_clusters
    names(annotations) <- levels(obj_reduced)
    obj_reduced <- RenameIdents(obj_reduced, annotations)
    Idents(obj_reduced) <- factor(Idents(obj_reduced),
                                  levels = sort(levels(obj_reduced))) # alphabetize the cell types
    obj_reduced$annotated_clusters <- Idents(obj_reduced) # useful metadata
  }

  return(obj_reduced)
}
