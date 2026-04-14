#' Reformat VDJ barcodes by adding sample names and removing suffixes
#'
#' @description
#' This function adds the sample name to the beginning of barcodes and removes
#' the "-1" suffixes from the ends to create unique cell identifiers.
#'
#' @param data The input dataset which contains a barcode column.
#' @param sample_name The sample name to add to make the barcode unique.
#' @param barcode_col The column in the data that contains the barcodes.
#'
#' @returns A character vector of reformatted barcodes.
#' @export
reformat_vdj_barcode_sample <- function(data, sample_name,
                                        barcode_col = "barcode") {
  unlist(lapply(strsplit(as.character(data[[barcode_col]]), "-"),
                function(x) paste0(sample_name, "_", x[1])))
}

#' Reformat VDJ barcodes to make them unique across samples
#'
#' @description
#' This function processes barcode data by separating cell IDs from suffixes and
#' creating unique cell identifiers by combining sample names with barcodes.
#'
#' @param data Input data frame containing barcode information.
#' @param col_samples Column name containing sample names.
#' @param col_barcodes Column name containing cell IDs/barcodes.
#' @param col_output Name for the output column.
#'
#' @returns A character vector of unique cell identifiers.
#' @export
reformat_vdj_barcode <- function(data, col_samples = "sample_id",
                                 col_barcodes = "cell_id",
                                 col_output = "cell_id") {
  data <- data %>%
    select(all_of(col_samples), all_of(col_barcodes)) %>%
    separate(!!sym(col_barcodes), sep = "-",
             into = "barcodes", remove = FALSE, extra = "drop") %>%
    unite(col_output, c(!!sym(col_samples), "barcodes"), sep = "_",
          remove = FALSE)

  return(data$col_output)
}