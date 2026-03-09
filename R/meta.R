#' Give an overview of the data in terms of files, samples and datasets.
#'
#' @description
#' This function provides a summary of metadata files, printing information about the number of subjects, samples, and total files for each disease in the dataset.
#' Can handle metadata files representing multiple diseases and datasets.
#'
#' @param meta_file The analyst-created metadata csv file containing disease,
#'   sample, and subject information.
#'
#' @returns A text description of the provided metadata file.
#' @export
print_metadata_summary <- function(meta_file) {
  for (disease in unique(meta_file$Disease)) {
    meta_disease <- meta_file %>% filter(Disease == disease)
    subjects_list <- paste(unique(meta_disease$SampleNameOriginal),
                           collapse = ", ")

    cat(paste("There are", dplyr::n_distinct(meta_disease$SampleNameOriginal),
              paste0("subjects (", subjects_list, ")"), "and",
              dplyr::n_distinct(meta_disease$SampleNameOriginal),
              "samples for a total of", nrow(meta_disease),
              "total files in the", disease, "data.\n"))
  }
}
