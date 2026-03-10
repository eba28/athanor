#' Calculate the proportion of a cell's neighbors that are mutated or unmutated
#'
#' @description
#' This function calculates the proportion of neighbors that are mutated or unmutated for each cell in a Seurat object.
#' It uses the neighbor information stored in the Seurat object to determine the mutation status of neighboring cells and computes the fraction accordingly.
#' The function also generates a plot to visualize the distribution of mutated or unmutated neighbors across different cell types.
#'
#' @details
#' Assumes that WNN neighbors are stored in "w.nn".
#' Defines mutated as anything about 0% SHM.
#' `setNames()` works better than a for loop checking all of the neighbors.
#'
#' @param seurat_obj A Seurat object with neighbors calculated.
#' @param assay One of: RNA, BCR, WNN
#' @param category One of: Mutated, Unmutated
#' @param plot_title A string to add to the plot title, e.g. "All Cells" or "Memory B Cells Only".
#'
#' @returns A data.frame and a plot.
#' @export
calc_nn_frac <- function(seurat_obj, assay = "WNN", plot_title = "",
                         category = "mutated") {
  embedding <- seurat_obj@misc$embedding_type

  # use Seurat's neighbor calculations
  # each column is a k, each row is a cell
  if (assay == "WNN") {
    neighbors <- seurat_obj@neighbors$w.nn@nn.idx
    rownames(neighbors) <- colnames(seurat_obj)
  }

  # more setup
  k <- ncol(neighbors) # usually 20 but could be different
  new_col <- paste("nn_frac", tolower(category), assay, sep = "_")

  if (category == "Mutated") {
    # only keep the mutated cells regardless of label
    mutation_status <-
      setNames(seurat_obj[[]]$mu_freq > 0, seurat_obj[[]]$cell_id)

    nn_frac <-
      seurat_obj[[]] %>%
      select(cell_id, mu_freq) %>%
      filter(mu_freq > 0) %>%
      rowwise() %>%
      mutate(cell_nn = list(colnames(seurat_obj)[neighbors[cell_id, ]]),
             !!sym(new_col) := sum(mutation_status[cell_nn]) / k) %>%
      ungroup() %>%
      select(-mu_freq, -cell_nn)
  } else {
    # only keep the unmutated cells regardless of label
    mutation_status <-
      setNames(seurat_obj[[]]$mu_freq == 0, seurat_obj[[]]$cell_id)

    nn_frac <-
      seurat_obj[[]] %>%
      select(cell_id, mu_freq) %>%
      filter(mu_freq == 0) %>%
      rowwise() %>%
      mutate(cell_nn = list(colnames(seurat_obj)[neighbors[cell_id, ]]),
             !!sym(new_col) := sum(mutation_status[cell_nn]) / k) %>%
      ungroup() %>%
      select(-mu_freq, -cell_nn)
  }

  # plot the data
  p <- ggplot(left_join(nn_frac,
                        seurat_obj[[]] %>% select(cell_id, annotated_clusters),
                        by = join_by(cell_id)),
              aes(x = annotated_clusters, y = !!sym(new_col),
                  color = annotated_clusters)) +
    geom_jitter(width = 0.5, height = 0.01, size = 0.6) +
    geom_hline(yintercept = 0.5, linewidth = 0.2, color = "black") +
    labs(title = paste(plot_title, category, "Cells with", category,
                       "Neighbors"),
         subtitle = paste(embedding, assay, paste("k =", k),
                          paste0(nrow(nn_frac), " cells (", ncol(seurat_obj),
                                 " total cells)"),
                          sep = " | "),
         x = "Cell Type", y = paste(category, "Cell Neighbors Percent")) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 1/k)) +
    scale_color_manual(values = seurat_obj@misc$colors_annotated) +
    theme_bw + labels_standard +
    theme(legend.position = "none", panel.grid.major.y = element_blank(),
          # bit more visible
          panel.grid.minor.y =
            element_line(linewidth = 0.2, color = "gray85"))

  return(list(df = nn_frac, plot = p))
}
