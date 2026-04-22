# minimal stubs for package-global ggplot2 theme objects used in plotting.R
# these are normally defined by the user in their analysis environment

library(ggplot2)

theme_bw <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey75",
                                               linewidth = 0.2)
  )

labels_standard <- theme(plot.title = element_text(size = 12, hjust = 0.5),
                         plot.subtitle = element_text(size = 9, hjust = 0.5),
                         axis.title = element_text(size = 9),
                         axis.text = element_text(size = 8),
                         legend.title = element_text(size = 10),
                         legend.text = element_text(size = 8))

labels_standard_vln <- theme(plot.title = element_text(size = 14, hjust = 0.5),
                              axis.title.x = element_blank(),
                              axis.title.y = element_text(size = 10),
                              axis.text = element_text(size = 10),
                              legend.position = "none")

labels_standard_vln_rotate <- labels_standard_vln +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

labels_rotate_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

clean_umap <- theme(plot.title = element_text(hjust = 0.5),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks = element_blank(),
                    aspect.ratio = 1)

plot_anno <- patchwork::plot_annotation(tag_levels = "I")
