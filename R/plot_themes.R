# internal ggplot2 theme constants used across plotting functions

# TODO: properly incorporate these into the package
# TODO: give this a better name so it doesn't look like the base theme is just missing parentheses
theme_bw <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey75",
                                               linewidth = 0.2)
  )

labels_standard <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
  plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5),
  plot.caption = ggplot2::element_text(size = 8),
  axis.title = ggplot2::element_text(size = 9),
  axis.text = ggplot2::element_text(size = 8),
  legend.title = ggplot2::element_text(size = 10),
  legend.text = ggplot2::element_text(size = 8)
)

labels_standard_vln <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
  plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
  axis.title.x = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_text(size = 10),
  axis.text = ggplot2::element_text(size = 10),
  legend.position = "none"
)

labels_standard_vln_rotate <- labels_standard_vln +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

labels_rotate_x <- ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
)

clean_umap <- ggplot2::theme(
  plot.title = ggplot2::element_text(hjust = 0.5),
  axis.text.x = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  aspect.ratio = 1
)

plot_anno <- patchwork::plot_annotation(tag_levels = "I")
