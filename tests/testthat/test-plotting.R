seurat_obj <- readRDS(testthat::test_path("fixtures", "manual.rds"))

# calc_pcts ####
test_that("calc_pcts returns a data frame with expected columns", {
  data <- data.frame(
    sample_id = rep(c("S1", "S2"), each = 6),
    Dataset = rep(c("D1", "D2"), each = 6),
    Cell_Type = rep(c("B", "T", "NK"), times = 4)
  )

  result <- calc_pcts(data, meta_group_by = "sample_id",
                      focus_group = "Cell_Type")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample_id", "Cell_Type", "Count", "Percent") %in%
                    colnames(result)))
})

test_that("calc_pcts percentages sum to 100 per group", {
  data <- data.frame(
    sample_id = rep(c("S1", "S2"), each = 9),
    Cell_Type = rep(c("B", "T", "NK"), times = 6)
  )

  result <- calc_pcts(data, meta_group_by = "sample_id",
                      focus_group = "Cell_Type")

  sums <- result %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(total = sum(Percent))

  expect_equal(sums$total, c(100, 100))
})

test_that("calc_pcts fills missing group combinations with 0", {
  # S1 has B and T, S2 has only NK
  data <- data.frame(
    sample_id = c("S1", "S1", "S2"),
    Cell_Type = c("B", "T", "NK")
  )

  result <- calc_pcts(data, meta_group_by = "sample_id",
                      focus_group = "Cell_Type")

  # all 3 cell types should appear for both samples
  expect_equal(nrow(result), 6)

  missing_row <- result[result$sample_id == "S2" &
                          result$Cell_Type == "B", ]
  expect_equal(missing_row$Count, 0)
  expect_equal(missing_row$Percent, 0)
})

test_that("calc_pcts order_by reorders the grouping variable", {
  data <- data.frame(
    sample_id = rep(c("S1", "S2", "S3"), each = 3),
    Cell_Type = rep(c("B", "T", "NK"), times = 3)
  )
  data$Count_extra <- c(3, 0, 0, 1, 1, 0, 0, 0, 3)

  result <- calc_pcts(data, meta_group_by = "sample_id",
                      focus_group = "Cell_Type", order_by = "B")

  expect_s3_class(result$sample_id, "factor")
  # S1 has most B cells so it should come first
  expect_equal(levels(result$sample_id)[1], "S1")
})

test_that("calc_pcts handles a single group", {
  data <- data.frame(
    sample_id = rep("S1", 4),
    Cell_Type = c("B", "B", "T", "T")
  )

  result <- calc_pcts(data, meta_group_by = "sample_id",
                      focus_group = "Cell_Type")

  expect_equal(sum(result$Percent), 100)
  expect_equal(nrow(result), 2)
})


# plot_color_scale ####
test_that("plot_color_scale returns a character vector when data is provided", {
  values <- c(-3, -1, 0, 1, 3)
  result <- plot_color_scale(data = values)

  expect_type(result, "character")
  expect_gt(length(result), 0)
})

test_that("plot_color_scale errors when both plot and data are provided", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  p$data[["avg.exp.scaled"]] <- 1

  expect_error(plot_color_scale(plot = p, data = c(1, 2, 3)))
})

test_that("plot_color_scale errors when neither plot nor data are provided", {
  expect_error(plot_color_scale())
})

test_that("plot_color_scale returns more colors when range is wider", {
  narrow <- plot_color_scale(data = c(-1, 0, 1))
  wide <- plot_color_scale(data = c(-10, 0, 10))

  # both should return valid color vectors
  expect_type(narrow, "character")
  expect_type(wide, "character")
})

test_that("plot_color_scale handles positive-only values", {
  result <- plot_color_scale(data = c(0, 1, 2, 3))

  expect_type(result, "character")
  expect_gt(length(result), 0)
})

test_that("plot_color_scale handles negative-only values", {
  result <- plot_color_scale(data = c(-3, -2, -1, 0))

  expect_type(result, "character")
  expect_gt(length(result), 0)
})


# add_info_bar ####
test_that("add_info_bar returns a ggplot object", {
  df <- data.frame(x = letters[1:4], y = 1:4,
                   Group = c("A", "A", "B", "B"))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point()

  result <- add_info_bar(p, method = "add", info_type = "Group",
                         info = df["Group"])

  expect_s3_class(result, "gg")
})

test_that("add_info_bar adds column facets by default (top_side = TRUE)", {
  df <- data.frame(x = letters[1:4], y = 1:4,
                   Group = c("A", "A", "B", "B"))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point()

  result <- add_info_bar(p, method = "add", info_type = "Group",
                         info = df["Group"])

  # facet_grid creates a FacetGrid
  expect_s3_class(result$facet, "FacetGrid")
})

test_that("add_info_bar adds row facets when top_side = FALSE", {
  df <- data.frame(x = letters[1:4], y = 1:4,
                   Group = c("A", "A", "B", "B"))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point()

  result <- add_info_bar(p, method = "add", info_type = "Group",
                         info = df["Group"], top_side = FALSE)

  expect_s3_class(result$facet, "FacetGrid")
})


# plot_pcts ####
test_that("plot_pcts returns a ggplot object", {
  data <- data.frame(
    sample_id = rep(c("S1", "S2"), each = 3),
    Cell_Type = rep(c("B", "T", "NK"), times = 2)
  )
  pcts <- calc_pcts(data, meta_group_by = "sample_id",
                    focus_group = "Cell_Type")

  result <- plot_pcts(pcts, data_source = "Blood",
                      fill_type = "Cell_Type",
                      x_axis = "sample_id")

  expect_s3_class(result, "gg")
})

test_that("plot_pcts title includes data_source and plot_value", {
  data <- data.frame(
    sample_id = rep(c("S1", "S2"), each = 3),
    Cell_Type = rep(c("B", "T", "NK"), times = 2)
  )
  pcts <- calc_pcts(data, meta_group_by = "sample_id",
                    focus_group = "Cell_Type")

  result <- plot_pcts(pcts, data_source = "Blood",
                      plot_value = "Cell Type",
                      fill_type = "Cell_Type",
                      x_axis = "sample_id")

  expect_true(grepl("Blood", result$labels$title))
  expect_true(grepl("Cell Type", result$labels$title))
})

test_that("plot_pcts errors for invalid plot_type", {
  data <- data.frame(
    sample_id = rep("S1", 3),
    Cell_Type = c("B", "T", "NK")
  )
  pcts <- calc_pcts(data, meta_group_by = "sample_id",
                    focus_group = "Cell_Type")

  expect_error(
    plot_pcts(pcts, data_source = "Blood", fill_type = "Cell_Type",
              x_axis = "sample_id", plot_type = "Invalid")
  )
})

test_that("plot_pcts drop_zeroes removes rows with Percent == 0", {
  # S2 is missing NK so it will have Percent = 0 for NK after calc_pcts
  data <- data.frame(
    sample_id = c("S1", "S1", "S1", "S2", "S2"),
    Cell_Type = c("B", "T", "NK", "B", "T")
  )
  pcts <- calc_pcts(data, meta_group_by = "sample_id",
                    focus_group = "Cell_Type")

  result_drop <- plot_pcts(pcts, data_source = "Blood",
                           fill_type = "Cell_Type", x_axis = "sample_id",
                           drop_zeroes = TRUE)
  result_keep <- plot_pcts(pcts, data_source = "Blood",
                           fill_type = "Cell_Type", x_axis = "sample_id",
                           drop_zeroes = FALSE)

  expect_lt(nrow(result_drop$data), nrow(result_keep$data))
})


# plot_dimplot ####
test_that("plot_dimplot returns a ggplot object", {
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "RNA",
                         reduc = "umap")

  expect_s3_class(result, "gg")
})

test_that("plot_dimplot sets the plot title to title", {
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "GEX",
                         reduc = "umap")

  expect_equal(result$labels$title, "GEX")
})

test_that("plot_dimplot sets subtitle to data_source", {
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "RNA",
                         data_source = "TestData", reduc = "umap")

  expect_equal(result$labels$subtitle, "TestData")
})

test_that("plot_dimplot details overrides the subtitle", {
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "RNA",
                         data_source = "TestData",
                         details = "Custom subtitle", reduc = "umap")

  expect_equal(result$labels$subtitle, "Custom subtitle")
})

test_that("plot_dimplot respects include_legend = FALSE", {
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "RNA",
                         include_legend = FALSE, reduc = "umap")

  # NoLegend() sets legend.position to "none"
  expect_equal(result$theme$legend.position, "none")
})

test_that("plot_dimplot highlights specific clusters via meta_col", {
  first_cluster <- as.character(levels(seurat_obj$seurat_clusters)[1])
  clrs <- setNames("red", first_cluster)
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "RNA",
                         reduc = "umap",
                         meta_col = "seurat_clusters",
                         highlight = first_cluster,
                         clrs_specific = clrs)

  expect_s3_class(result, "gg")
})

test_that("plot_dimplot applies legend_label with highlight", {
  first_cluster <- as.character(levels(seurat_obj$seurat_clusters)[1])
  clrs <- setNames("red", first_cluster)
  result <- plot_dimplot(seurat_obj = seurat_obj, plot_title = "RNA",
                         reduc = "umap",
                         meta_col = "seurat_clusters",
                         highlight = first_cluster,
                         clrs_specific = clrs,
                         legend_label = "Data Type")

  expect_equal(result$labels$colour, "Data Type")
})


# plot_vln_feat ####
test_that("plot_vln_feat returns a patchwork object", {
  feature <- rownames(seurat_obj)[1]
  result <- plot_vln_feat(seurat_obj = seurat_obj, feature = feature)

  expect_s3_class(result, "patchwork")
})

test_that("plot_vln_feat works with a specified group_col", {
  feature <- rownames(seurat_obj)[1]
  result <- plot_vln_feat(seurat_obj = seurat_obj, feature = feature,
                          group_col = "seurat_clusters")

  expect_s3_class(result, "patchwork")
})
