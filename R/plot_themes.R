#' ggplot2 theme with horizontal grid lines only
#'
#' @export
theme_bw_custom <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey75",
                                               linewidth = 0.2)
  )

#' Standard label sizes for athanor plots
#'
#' @export
labels_standard <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
  plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5),
  plot.caption = ggplot2::element_text(size = 8),
  axis.title = ggplot2::element_text(size = 9),
  axis.text = ggplot2::element_text(size = 8),
  legend.title = ggplot2::element_text(size = 10),
  legend.text = ggplot2::element_text(size = 8)
)

#' Standard label sizes for violin plots
#'
#' @details
#' Removes "Identity".
#'
#' @export
labels_standard_vln <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
  plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
  axis.title.x = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_text(size = 10),
  axis.text = ggplot2::element_text(size = 10),
  legend.position = "none"
)

#' Violin plot labels with horizontal x-axis text
#'
#' @export
labels_standard_vln_rotate <- labels_standard_vln +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

#' Rotate x-axis labels 45 degrees
#'
#' @export
labels_rotate_x <- ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
)

#' Clean theme for DimPlots
#'
#' @export
clean_dimplot <- ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                axis.text.x = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(),
                                axis.ticks = ggplot2::element_blank(),
                                aspect.ratio = 1)

#' Clean theme for DimPlots
#'
#' @details
#' Same as clean_dimplot but with a square aspect ratio.
#'
#' @export
clean_dimplot2 <- clean_dimplot + ggplot2::theme(aspect.ratio = 1)

#' Patchwork annotation with Roman numeral panel labels
#'
#' @export
plot_anno <- patchwork::plot_annotation(tag_levels = "I")


#' Named color palettes for athanor plots
#'
#' A list of named character vectors mapping categorical values to hexadecimal colors, typically used as the `clrs_specific` argument in plotting functions. These color mappings are used to ensure consistent and meaningful coloring across different plots in the athanor package, especially for categories that are commonly visualized (e.g. cell types, isotypes, SHM frequencies).
#'
#' @details
#' The `named_colors` list contains color mappings for various categories relevant to BCR and single-cell analysis, including:
#' - `c_call`: Colors for IG constant region calls (e.g. IGHA, IGHD, etc.)
#' - `cdr3`: Colors for CDR3 lengths from 4 to 41
#' - `cell_types_celltypist`, `cell_types_simpler`: Colors for cell type annotations from CellTypist and a simpler scheme
#' - `datatype`: Colors for different data types (e.g. ADT, BCR, GEX, etc.)
#' - `doublet`: Colors for singlets vs doublets in doublet detection
#' - `embeddings`: Colors for different embedding methods (e.g. antiberta2, balm-paired, etc.)
#' - `isotype_stage`: Colors for isotype switching stages (i.e. Unswitched, Switched)
#' - `isotype`: Colors for immunoglobulin isotypes (e.g. IgA, IgD, etc.)
#' - `light`: Colors for light chain types (i.e. IGK, IGL, IGK/IGL)
#' - `mu_freq_bins`, `mu_freq_bins_fewer`, `mu_freq_bins_binary`: Colors for binned mutation frequencies
#' - `mu_freq_iso`: Colors for combinations of mutation frequency bins and isotype switching stages
#' - `v_call_family`, `d_call_family`, `j_call_family`: Colors for V, D, J gene families respectively
#' - `weights`: Colors for WNN weights from 0 to 1
#'
#' When possible, colors are chosen to be colorblind-friendly.
#'
#'
#' @export
named_colors <- list()

# general colors
named_colors$datatype <- c("ADT" = "#00853C", "BCR" = "#00826e",
                           "GEX" = "#38d77f", "TCR" = "#004036",
                           "WNN" = "#00c7d9", "GEX_BCR" = "#2276a3") # shades of green/teal

named_colors$weight_assay <- c("RNA" = named_colors$datatype[["GEX"]],
                               "BCR" = named_colors$datatype[["BCR"]],
                               "Tie" = named_colors$datatype[["GEX_BCR"]])

named_colors$isotype <- c("IgA" = alakazam::IG_COLORS[["IGHA"]], # blue
                          "IgD" = alakazam::IG_COLORS[["IGHD"]], # orange
                          "IgE" = alakazam::IG_COLORS[["IGHE"]], # red
                          "IgG" = alakazam::IG_COLORS[["IGHG"]], # green
                          "IgM" = alakazam::IG_COLORS[["IGHM"]]) # purple
# based on IG_COLORS from alakazam (+ base colors blended with black)
# doesn't have everything (e.g. there are 7 IGLCs)
named_colors$c_call <- c("IGHA" = alakazam::IG_COLORS[["IGHA"]],
                         "IGHA1" = "#25547B", "IGHA2" = "#122A3D",
                         "IGHD" = alakazam::IG_COLORS[["IGHD"]],
                         "IGHE" = alakazam::IG_COLORS[["IGHE"]],
                         "IGHG" = alakazam::IG_COLORS[["IGHG"]],
                         "IGHG1" = "#40923E", "IGHG2" = "#337531",
                         "IGHG3" = "#275825", "IGHG4" = "#1A3A19",
                         "IGHM" = alakazam::IG_COLORS[["IGHM"]],
                         "IGHK" = alakazam::IG_COLORS[["IGHK"]],
                         "IGKC" = alakazam::IG_COLORS[["IGHK"]],
                         "IGHL" = alakazam::IG_COLORS[["IGHL"]],
                         "IGL" = alakazam::IG_COLORS[["IGHL"]],
                         "IGLC1" = "#ffd208", "IGLC2" = "#e0b700")
named_colors$isotype_stage <-
  c("Unswitched" = "#aaaaff", "Switched" = "#51518a")
named_colors$light <- c("IGK" = "#61b498", "IGL" = pals::ocean.dense(n = 4)[2],
                        "IGK, IGL" = "#fcce00")

# gene family colors
named_colors$v_call_igh_family <-
  setNames(c("#7690c7", "#8d4bca", "#90b648", "#d16099",
             "#57865f", "#62396e", "#bc9149"),
           nm = str_c("IGHV", 1:7)) # no IGHV8
named_colors$v_call_igk_family <-
  setNames(pals::ocean.solar(n = 7), nm = str_c("IGKV", 1:7))
named_colors$v_call_igl_family <-
  setNames(pals::ocean.haline(n = 10), nm = str_c("IGLV", 1:10))
named_colors$v_call_family <-
  c(named_colors$v_call_igh_family, named_colors$v_call_igk_family,
    named_colors$v_call_igl_family)
named_colors$d_call_family <-
  setNames(pals::ocean.matter(n = 7), nm = str_c("IGHD", 1:7))
named_colors$j_call_igh_family <-
  setNames(pals::ocean.amp(n = 6), nm = str_c("IGHJ", 1:6))
named_colors$j_call_igk_family <-
  setNames(pals::ocean.algae(n = 5), nm = str_c("IGKJ", 1:5))
named_colors$j_call_igl_family <-
  setNames(pals::ocean.dense(n = 7), nm = str_c("IGLJ", 1:7))
named_colors$j_call_family <-
  c(named_colors$j_call_igh_family, named_colors$j_call_igk_family,
    named_colors$j_call_igl_family)

# for doublet plots
named_colors$doublet <- c("singlet" = "#cbcae3", "doublet" = "#807DBA")

# for plotting cell types
# TODO: make sure these don't conflict
named_colors$cell_types_celltypist <-
  c("Age-associated B cells" = "#a83ab8",
    "C1 non-classical monocytes" = "#6b4c93",
    "CD16- NK cells" = "#4fb3a9",
    "CD16+ NK cells" = "#32c9a6",
    "Classical monocytes" = "#6b1650",
    "Cycling immune mix" = "#05b58c",
    "DC1" = "#f4a261",
    "DC2" = "#E39C37",
    "gdT" = "#99e9ff",
    "HSC/MPP" = "#e76f51",
    "MAIT cells" = "#0b607c",
    "Megakaryocytes/platelets" = "#f07fb0",
    "Memory B cells" = "#d95f4b",
    "Naive B cells" = "#ffd3b6",
    "Non-classical monocytes" = "#8f1468",
    "Non-switched memory B cells" = "#0b8546",
    "pDC" = "#bf7e21",
    "Plasma cells" = "#6b2342",
    "Regulatory T cells" = "#264653",
    "Switched memory B cells" = "#054d27",
    # these are similar but it doesn't matter since we're focused on B cells
    "Tcm/Naive cytotoxic T cells" = "#5cc9ed",
    "Tcm/Naive helper T cells" = "#4e7380",
    "Tem/Effector helper T cells" = "#11a5d6",
    "Tem/Temra cytotoxic T cells" = "#118ab2",
    "Tem/Temra helper T cells" = "#0e6986",
    "Tem/Trm cytotoxic T cells" = "#0b607c")
# blue shifted to see them better
named_colors$cell_types_simpler <-
  c("Naive B cells" = "#63CC9E", "Non-Naive B cells" = "#076659")

# for the WNN weights
named_colors$weights <- setNames(pals::ocean.algae(n = 11), nm = seq(0, 1, 0.1))

# for the binned mutation frequencies (gives FFs for full alphas at the ends)
# can't set the names because the last one will vary depending on the max freq
named_colors$mu_freq_bins <-
  scales::viridis_pal(direction = -1, option = "C")(n = 5) %>% stringr::str_sub(1, 7)
named_colors$mu_freq_bins_fewer <-
  setNames(scales::viridis_pal(direction = -1, option = "C")(n = 3) %>% stringr::str_sub(1, 7),
           nm = c("0%", "0% to 1%", ">1%"))
named_colors$mu_freq_bins_fewer[["0% to 1%"]] <- named_colors$mu_freq_bins[[2]]
named_colors$mu_freq_bins_binary <- c("0% to 1%" = "#F4C731", ">1%" = "#A52590")

# manually crossed between the existing color vectors (mu_freq and iso_stage)
named_colors$mu_freq_iso <- c("0% Unswitched" = "#CDD290",
                              "0% Switched" = "#A1A556",
                              "0% to 1% Unswitched" = "#D19FA0",
                              "0% to 1% Switched" = "#A57366",
                              ">1% Unswitched" = "#5C59C3",
                              ">1% Switched" = "#2F2D89")

# adjust the lengths as needed
named_colors$cdr3 <- setNames(pals::viridis(length(4:41)), nm = 4:41)

# embedding methods
# https://coolors.co/palette/003049-6b2c39-d62828-f77f00-fcbf49-eae2b7
named_colors$embeddings <-
  setNames(c("#003049", "#6B2C39", "#D62828", "#F77F00", "#FCBF49", "#EAE2B7"),
           nm = c("antiberta2", "antiberty", "balm-paired", "esm2",
                  "immune2vec", "simulated"))

# TODO: rearrange in alphabetical order
