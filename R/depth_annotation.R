# Per-cell sequencing depth as a row annotation on the numbat panels (github #43).
#
# Why a bar and not another distribution plot
# -------------------------------------------
# plot_cell_depth_panel() already shows each sample's depth distribution, but a
# distribution cannot say WHICH cells are shallow. The finding that motivates
# this is per-cell, not per-sample: spearman(nCount_gene, p_opt) is positive in
# 36/36 samples (median +0.319), so cells in the shallow tail carry clone labels
# numbat itself is unsure of. Reading that off a heatmap needs depth aligned to
# the same rows as the clone calls.
#
# Fixed bins, not quantiles
# -------------------------
# Within-sample quantiles would make every sample's bar look identical and hide
# the 23-fold spread across the cohort (median nCount_gene 831 to 18,980). The
# bin edges below are absolute, so bars are comparable between samples, and the
# lowest bin is exactly the complement of the pipeline's own filter
# (`nCount_gene > 1000` in filter_sample_qc() and generate_filtering_cell_counts()) --
# cells numbat used but the Seurat analysis later discards.

#' Bin per-cell counts into fixed, cross-sample-comparable depth classes
#'
#' @param n Numeric vector of per-cell counts (`nCount_gene`).
#' @return An ordered factor with levels `<=1k`, `1-2.5k`, `2.5-5k`, `5-10k`, `>10k`.
#' @export
depth_bins <- function(n) {
  cut(as.numeric(n),
      breaks = c(-Inf, 1000, 2500, 5000, 10000, Inf),
      labels = c("<=1k", "1-2.5k", "2.5-5k", "5-10k", ">10k"),
      right  = TRUE, ordered_result = TRUE)
}

#' Palette for [depth_bins()]
#'
#' The lowest bin is red because it is the failing side of the pipeline's own
#' `nCount_gene > 1000` filter; the rest is a light-to-dark blue ramp so deeper
#' reads darker.
#'
#' @return Named character vector of hex colours.
#' @export
depth_bin_pal <- function() {
  c("<=1k"   = "#bd4034",
    "1-2.5k" = "#cfe0ee",
    "2.5-5k" = "#8fb8dc",
    "5-10k"  = "#4a86bd",
    ">10k"   = "#1f4e79")
}

#' Attach a depth bin column to a numbat annotation frame
#'
#' numbat's `plot_phylo_heatmap()` melts every non-`cell` column of `annot` into
#' one value column with `value.factor = TRUE`, so a numeric depth would be
#' coerced to character and ordered lexicographically. The column therefore has
#' to be discrete before it goes in.
#'
#' @param annot Data frame with a `cell` column (e.g. cell + clone_opt).
#' @param counts Named numeric vector of per-cell counts, names = cell barcodes.
#' @return `annot` with a `depth` column appended, or unchanged when no cell matches.
#' @export
add_depth_annot <- function(annot, counts) {
  if (is.null(counts) || length(counts) == 0) return(annot)
  n <- counts[as.character(annot$cell)]
  if (all(is.na(n))) return(annot)
  annot$depth <- as.character(depth_bins(n))
  annot
}

#' Combined palette for a clone bar plus a depth bar
#'
#' `annot_bar()` applies ONE `scale_fill_manual()` to every annotation bar, so
#' both bars have to be served by a single palette. Names disambiguate them:
#' clone ids are "1", "2", ... and depth bins are "<=1k", "1-2.5k", ...
#'
#' @param clone_pal Named character vector of clone colours.
#' @return Named character vector covering clones and depth bins.
#' @export
clone_depth_pal <- function(clone_pal) c(clone_pal, depth_bin_pal())


# ---------------------------------------------------------------------------
# Expression panel
# ---------------------------------------------------------------------------

#' numbat's rolling-expression panel with a sequencing-depth bar
#'
#' A fork of `numbat:::plot_exp_roll()`, which takes no annotation argument. Two
#' changes beyond the bar itself:
#'
#'  * **The cell sample is seeded.** numbat draws `n_sample` cells with an
#'    unseeded `sample()`, so the panel differs on every rebuild and no bar could
#'    be aligned to it. That also means the current panel is not reproducible
#'    run to run, independently of this annotation.
#'  * **The bar shares the heatmap's row faceting.** `facet_grid(cluster ~ .)`
#'    with `space = "free"` gives both panels row heights proportional to the
#'    same cell counts, which is what keeps them aligned.
#'
#' @param gexp_roll_wide Cells x genes rolling-mean expression matrix.
#' @param hc Hierarchical clustering of cells (`nb$hc`).
#' @param k Number of clusters to cut `hc` into.
#' @param gtf Gene position table (`nb$gtf`).
#' @param lim Colour scale limit, passed through to the fill scale.
#' @param n_sample Number of cells to draw.
#' @param annot Optional data frame with `cell` and `depth` columns; when NULL
#'   the function reproduces numbat's own panel (seeded).
#' @param seed Integer seed for the cell draw.
#' @param plot_tree Whether to prepend the dendrogram.
#' @param annot_width Relative width of the depth bar.
#' @return A patchwork object.
#' @export
plot_exp_roll_annotated <- function(gexp_roll_wide, hc, k, gtf, lim = 0.8,
                                    n_sample = 300, annot = NULL, seed = 1L,
                                    plot_tree = TRUE, annot_width = 0.45) {

  gexp_norm_long <- gexp_roll_wide %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    data.table::as.data.table() %>%
    data.table::melt(id.var = "cell", variable.name = "gene",
                     value.name = "exp_rollmean") %>%
    dplyr::left_join(gtf, by = "gene") %>%
    dplyr::arrange(CHROM, gene_start) %>%
    dplyr::mutate(gene_index = as.integer(factor(gene, unique(gene))))

  cells <- unique(gexp_norm_long$cell)
  # Seed locally and put the stream back, so callers' RNG is untouched.
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  set.seed(seed)
  cell_sample <- sample(cells, min(n_sample, length(cells)), replace = FALSE)
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())

  p_tree <- ggtree::ggtree(hc, size = 0.2)
  cell_order <- as.data.frame(p_tree$data) %>%
    dplyr::filter(isTip) %>%
    dplyr::arrange(y) %>%
    dplyr::pull(label)
  clusters <- stats::cutree(hc, k = k)

  chrom_labeller <- function(chr) { chr[chr %in% c(18, 21)] <- ""; chr }

  d <- gexp_norm_long %>%
    dplyr::filter(cell %in% cell_sample) %>%
    dplyr::mutate(cell = factor(cell, cell_order)) %>%
    dplyr::mutate(cluster = clusters[as.character(cell)]) %>%
    dplyr::arrange(cell) %>%
    dplyr::mutate(cluster = factor(cluster, rev(unique(cluster))))

  p_heatmap <- ggplot2::ggplot(d) +
    ggplot2::geom_raster(ggplot2::aes(x = gene_index, y = cell, fill = exp_rollmean)) +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  midpoint = 0, limits = c(-lim, lim),
                                  oob = scales::squish) +
    ggplot2::theme_void() +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      legend.position = "top",
      panel.spacing = ggplot2::unit(0, "mm"),
      panel.border = ggplot2::element_rect(linewidth = 0.5, color = "wheat4", fill = NA),
      strip.text.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(cluster ~ CHROM, scales = "free", space = "free",
                        labeller = ggplot2::labeller(CHROM = chrom_labeller)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Expression\nmagnitude")) +
    ggplot2::xlab("Gene index") + ggplot2::ylab("Cell")

  if (is.null(annot) || !"depth" %in% colnames(annot)) {
    return(if (plot_tree) (p_tree | p_heatmap) + patchwork::plot_layout(widths = c(1, 10))
           else p_heatmap)
  }

  # One row per drawn cell, carrying the same factor levels and the same row
  # facet as the heatmap -- that is what makes the two line up.
  bar_d <- d %>%
    dplyr::distinct(cell, cluster) %>%
    dplyr::left_join(dplyr::select(annot, cell, depth) %>%
                       dplyr::mutate(cell = as.character(cell)),
                     by = c("cell" = "cell")) %>%
    dplyr::mutate(depth = factor(depth, names(depth_bin_pal())))

  p_annot <- ggplot2::ggplot(bar_d) +
    ggplot2::geom_raster(ggplot2::aes(x = 1L, y = cell, fill = depth)) +
    ggplot2::scale_fill_manual(values = depth_bin_pal(), na.value = "gray90",
                               drop = FALSE, name = "Depth\n(nCount_gene)") +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      legend.position = "top",
      panel.spacing = ggplot2::unit(0, "mm"),
      panel.border = ggplot2::element_rect(linewidth = 0.5, color = "wheat4", fill = NA),
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(cluster ~ ., scales = "free", space = "free")

  if (plot_tree) {
    (p_tree | p_annot | p_heatmap) +
      patchwork::plot_layout(widths = c(1, annot_width, 10), guides = "collect")
  } else {
    (p_annot | p_heatmap) +
      patchwork::plot_layout(widths = c(annot_width, 10), guides = "collect")
  }
}
