# Per-bar enrichment test backing `signif = TRUE` in
# plot_distribution_of_clones_across_clusters().
#
# For each level of `var_y` (a cluster), build the 2xK table of `var_x` counts
# in that level vs every other cell, and test it with Fisher's exact. Returns
# one row per level -- always every level, including the untested ones, so a
# missing mark is never ambiguous:
#
#   p_value NA + label "... (n<min)" -> too few cells, no test attempted
#   p_value NA + label "... (1 grp)" -> only one var_x group present, untestable
#   p_value set                      -> tested; q_value is BH across tested levels
#
# BH is applied within the panel only (one panel = one sample x SCNA x
# resolution). That scope is stated in the plot caption; it is deliberately NOT
# corrected across panels, which would require knowing the whole sweep here.
.bar_enrichment_stats <- function(seu_meta, var_x, var_y, min_cells = 20) {
  x <- as.character(seu_meta[[var_x]])
  y <- as.character(seu_meta[[var_y]])
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]

  levs <- levels(seu_meta[[var_y]])
  levs <- if (is.null(levs)) sort(unique(y)) else levs[levs %in% y]
  levs <- setdiff(levs, "all")

  res <- purrr::map_dfr(levs, function(lv) {
    inside <- y == lv
    n <- sum(inside)
    tab <- table(factor(inside, levels = c(TRUE, FALSE)), x)
    # Drop var_x groups absent from the whole object; a structural zero column
    # carries no information and only slows the exact test down.
    tab <- tab[, colSums(tab) > 0, drop = FALSE]

    if (n < min_cells) {
      return(tibble::tibble(level = lv, n_cells = n, p_value = NA_real_,
                            note = paste0("n<", min_cells)))
    }
    if (ncol(tab) < 2) {
      return(tibble::tibble(level = lv, n_cells = n, p_value = NA_real_,
                            note = "1 grp"))
    }
    # 2x2 is exact and cheap; wider tables use Monte Carlo, which is why the
    # seed is fixed -- otherwise the same collage rebuilt twice yields different
    # stars.
    p <- tryCatch({
      if (ncol(tab) == 2) {
        stats::fisher.test(tab)$p.value
      } else {
        # Restore the caller's RNG state afterwards -- a collage build should not
        # silently reseed the session it runs in.
        old_seed <- if (exists(".Random.seed", .GlobalEnv)) {
          get(".Random.seed", .GlobalEnv)
        } else {
          NULL
        }
        on.exit({
          if (is.null(old_seed)) {
            suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
          } else {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
          }
        }, add = TRUE)
        set.seed(1L)
        stats::fisher.test(tab, simulate.p.value = TRUE, B = 1e4)$p.value
      }
    }, error = function(e) NA_real_)
    tibble::tibble(level = lv, n_cells = n, p_value = p, note = NA_character_)
  })

  if (nrow(res) == 0) {
    return(tibble::tibble(level = character(0), n_cells = integer(0),
                          p_value = numeric(0), q_value = numeric(0),
                          note = character(0), stars = character(0),
                          label = character(0)))
  }

  res$q_value <- stats::p.adjust(res$p_value, method = "BH")
  res$stars <- dplyr::case_when(
    is.na(res$q_value)  ~ NA_character_,
    res$q_value < 0.001 ~ "***",
    res$q_value < 0.01  ~ "**",
    res$q_value < 0.05  ~ "*",
    res$q_value < 0.1   ~ ".",
    TRUE                ~ "ns"
  )
  res$label <- ifelse(is.na(res$stars),
                      paste0(res$level, " (", res$note, ")"),
                      paste0(res$level, " ", res$stars))
  res
}

#' plot distribution of clones across clusters
#'
#' @param signif Annotate each `var_y` bar with the significance of its `var_x`
#'   composition against all other cells (Fisher's exact; Monte Carlo when more
#'   than two `var_x` groups; BH-adjusted across the tested bars in this panel
#'   only). The mark rides in the axis label, since `position = "fill"` bars
#'   leave no headroom. The full table is attached to the returned plot as
#'   `attr(., "enrichment_stats")`.
#' @param signif_min_cells Bars with fewer cells than this are not tested and are
#'   labelled `(n<min)`, so an absent star is never mistaken for a
#'   non-significant one (default 20). Ignored when `signif = FALSE`.
#' @export
plot_distribution_of_clones_across_clusters <- function(seu, seu_name, var_x = "scna", var_y = "SCT_snn_res.0.6", plot_type = c("both", "clone", "cluster"), avg_line = NULL, signif = FALSE, signif_min_cells = 20, integrated = FALSE, reverse_fill = FALSE) {
  plot_type <- match.arg(plot_type)

  # reverse_fill flips the stacking direction of the position="fill" bars so the
  # fill levels read in ascending order (e.g. clones 1,2,3 left-to-right after
  # coord_flip). Default FALSE preserves the original order for existing callers.
  fill_pos <- ggplot2::position_fill(reverse = reverse_fill)

  seu_meta <- seu@meta.data
  summarized_clones <- dplyr::mutate(dplyr::select(seu_meta, .data[[var_x]], .data[[var_y]]), scna = "all")
  cluster_plot <- ggplot(seu_meta) +
    geom_bar(position = fill_pos, aes(x = .data[[var_x]], fill = .data[[var_y]])) +
    geom_bar(data = summarized_clones, position = fill_pos, aes(x = .data[[var_x]], fill = .data[[var_y]])) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20), limits = rev) +
    coord_flip()

  summarized_clusters <- dplyr::mutate(dplyr::select(seu_meta, .data[[var_x]], .data[[var_y]]), clusters = "all")

  if (signif) {
    # Per bar (level of var_y, normally a cluster), test whether the composition
    # over var_x -- clone, or on the two-clone collages the SCNA-of-interest
    # status -- differs from the rest of the object: a 2xK contingency of
    # "this bar" vs "all other cells". Fisher's exact handles the small counts
    # the degenerate hypoxia/two-clone subsets routinely produce; K > 2 falls
    # back to a Monte Carlo p-value since the exact network algorithm blows up.
    #
    # Everything below keys on var_x/var_y rather than hardcoded "scna"/
    # "clusters", so the same marks work on the full multi-clone collages and on
    # the two-clone ones.
    stats_tbl <- .bar_enrichment_stats(seu_meta, var_x, var_y, min_cells = signif_min_cells)

    # Marks ride in the axis label rather than as text above the bar: these are
    # position="fill" stacked bars, so there is no free y-space above a bar, and
    # every bar is the same height. `all` is the pooled summary bar and is not a
    # test, so it never carries a mark.
    lab_map <- stats::setNames(stats_tbl$label, as.character(stats_tbl$level))
    relabel <- function(x) {
      x <- as.character(x)
      ifelse(x %in% names(lab_map), unname(lab_map[x]), x)
    }
    y_levels <- c("all", relabel(stats_tbl$level))

    # Own pooled-summary frame keyed on var_y. The shared `summarized_clusters`
    # above hardcodes a `clusters` column, which only lands on var_y when the
    # caller happens to pass var_y = "clusters".
    summary_bar <- seu_meta %>%
      dplyr::select(.data[[var_x]], .data[[var_y]]) %>%
      dplyr::mutate(!!var_y := "all")

    clone_input <- seu_meta %>%
      dplyr::mutate(!!var_y := as.character(.data[[var_y]])) %>%
      dplyr::bind_rows(summary_bar) %>%
      dplyr::mutate(!!var_y := factor(relabel(.data[[var_y]]), levels = y_levels))

    clone_plot <- ggplot(clone_input) +
      geom_bar(position = fill_pos, aes(x = .data[[var_y]], fill = .data[[var_x]])) +
      scale_x_discrete(limits = rev) +
      coord_flip() +
      labs(caption = paste0(
        "Fisher's exact, each ", var_y, " vs all other cells; BH-adjusted across the ",
        sum(!is.na(stats_tbl$p_value)), " tested ", var_y,
        ".  *** q<0.001  ** q<0.01  * q<0.05  . q<0.1;  ns = not significant, ",
        "n<", signif_min_cells, " = too few cells to test")) +
      theme_minimal()

  } else {
    stats_tbl <- NULL
    clone_plot <- ggplot(seu_meta) +
      geom_bar(position = fill_pos, aes(x = .data[[var_y]], fill = .data[[var_x]])) +
      geom_bar(data = summarized_clusters, position = fill_pos, aes(x = .data[[var_y]], fill = .data[[var_x]])) +
      scale_x_discrete(limits = rev) +
      coord_flip() +
      theme_minimal()
  }

  if (!is.null(avg_line)) {
    clone_plot <- clone_plot + geom_hline(yintercept = avg_line)
  }

  plot_return <- switch(plot_type,
    clone = clone_plot,
    cluster = cluster_plot,
    both = (clone_plot / cluster_plot) +
      plot_layout(ncol = 1) +
      plot_annotation(title = seu_name)
  )

  # Hand the full table back on whatever object is actually returned, so the
  # stars in the PDF are never the only record of the test. It has to be attached
  # here rather than to `clone_plot`: under plot_type = "both" the return value is
  # a fresh patchwork built from clone_plot, and `/` does not carry attributes
  # across.
  if (!is.null(stats_tbl)) {
    attr(plot_return, "enrichment_stats") <- stats_tbl
  }

  return(plot_return)
}
