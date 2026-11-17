#' Cohort summary of what each filter removes, by numbat genotype (issue #54)
#'
#' Reads the per-cell tibbles from `extract_filter_metadata()` (the
#' `filter_inspection_metadata` target) and plots, per sample, how many cells each filter
#' removes from tumor, diploid and un-genotyped cells, with a page on the per-cell cell-type
#' exclusion (M2 call + marker-panel or cluster-majority support).
#'
#' @param filter_meta List (or data.frame) of per-sample tibbles from `extract_filter_metadata()`.
#' @param out_pdf Output PDF.
#' @param out_csv Output CSV: sample_id, role, fate, n, pct_of_role.
#' @return Character vector c(out_pdf, out_csv).
#' @export
plot_filtering_impact <- function(filter_meta,
                                  out_pdf = "results/filtering_impact.pdf",
                                  out_csv = "results/filtering_impact.csv") {
  if (!is.data.frame(filter_meta)) {
    filter_meta <- Filter(is.data.frame, filter_meta)
    filter_meta <- dplyr::bind_rows(filter_meta)
  }
  need <- c("sample_id", "filter_reason", "GT_opt")
  miss <- setdiff(need, names(filter_meta))
  if (length(miss)) stop("plot_filtering_impact: filter metadata lacks ", paste(miss, collapse = ", "))
  if (!"celltype_m2_group" %in% names(filter_meta)) filter_meta$celltype_m2_group <- NA_character_
  if (!"celltype_exclude" %in% names(filter_meta)) filter_meta$celltype_exclude <- FALSE

  fate_levels <- c("kept", "clone_opt_na", "qc_fail", "manual_exclude", "celltype_m2")
  role_levels <- c("tumor", "diploid", "no_numbat")
  d <- filter_meta |>
    dplyr::mutate(
      role = factor(dplyr::case_when(is.na(GT_opt) ~ "no_numbat", GT_opt == "" ~ "diploid",
                                     TRUE ~ "tumor"), levels = role_levels),
      fate = factor(tidyr::replace_na(as.character(filter_reason), "kept"), levels = fate_levels),
      excl_group = celltype_m2_group %in% celltype_m2_excluded_groups(),
      celltype_exclude = celltype_exclude %in% TRUE)

  counts <- d |>
    dplyr::count(sample_id, role, fate, name = "n", .drop = FALSE) |>
    dplyr::group_by(sample_id, role) |>
    dplyr::mutate(pct_of_role = if (sum(n) > 0) round(100 * n / sum(n), 2) else NA_real_) |>
    dplyr::ungroup()
  fs::dir_create(dirname(out_csv))
  readr::write_csv(counts, out_csv)

  fate_cols <- c(kept = "grey80", clone_opt_na = "#8c564b", qc_fail = "#ff7f0e",
                 celltype_m2 = "#d62728",
                 manual_exclude = "#17becf")
  # order samples by tumor-cell count
  ord <- d |> dplyr::filter(role == "tumor") |> dplyr::count(sample_id) |>
    dplyr::arrange(dplyr::desc(n)) |> dplyr::pull(sample_id)
  ord <- c(ord, setdiff(unique(d$sample_id), ord))
  counts$sample_id <- factor(counts$sample_id, levels = rev(ord))

  base_theme <- ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", panel.grid.major.y = ggplot2::element_blank())

  p_n <- ggplot2::ggplot(dplyr::filter(counts, n > 0),
                         ggplot2::aes(x = n, y = sample_id, fill = fate)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~role, nrow = 1, scales = "free_x") +
    ggplot2::scale_fill_manual(values = fate_cols, drop = FALSE) +
    ggplot2::labs(x = "cells", y = NULL, fill = "fate",
                  title = "Cells by filter outcome, per sample and numbat genotype",
                  subtitle = "Unfiltered object; each cell takes its first applicable filter reason") +
    base_theme

  p_pct <- ggplot2::ggplot(dplyr::filter(counts, n > 0),
                           ggplot2::aes(x = pct_of_role, y = sample_id, fill = fate)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~role, nrow = 1) +
    ggplot2::scale_fill_manual(values = fate_cols, drop = FALSE) +
    ggplot2::labs(x = "% of the sample's cells in that genotype", y = NULL, fill = "fate",
                  title = "Filter outcome as a share of each genotype") +
    base_theme

  # cell-type pages
  ct_tot <- d |> dplyr::group_by(sample_id, role) |> dplyr::summarise(n_role = dplyr::n(), .groups = "drop")
  ct <- d |>
    dplyr::filter(excl_group) |>
    dplyr::mutate(status = ifelse(fate == "celltype_m2", "removed (call + support)",
                           ifelse(celltype_exclude, "excluded call + support, removed earlier",
                                  "no support (kept by this rule)"))) |>
    dplyr::count(sample_id, role, celltype_m2_group, status, name = "n") |>
    dplyr::left_join(ct_tot, by = c("sample_id", "role")) |>
    dplyr::mutate(pct = 100 * n / n_role,
                  sample_id = factor(sample_id, levels = rev(ord)))
  tot <- d |> dplyr::filter(fate == "celltype_m2") |> dplyr::count(role)
  tot_txt <- paste0("Removed as celltype_m2: ",
                    paste(tot$role, tot$n, sep = " ", collapse = ", "),
                    " (tumor ", round(100 * sum(tot$n[tot$role == "tumor"]) /
                                        max(1, sum(d$role == "tumor")), 2), "% of all tumor cells)")

  group_cols <- c(rod = "#1f77b4", bipolar = "#17becf", rbc = "#d62728", microglia = "#2ca02c",
                  other_immune = "#98df8a", other_nonretinal = "#9467bd", rpe = "#8c564b")
  mk_ct <- function(dd, title, subtitle) {
    if (!nrow(dd)) return(ggplot2::ggplot() + ggplot2::labs(title = title, subtitle = "none") + base_theme)
    ggplot2::ggplot(dd, ggplot2::aes(x = n, y = sample_id, fill = celltype_m2_group)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::facet_wrap(~role, nrow = 1, scales = "free_x", drop = FALSE) +
      ggplot2::scale_fill_manual(values = group_cols) +
      ggplot2::geom_text(data = dd |> dplyr::group_by(sample_id, role) |>
                           dplyr::summarise(n = sum(n), pct = sum(pct), .groups = "drop"),
                         ggplot2::aes(x = n, y = sample_id, label = sprintf("%.1f%%", pct)),
                         inherit.aes = FALSE, hjust = -0.1, size = 2.3) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
      ggplot2::labs(x = "cells", y = NULL, fill = "M2 group", title = title, subtitle = subtitle) +
      base_theme
  }
  p_rm <- mk_ct(dplyr::filter(ct, status == "removed (call + support)"),
                "Cells removed by the cell-type rule (M2 call + marker panel or >=60% of its cluster)",
                paste0(tot_txt, ". Label: % of that sample's cells in the genotype."))
  p_spared <- mk_ct(dplyr::filter(ct, status == "no support (kept by this rule)"),
                    "M2 calls an excluded type but no marker or cluster-majority support: not removed by the cell-type rule",
                    "These cells may still be removed by another filter. Label: % of that sample's cells in the genotype.")

  # fine M2 types (all 30 plae types) ----------------------------------------
  ct_pages <- list(); bar_pages <- list()
  ct_csv <- sub("\\.csv$", "_celltypes.csv", out_csv)
  if ("celltype_m2" %in% names(d)) {
    group_order <- c("cone", "pr_precursor", "rpc", "neurogenic", "other_retinal", "rod", "bipolar",
                     "rbc", "microglia", "other_immune", "other_nonretinal", "rpe")
    dd <- d |> dplyr::filter(role %in% c("tumor", "diploid"), !is.na(celltype_m2)) |>
      dplyr::mutate(role = droplevels(role))
    fine <- dd |>
      dplyr::group_by(sample_id, role, celltype_m2, celltype_m2_group) |>
      dplyr::summarise(n_all = dplyr::n(), n_kept = sum(fate %in% "kept"),
                       n_removed_celltype_m2 = sum(fate %in% "celltype_m2"),
                       n_removed_other = sum(!fate %in% c("kept", "celltype_m2")), .groups = "drop")
    readr::write_csv(fine, ct_csv)
    type_tab <- dplyr::distinct(fine, celltype_m2, celltype_m2_group) |>
      dplyr::mutate(g = factor(celltype_m2_group, levels = group_order)) |>
      dplyr::arrange(g, celltype_m2) |>
      dplyr::mutate(lab = ifelse(celltype_m2_group %in% celltype_m2_excluded_groups(),
                                 paste0(celltype_m2, " *"), celltype_m2))
    mk_heat <- function(col, title) {
      h <- fine |>
        dplyr::mutate(n = .data[[col]]) |>
        dplyr::group_by(sample_id, role) |>
        dplyr::mutate(pct = 100 * n / sum(n)) |>
        dplyr::ungroup() |>
        dplyr::filter(n > 0) |>
        dplyr::mutate(type = factor(type_tab$lab[match(celltype_m2, type_tab$celltype_m2)],
                                    levels = type_tab$lab),
                      group = factor(celltype_m2_group, levels = group_order),
                      sample_id = factor(sample_id, levels = rev(ord)))
      ggplot2::ggplot(h, ggplot2::aes(x = type, y = sample_id, fill = pct)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = n), size = 1.9) +
        ggplot2::facet_grid(role ~ group, scales = "free", space = "free") +
        ggplot2::scale_fill_gradient(low = "#e8f1fb", high = "#d62728", trans = "sqrt",
                                     limits = c(0, 100), name = "% of sample's\ncells in genotype") +
        ggplot2::labs(x = NULL, y = NULL, title = title,
                      subtitle = "M2 call (plae, 30 types) per cell; text = cells. * = type the cell-type rule can exclude (needs marker or cluster-majority support).") +
        ggplot2::theme_minimal(base_size = 8) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                       strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, size = 6),
                       strip.text.y = ggplot2::element_text(face = "bold"),
                       panel.grid = ggplot2::element_blank(), panel.spacing = grid::unit(1, "pt"))
    }
    pal <- celltype_m2_palette()
    bar_lab <- stats::setNames(type_tab$lab, type_tab$celltype_m2)
    bar_cols <- stats::setNames(pal[type_tab$celltype_m2], type_tab$lab)
    bar_cols[is.na(bar_cols)] <- "grey50"
    mk_bar <- function(col, title) {
      b <- fine |>
        dplyr::mutate(n = .data[[col]]) |>
        dplyr::group_by(sample_id, role) |>
        dplyr::mutate(pct = if (sum(n) > 0) 100 * n / sum(n) else 0) |>
        dplyr::ungroup() |>
        dplyr::filter(n > 0) |>
        dplyr::mutate(type = factor(bar_lab[celltype_m2], levels = rev(type_tab$lab)),
                      sample_id = factor(sample_id, levels = rev(ord)))
      ggplot2::ggplot(b, ggplot2::aes(x = pct, y = sample_id, fill = type)) +
        ggplot2::geom_col(width = 0.85) +
        ggplot2::facet_wrap(~role, nrow = 1) +
        ggplot2::scale_fill_manual(values = bar_cols, breaks = type_tab$lab, drop = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
        ggplot2::labs(x = "% of the sample's cells in that genotype", y = NULL, fill = "M2 type",
                      title = title,
                      subtitle = "All 30 plae M2 types, in lineage order. * = type the cell-type rule can exclude.") +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 6),
                       legend.key.size = grid::unit(8, "pt"),
                       panel.grid.major.y = ggplot2::element_blank())
    }
    bar_pages <- list(mk_bar("n_all", "M2 cell type of every genotyped cell (before filtering)"),
                      mk_bar("n_kept", "M2 cell type of kept cells (after filtering)"))
    ct_pages <- list(mk_heat("n_all", "M2 cell type of every genotyped cell (before filtering)"),
                     mk_heat("n_kept", "M2 cell type of kept cells (what remains in the filtered objects)"))
  }

  fs::dir_create(dirname(out_pdf))
  h <- max(5, 0.18 * length(ord) + 2)
  grDevices::pdf(out_pdf, width = 13, height = h)
  on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
  print(p_n); print(p_pct)
  for (p in bar_pages) print(p)
  print(p_rm); print(p_spared)
  if (length(ct_pages)) {
    grDevices::dev.off()
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp, width = 16, height = max(8, 0.36 * length(ord) + 3))
    for (p in ct_pages) print(p)
    grDevices::dev.off()
    main <- tempfile(fileext = ".pdf"); file.copy(out_pdf, main)
    qpdf::pdf_combine(c(main, tmp), out_pdf)
    return(c(out_pdf, out_csv, ct_csv))
  }
  grDevices::dev.off()
  c(out_pdf, out_csv)
}

#' Per-sample filtering steps: cells remaining after each filter, by numbat genotype (#54)
#'
#' Replaces the unlabelled `plot_filtering_timeline()` bars. Steps follow the precedence in
#' `annotate_filter_reason()` (each cell is counted under the first filter that removes it):
#' no numbat genotype, QC, manual list, cell type (M2 call + marker or cluster-majority support). One panel per genotype (tumor / diploid / no genotype), each
#' with its own y axis; the label above each bar gives the cells that step removed.
#'
#' @param filter_meta List or data.frame of per-cell metadata with `sample_id`,
#'   `filter_reason` and `GT_opt` (e.g. the `filter_inspection_metadata` target).
#' @param out_dir Directory for per-sample PDFs (`<sample>_filtering_steps.pdf`).
#' @param out_pdf Combined PDF (one page per sample).
#' @return Character vector: the combined PDF followed by the per-sample PDFs.
#' @export
plot_filtering_steps <- function(filter_meta,
                                 out_dir = "results/filtering_steps",
                                 out_pdf = "results/filtering_timelines.pdf") {
  if (!is.data.frame(filter_meta)) filter_meta <- dplyr::bind_rows(Filter(is.data.frame, filter_meta))
  miss <- setdiff(c("sample_id", "filter_reason", "GT_opt"), names(filter_meta))
  if (length(miss)) stop("plot_filtering_steps: filter metadata lacks ", paste(miss, collapse = ", "))

  steps <- c(clone_opt_na = "no numbat\ngenotype", qc_fail = "QC", manual_exclude = "manual\nlist",
             celltype_m2 = "cell type\n(M2+support)")
  step_lab <- paste0("- ", steps)
  x_levels <- c("all", step_lab, "kept")
  role_cols <- c(tumor = "#1f77b4", diploid = "#2ca02c", `no genotype` = "grey60")

  d <- filter_meta |>
    dplyr::transmute(
      sample_id,
      role = factor(dplyr::case_when(is.na(GT_opt) ~ "no genotype", GT_opt == "" ~ "diploid",
                                     TRUE ~ "tumor"), levels = names(role_cols)),
      reason = as.character(filter_reason))

  one_sample <- function(ds, sid) {
    bars <- dplyr::bind_rows(lapply(0:length(steps), function(k) {
      gone <- names(steps)[seq_len(k)]
      ds |> dplyr::filter(!(reason %in% gone)) |> dplyr::count(role, name = "n", .drop = FALSE) |>
        dplyr::mutate(step = if (k == 0) "all" else step_lab[k])
    }))
    kept <- ds |> dplyr::filter(is.na(reason)) |> dplyr::count(role, name = "n", .drop = FALSE) |>
      dplyr::mutate(step = "kept")
    rm_n <- ds |> dplyr::filter(!is.na(reason)) |> dplyr::count(reason, role, name = "removed", .drop = FALSE) |>
      dplyr::mutate(step = step_lab[match(reason, names(steps))]) |> dplyr::select(-reason)
    n_all <- bars |> dplyr::filter(step == "all") |> dplyr::select(role, n_all = n)
    bars <- dplyr::bind_rows(bars, kept) |>
      dplyr::left_join(rm_n, by = c("role", "step")) |>
      dplyr::left_join(n_all, by = "role") |>
      dplyr::mutate(
        removed = tidyr::replace_na(removed, 0L),
        label = dplyr::case_when(
          step == "all" ~ format(n, big.mark = ",", trim = TRUE),
          step == "kept" ~ paste0(format(n, big.mark = ",", trim = TRUE), "\n(",
                                  ifelse(n_all > 0, round(100 * n / n_all, 1), 0), "%)"),
          removed == 0 ~ "0",
          TRUE ~ paste0("-", format(removed, big.mark = ",", trim = TRUE))),
        step = factor(step, levels = x_levels))
    ggplot2::ggplot(bars, ggplot2::aes(x = step, y = n, fill = role)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = label,
                                      colour = step != "all" & step != "kept" & removed > 0),
                         vjust = -0.25, size = 2.6, lineheight = 0.9, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "#c0392b")) +
      ggplot2::facet_wrap(~role, ncol = 1, scales = "free_y", drop = FALSE) +
      ggplot2::scale_fill_manual(values = role_cols, drop = FALSE, guide = "none") +
      ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.35))) +
      ggplot2::labs(x = NULL, y = "cells remaining",
                    title = paste0(sid, ": cells remaining after each filter, by numbat genotype"),
                    subtitle = "Each cell is counted under the first filter that removes it (left to right). Red: cells removed at that step.") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                     plot.subtitle = ggplot2::element_text(size = 8))
  }

  fs::dir_create(out_dir); fs::dir_create(dirname(out_pdf))
  sids <- sort(unique(d$sample_id))
  plots <- lapply(sids, function(s) one_sample(dplyr::filter(d, sample_id == s), s))
  per <- file.path(out_dir, paste0(sids, "_filtering_steps.pdf"))
  for (i in seq_along(sids)) ggplot2::ggsave(per[i], plots[[i]], width = 8, height = 8)
  grDevices::pdf(out_pdf, width = 8, height = 8)
  for (p in plots) print(p)
  grDevices::dev.off()
  c(out_pdf, per)
}
