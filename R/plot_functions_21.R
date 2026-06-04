# Plot Functions (114)

#' Add annotations to data
#'
#' @param seu Seurat object
#' @return Modified Seurat object
#' @export
annotate_seu_with_rb_subtype_gene_expression <- function(seu) {
  hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

  subtype_genes <- list(
    gp1 = c("EGF", "TPBG", "GUCA1C", "GUCA1B", "GUCA1A", "GNAT2", "GNGT2", "ARR3", "PDE6C", "PDE6H", "OPN1SW"),
    gp2 = c("TFF1", "CD24", "EBF3", "GAP43", "STMN2", "POU4F2", "SOX11", "EBF1", "DCX", "ROBO1", "PCDHB10")
  )

  seu <- AddModuleScore(
    object = seu,
    features = subtype_genes,
    ctrl = 5,
    name = "exprs_gp"
  )

  subtype_hallmarks <- c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME", "HALLMARK_BILE_ACID_METABOLISM",
    "HALLMARK_PROTEIN_SECRETION", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_MYC_TARGETS_V2"
  )

  
#' Extract or pull specific data elements
#'
#' @param hallmark_gene_set Gene names or identifiers
#' @return Extracted data elements
#' @export
pull_hallmark_genes <- function(hallmark_gene_set) {
    hallmark_gene_sets$gene_symbol[hallmark_gene_sets$gs_name == hallmark_gene_set]
  }

  hallmark_genes <- lapply(subtype_hallmarks, function(hm) intersect(pull_hallmark_genes(hm), rownames(seu)))
  names(hallmark_genes) <- subtype_hallmarks

  seu <- AddModuleScore(
    object = seu,
    features = hallmark_genes,
    ctrl = 5,
    name = "hallmark"
  )

  meta_names <- names(seu@meta.data)
  hallmark_idx <- which(meta_names %in% paste0("hallmark", seq_along(hallmark_genes)))
  names(seu@meta.data)[hallmark_idx] <- names(hallmark_genes)

  return(seu)
}
#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param seu_name Parameter for seu name
#' @param var_x Character string (default: "scna")
#' @param var_y Character string (default: "SCT_snn_res.0.6")
#' @param avg_line Parameter for avg line
#' @return ggplot2 plot object
#' @export
plot_distribution_of_clones_pearls <- function(seu, seu_name, var_x = "scna", var_y = "SCT_snn_res.0.6", avg_line = NULL) {
  seu_meta <- seu@meta.data
  phase_levels <- seu_meta$phase_level

  df <- seu_meta %>%
    group_by(.data[[var_x]], .data[[var_y]]) %>%
    dplyr::count() %>%
    group_by(.data[[var_x]]) %>%
    mutate(
      per = proportions(n) * 100,
      phase = factor(str_remove(.data[[var_y]], "_[0-9]*$"), levels = levels(phase_levels)),
      label = paste0(.data[[var_y]], " Fraction: \n", round(per, 2), "%"),
      label_lag = dplyr::lag(per, default = 0),
      label_position = cumsum(per) - 0.5 * per,
      scna = factor(scna, levels = c("all", levels(seu_meta$scna)))
    )

  y_setting <- df %>%
    group_by(scna, phase) %>%
    summarize(n = n(), .groups = "drop")

  pearls_plots <- make_pearls_plot(df, y_setting = max(y_setting$n) + 2, var_y = var_y)
  return(pearls_plots)
}

#' Create a plot visualization
#'
#' @param df Input data frame or dataset
#' @param y_setting Parameter for y setting
#' @param var_y Character string (default: "clusters")
#' @return ggplot2 plot object
#' @export
make_pearls_plot <- function(df, y_setting = 4, var_y = "clusters") {
  
#' Perform shifter operation
#'
#' @param x Parameter for x
#' @param n Parameter for n
#' @return Function result
#' @export
shifter <- function(x, n = 1) {
    if (n == 0) x else c(tail(x, -n), head(x, n))
  }

  only_phases_df <-
    df %>%
    mutate(cluster_number = str_extract(.data[[var_y]], "[0-9]*$")) |>
    dplyr::filter(!phase %in% c("hsp", "hypoxia", "other", "s_star")) |>
    group_by(scna, phase) |>
    mutate(cluster_sequence = dplyr::row_number()) |>
    mutate(cluster_sequence = scale(cluster_sequence, scale = FALSE) + y_setting) |>
    dplyr::ungroup() |>
    tidyr::nest(data = -phase) |>
    dplyr::mutate(to_x = dplyr::lead(phase)) |>
    dplyr::mutate(to_x = dplyr::coalesce(to_x, phase)) |>
    # dplyr::mutate(to_x = shifter(phase)) |>
    tidyr::unnest(cols = c(data)) |>
    arrange(scna, .data[[var_y]]) |>
    identity()

  non_phases_df <-
    df |>
    dplyr::filter(phase %in% c("hsp", "hypoxia", "other", "s_star")) |>
    mutate(cluster_number = str_extract(.data[[var_y]], "[0-9]*$")) |>
    group_by(scna, phase) |>
    # mutate(cluster_sequence = runif(dplyr::n())) |>
    mutate(cluster_sequence = dplyr::row_number()) |>
    mutate(cluster_sequence = scale(cluster_sequence, scale = FALSE) + 1) |>
    group_by(scna) |>
    mutate(phase = sample(unique(only_phases_df$phase), dplyr::n(), replace = TRUE)) |>
    mutate(cluster_number = .data[[var_y]]) |>
    identity()

  n_groups <- n_distinct(only_phases_df$phase)

  scale_breaks <- (0:floor(max(only_phases_df$per) / 10) * 10)

  ggplot(only_phases_df, aes(x = phase, y = cluster_sequence, size = per)) +
    geom_rect(aes(xmin = 0, xmax = n_groups + 0.5, ymin = 0, ymax = 1.9), fill = "white") +
    geom_segment(data = only_phases_df, aes(x = phase, y = y_setting, xend = to_x, yend = y_setting), size = 1, arrow = arrow(angle = 15, ends = "last", length = unit(0.03, "npc"), type = "closed"), alpha = 0.5) +
    # geom_point(aes(color = cluster_number)) +
    geom_jitter(aes(color = cluster_number), width = 0.05) +
    # geom_text_repel(aes(x = phase, y=cluster_sequence, label = cluster_number), size = 4, color = "black") +
    geom_point(data = non_phases_df, aes(x = phase, y = cluster_sequence, color = cluster_number)) +
    geom_text_repel(data = non_phases_df, aes(x = phase, y = cluster_sequence, label = cluster_number), size = 4, color = "black") +
    coord_polar() +
    scale_size(range = c(0, max(scale_breaks) / 6)) +
    facet_wrap(~scna, nrow = 1) +
    ylim(0, y_setting + 2) + # adjust as you like
    theme_minimal() +
    guides(
      color = "none",
      size = guide_legend(override.aes = list(alpha = 0.2), theme = theme(legend.key = element_rect(colour = NA, fill = NA)))
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      text = element_text(size = 14),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      legend.spacing.y = unit(0.03, "npc"),
      legend.position = "top",
      panel.spacing.x = unit(1, "cm")
    ) +
    labs(size = "%") +
    NULL

  # ggsave("~/tmp/test.pdf", height = 8, width = 8) |>
  # 	browseURL()
}