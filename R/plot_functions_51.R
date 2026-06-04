# Plot Functions (150)

#' Merge or join datasets
#'
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
merge_orders <- function(...) {
  overall_orders <- c(...)

  overall_orders <- overall_orders[!duplicated(names(overall_orders))]

  overall_orders <- overall_orders[order(names(overall_orders))]
}

#' Perform differential expression analysis
#'
#' @param diffex Parameter for diffex
#' @param p_val_limit Parameter for p val limit
#' @param log2fc_limit Parameter for log2fc limit
#' @return Differential expression results
#' @export
clean_diffex <- function(diffex, p_val_limit = 0.05, log2fc_limit = 1) {
  diffex %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::select(description, everything()) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::filter(abs(avg_log2FC) > log2fc_limit) %>%
    dplyr::filter(p_val_adj < p_val_limit) %>%
    dplyr::arrange(desc(avg_log2FC), p_val_adj) %>%
    identity()
}

#' Create a plot visualization
#'
#' @param ident.1 Cell identities or groups
#' @param ident.2 Cell identities or groups
#' @param scna_of_interest Parameter for scna of interest
#' @param seu Seurat object
#' @param sample_id Parameter for sample id
#' @return List object
#' @export
table_and_plot_enrichment <- function(ident.1, ident.2, scna_of_interest, seu, sample_id) {
  #
  my_diffex <- Seurat::FindMarkers(seu, group.by = "clusters", ident.1 = ident.1, ident.2 = ident.2)

  my_enrichment <-
    my_diffex %>%
    # clean_diffex() %>%
    enrichment_analysis(gene_set = "hallmark") %>%
    DOSE::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

  core_enriched_genes <- my_enrichment %>%
    as_tibble() %>%
    dplyr::mutate(core_enrichment = str_split(core_enrichment, "\\/")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(core_enrichment = list(core_enrichment[core_enrichment %in% rownames(my_diffex)])) %>%
    dplyr::mutate(core_enrichment = list(paste(core_enrichment, collapse = "/"))) %>%
    dplyr::pull(core_enrichment) %>%
    identity()

  simplified_enrichment <- my_enrichment

  simplified_enrichment@result$core_enrichment <- core_enriched_genes

  my_enrichment_plot <-
    my_enrichment %>%
    plot_enrichment() +
    labs(title = glue("{scna_of_interest}: {ident.1} v. {ident.2}"))

  netplot <- simplified_enrichment %>%
    enrichplot::cnetplot(
      node_label = "all",
      showCategory = 10,
      # cex.params = list(gene_node = 2),
      cex_category = 0.5,
      cex_gene = 0.5,
      cex_label_category = 0.6,
      cex_label_gene = 0.3
    ) +
    labs(title = glue("{scna_of_interest}: {ident.1} v. {ident.2}"))

  return(
    list(
      "diffex" = clean_diffex(my_diffex),
      "enrichment_table" = as_tibble(my_enrichment),
      "dotplot" = my_enrichment_plot,
      "netplot" = netplot
    )
  )
}

#' Perform clustering analysis
#'
#' @param cluster_comparisons_file File path
#' @return Extracted data elements
#' @export
pull_cluster_comparisons <- function(cluster_comparisons_file) {
  cluster_comparisons <- read_csv(cluster_comparisons_file) %>%
    split(.$sample_id)
}

#' Perform clustering analysis
#'
#' @param cluster_comparison Cluster information
#' @param debranched_seus Parameter for debranched seus
#' @return Function result
#' @export
make_cluster_comparisons_by_phase_for_disctinct_clones <- function(cluster_comparison, debranched_seus) {
  sample_id <- unique(cluster_comparison[[1]][["sample_id"]])

  tumor_id <- str_extract(sample_id, "SR[RX][0-9]+")

  debranched_seus <- debranched_seus %>%
    set_names(str_extract(., "SR[RX][0-9]+.*(?=_filtered_seu.rds)"))

  seu_path <- debranched_seus[[sample_id]]

  test0 <-
    cluster_comparison[[1]] %>%
    split(.$scna_of_interest) %>%
    map(dplyr::select, all_of(c("ident.1", "ident.2", "scna_of_interest"))) %>%
    as.list() %>%
    identity()

  seu <- readRDS(seu_path)

  safe_table_and_plot_enrichment <- safely(table_and_plot_enrichment)

  plot_outcome <- map(test0, ~ pmap(.x, safe_table_and_plot_enrichment, seu, sample_id))

  outfile <- glue("results/cluster_comparisons_by_phase_for_disctinct_clones_{sample_id}.pdf")

  slugs <-
    test0 %>%
    map(dplyr::select, ident.1, ident.2) %>%
    map(tidyr::unite, "slug", everything(), sep = "-v-") %>%
    map(dplyr::pull, slug) %>%
    identity()

  plot_outcome <- map2(plot_outcome, slugs, ~ set_names(.x, .y))

  results <-
    plot_outcome %>%
    list_flatten() %>%
    map("result") %>%
    purrr::compact() %>%
    map(~ .x[c("diffex", "enrichment_table")]) %>%
    imap(~ write_xlsx(.x, glue("results/{sample_id}_{.y}.xlsx"))) %>%
    identity()

  pdf(outfile, height = 8, width = 8)
  print(plot_outcome)
  dev.off()

  return(outfile)
}

