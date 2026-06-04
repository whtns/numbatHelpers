# Plot Functions (147)

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param clone_simplifications Parameter for clone simplifications
#' @param label Character string (default: "_clone_tree")
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
save_cc_space_plot_from_path <- function(seu_path, clone_simplifications, label = "_clone_tree", ...) {
  # Handle NA inputs (e.g., when clone_post is NULL for a sample)
  if (is.na(seu_path)) {
    return(NA_character_)
  }
  
  seu <- readRDS(seu_path)
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  plot_cc_space_plot(seu, tumor_id = tumor_id, sample_id = sample_id, ...)

  plot_path <- ggsave(glue("results/{sample_id}{label}.pdf"), width = 4, height = 4)
  return(plot_path)
}

#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param tumor_id Parameter for tumor id
#' @param nb_path File path
#' @param clone_simplifications Parameter for clone simplifications
#' @param sample_id Parameter for sample id
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_clone_tree <- function(clone_df, tumor_id, nb_path, clone_simplifications = NULL, sample_id = NULL, show_distance = FALSE, ...) {
  if (!"clone_opt" %in% colnames(clone_df)) {
    warning("clone_opt not found for ", tumor_id, "; skipping clone tree")
    return(NULL)
  }

  mynb <- readRDS(nb_path)

  mynb$mut_graph <-
    mynb$mut_graph |>
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(clone %in% unique(clone_df$clone_opt)) %>%
    as.igraph() %>%
    identity()

  mynb$clone_post <- dplyr::filter(mynb$clone_post, cell %in% clone_df$cell)

  ## clone tree ------------------------------

  if(!is.null(clone_simplifications)){
  	rb_scnas <- clone_simplifications[[tumor_id]]
  	mynb <- simplify_gt(mynb, rb_scnas)
  }

  # Renumber clones in BFS tree order so clone numbers follow phylogenetic order
  g <- mynb$mut_graph
  root_v <- which(igraph::degree(g, mode = "in") == 0)
  bfs_order <- igraph::bfs(g, root = root_v, mode = "out")$order
  old_clones <- igraph::V(g)$clone[bfs_order]
  remap <- setNames(seq_along(old_clones), as.character(old_clones))

  igraph::V(mynb$mut_graph)$clone <- as.integer(remap[as.character(igraph::V(mynb$mut_graph)$clone)])
  mynb$clone_post <- mynb$clone_post %>%
    dplyr::mutate(clone_opt = as.integer(remap[as.character(clone_opt)]))
  clone_df$clone_opt <- as.integer(remap[as.character(clone_df$clone_opt)])

  nclones <- max(as.integer(unique(clone_df$clone_opt)), na.rm = TRUE)
  mypal <- scales::hue_pal()(nclones) %>%
    set_names(1:nclones)

  plot_title <- ifelse(is.null(sample_id), tumor_id, sample_id)

  clone_plot <- mynb$plot_mut_history(
    pal = mypal,
    show_distance = show_distance,
    ...
  ) +
    labs(title = plot_title) +
    theme(plot.title = element_text(hjust = 0.5))

  # Add background to edge labels using geom_label at edge midpoints
  # Extract edge data from ggplot_build, calculate midpoints, and wrap label text
  edge_data <- ggplot2::ggplot_build(clone_plot)$data[[1]]
  edge_labels <- edge_data %>%
    group_by(group) %>%
    summarise(
      x = mean(x),
      y = mean(y),
      label = unique(label)
    ) %>%
    mutate(label = str_wrap(sub('.*-> *', '', label), width = 10)) %>%
    filter(label != "")

  clone_plot <- clone_plot +
    geom_label(
      data = edge_labels,
      aes(x = x, y = y, label = label),
      fill = "white",
      color = "black",
      label.size = 0.2,
      label.padding = unit(0.15, "lines"),
      size = 3,
      na.rm = TRUE
    )

  # Remove the original label mapping from edge layers to avoid double labels
  for (i in seq_along(clone_plot$layers)) {
    if (inherits(clone_plot$layers[[i]]$geom, "GeomEdgePath")) {
      mapping <- clone_plot$layers[[i]]$mapping
      if (is.list(mapping)) {
        clone_plot$layers[[i]]$mapping <- mapping[base::setdiff(names(mapping), "label")]
      }
    }
  }

  clone_plot$data$clone <- clone_plot$data$id

  return(clone_plot)
}

#' Perform differential expression analysis
#'
#' @param to_SCT_snn_res. Parameter for to SCT snn res.
#' @param to_clust Character string (default: "1_10")
#' @param sample_id Parameter for sample id
#' @param tumor_id Parameter for tumor id
#' @param seu Seurat object
#' @param mynb Numbat object
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
find_diffex_from_clustree <- function(to_SCT_snn_res. = 1, to_clust = "1_10", sample_id, tumor_id, seu, mynb, ...) {
  to_clust <- str_split(to_clust, pattern = "_") %>%
    unlist()

  to_SCT_snn_res. <- glue("SCT_snn_res.{to_SCT_snn_res.}")

  divergent_diffex <- find_diffex_bw_divergent_clusters(sample_id, tumor_id, seu, mynb, to_SCT_snn_res., to_clust, ...)

  #

  divergent_diffex <-
    divergent_diffex %>%
    # compact() %>%
    map(dplyr::bind_rows, .id = "location") %>%
    bind_rows(.id = "clone_comparison") %>%
    dplyr::mutate(sample_id := {{ sample_id }}) %>%
    # dplyr::arrange(cluster, p_val_adj) %>%
    identity() %>%
    group_by(to_clust, clone_comparison, location)

  group_names <-
    divergent_diffex %>%
    group_keys() %>%
    dplyr::mutate(plot_label = glue("clusters: {to_clust}; clones: {clone_comparison}; location: {location}")) %>%
    dplyr::select(plot_label) %>%
    tibble::deframe()

  volcano_plots <-
    divergent_diffex %>%
    group_split() %>%
    set_names(group_names) %>%
    map(tibble::column_to_rownames, "symbol") %>%
    map(dplyr::mutate, diffex_comparison = to_clust) %>%
    imap(make_volcano_plots, sample_id = sample_id) %>%
    identity()

  pdf_path <- glue("results/divergent_cluster_diffex_{sample_id}_{to_SCT_snn_res.}_{paste(to_clust, collapse = '_')}.pdf")
  pdf(pdf_path)
  print(volcano_plots)
  dev.off()

  # enrichment_table <-
  # 	diffex %>%
  # 	dplyr::distinct(symbol, .keep_all = TRUE) %>%
  # 	tibble::column_to_rownames("symbol") %>%
  # 	dplyr::select(-any_of(colnames(annotables::grch38))) %>%
  # 	enrichment_analysis() %>%
  # 	setReadable(org.Hs.eg.db::org.Hs.eg.db, keyType = "auto")
  #
  # enrichment_plot <- ggplotify::as.ggplot(
  # 	plot_enrichment(enrichment_table)
  # ) +
  # 	labs(title = glue("{sample_id}_{unique(diffex$to_clust)}_{unique(diffex$to_SCT_snn_res.)}"))


  # return(list("diffex" = divergent_diffex, "enrichment_table" = enrichment_table, "enrichment_plot" = enrichment_plot))

  return(list("diffex" = divergent_diffex, "plot" = pdf_path))
}

#' Perform differential expression analysis
#'
#' @param table_set Parameter for table set
#' @param debranched_seus Parameter for debranched seus
#' @param ... Additional arguments passed to other functions
#' @return Differential expression results
#' @export
find_all_diffex_from_clustree <- function(table_set, debranched_seus, ...) {
  tumor_id <- str_extract(names(table_set), "SR[RX][0-9]+")
  message(tumor_id)

  sample_id <- names(table_set)
  message(sample_id)

  table_set <- table_set[[1]] %>%
    dplyr::distinct(to_SCT_snn_res., to_clust, .keep_all = TRUE)

  debranched_seus <-
    debranched_seus %>%
    unlist() %>%
    set_names(str_remove(fs::path_file(.), "_filtered_seu.*"))

  seu <- readRDS(debranched_seus[[sample_id]])

  mynb <- readRDS(glue("output/numbat_sridhar/{tumor_id}_numbat.rds"))

  
#' Perform check table set operation
#'
#' @param to_clust Parameter for to clust
#' @param to_SCT_snn_res. Parameter for to SCT snn res.
#' @param seu Seurat object
#' @return Function result
#' @export
check_table_set <- function(to_clust, to_SCT_snn_res., seu) {
    #
    idents <-
      to_clust %>%
      str_split(pattern = "_") %>%
      unlist()

    to_SCT_snn_res. <- glue("SCT_snn_res.{to_SCT_snn_res.}")

    all(idents %in% seu@meta.data[[to_SCT_snn_res.]])
  }

  table_set <- table_set %>%
    dplyr::rowwise() %>%
    dplyr::mutate(good_set = check_table_set(to_clust, to_SCT_snn_res., seu)) %>%
    identity()

  message("running comparison")
  test0 <- purrr::map2(table_set$to_SCT_snn_res., table_set$to_clust, find_diffex_from_clustree, sample_id, tumor_id, seu, mynb, ...)

  return(test0)
}