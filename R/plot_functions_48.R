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
  # Accept a Seurat object as well as a cell/clone_opt data frame. Several callers
  # (e.g. plot_seu_marker_heatmap) pass `seu` directly; colnames(seu) are cell
  # barcodes, so the clone_opt check below always failed and returned NULL, which
  # then crashed the collage's wrap_plots() with "Only know how to add
  # <ggplot>/<grob> objects". Extract the metadata table (cell + clone_opt) here,
  # normalising the "." -> "-" barcode format numbat's clone_post uses.
  if (inherits(clone_df, "Seurat")) {
    clone_df <- clone_df@meta.data %>% tibble::rownames_to_column("cell")
    clone_df$cell <- stringr::str_replace(clone_df$cell, "\\.", "-")
  }
  if (!"clone_opt" %in% colnames(clone_df)) {
    warning("clone_opt not found for ", tumor_id, "; skipping clone tree")
    return(NULL)
  }

  # nb_path is normally a path to a *_numbat.rds, but may also be an in-memory
  # numbat object. The per-round iteration summaries hold a Numbat built at one
  # specific consensus round, and re-reading the RDS from disk would silently
  # substitute the SELECTED round's tree for round k's -- the exact class of bug
  # that made SRX11133593's clone keys wrong.
  #
  # R6 objects are environments, so the mut_graph/clone_post assignments below
  # would mutate the CALLER's object in place. Deep-clone before touching it.
  mynb <- if (is.character(nb_path)) {
    readRDS(nb_path)
  } else if (inherits(nb_path, "R6")) {
    nb_path$clone(deep = TRUE)
  } else {
    nb_path
  }

  # Node count of the FULL numbat tree, captured before the filter below, so the
  # steps that follow can tell a whole-tree plot from a subset one.
  n_nodes_full <- igraph::vcount(mynb$mut_graph)

  mynb$mut_graph <-
    mynb$mut_graph |>
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(clone %in% unique(clone_df$clone_opt)) %>%
    igraph::as.igraph() %>%
    identity()

  mynb$clone_post <- dplyr::filter(mynb$clone_post, cell %in% clone_df$cell)

  # Did the filter drop any clone? The two-clone SCNA collages draw this panel
  # beside a clone UMAP, a clone heatmap annotation and a stacked bar captioned
  # "6p+ (clone 7)" / "preceding (clone 4)", all of which label cells by
  # `clone_opt`, so a subset tree has to keep speaking numbat clone ids.
  subset_tree <- igraph::vcount(mynb$mut_graph) != n_nodes_full

  ## clone tree ------------------------------

  if (!is.null(clone_simplifications)) {
    # Support both whole-dict (keyed by tumor_id) and already-extracted per-sample form
    rb_scnas <- if (tumor_id %in% names(clone_simplifications)) {
      clone_simplifications[[tumor_id]]
    } else {
      clone_simplifications
    }
    if (!is.null(rb_scnas) && length(rb_scnas) > 0) {
      mynb <- simplify_gt(mynb, rb_scnas)
    }
  }

  # No renumbering. numbat's clone ids ARE phylogenetic order already: the
  # mut_graph `clone` attribute is the DFS visit order from the root, and
  # clone_post$clone_opt uses the same ids, so nodes are labelled with the ids
  # every other panel (clone UMAP, heatmap annotation, clone_labels) uses. This
  # used to renumber clone_post in BFS order, which plot_mut_history()'s
  # label_genotype() then overwrote with DFS order on the graph side, so nodes
  # were sized by the wrong clone wherever BFS != DFS (12 of 33 t=1e-5 graphs).

  # One hue per DISPLAYED clone. plot_mut_history() colours by numbat's
  # label_genotype() numbering, which is always 1..n over the nodes it is given,
  # so the palette has to be named 1..n and sized to the nodes on show. For a
  # two-clone subset this yields the same two hues, in the same order, as the
  # two-level clone UMAP beside it instead of the first two of a 7-colour ramp.
  .clones_shown <- unique(stats::na.omit(as.integer(clone_df$clone_opt)))
  n_display <- if (subset_tree) length(.clones_shown) else max(.clones_shown)
  mypal <- scales::hue_pal()(n_display) %>%
    set_names(seq_len(n_display))

  plot_title <- ifelse(is.null(sample_id), tumor_id, sample_id)

  # plot_mut_history() str_trunc()s edge labels to 20 characters ("17p-,...18q+"),
  # so the edges carry short keys through it and the full labels are restored
  # below. Edge lengths are set from the real labels first: with show_distance,
  # numbat would otherwise count mutations in the keys.
  edge_full <- igraph::E(mynb$mut_graph)$to_label
  edge_keys <- sprintf("e%d", seq_along(edge_full))
  mynb$mut_graph <- mynb$mut_graph %>%
    igraph::set_edge_attr("length", value = lengths(str_split(edge_full, ","))) %>%
    igraph::set_edge_attr("to_label", value = edge_keys)

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
    mutate(label = unname(setNames(edge_full, edge_keys)[label])) %>%
    filter(!is.na(label)) %>%
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

  # numbat's plot_mut_history() runs the graph through label_genotype(), which
  # renumbers the nodes 1..n in DFS order; the plot's `clone` column is what the
  # node text is mapped to. For a whole tree that IS the numbat clone id (DFS
  # order is how numbat assigned the ids), so it is left alone. Never use `id`
  # here: it is the vertex index, which differs from the clone id in 25 of 33
  # t=1e-5 graphs.
  #
  # For a SUBSET tree the 1..n renumbering is wrong: clones 4 and 7 of
  # SRX14116946 would be drawn as "1 -> 2" beside panels labelled clone 4 and
  # clone 7. Relabel a subset tree from the graph's own clone ids.
  if (subset_tree) {
    clone_ids <- setNames(igraph::V(mynb$mut_graph)$clone,
                          as.character(igraph::V(mynb$mut_graph)$id))
    clone_plot$data$clone <- unname(clone_ids[as.character(clone_plot$data$id)])
  }

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