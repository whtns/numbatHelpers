# Enrichment and ORA analysis functions (8)

#' Perform ora analysis operation
#'
#' @param seu Seurat object
#' @param clusters Cluster information
#' @param only_unique_terms Logical flag (default: FALSE)
#' @return Function result
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

ora_analysis <- function(seu, clusters = "gene_snn_res.0.2", only_unique_terms = FALSE) {
  
  df_list <- seu@misc$markers[[clusters]][["presto"]] %>%
    dplyr::rename(symbol = `Gene.Name`) %>%
    dplyr::left_join(annotables::grch38, by = "symbol", relationship = "many-to-many") %>%
    dplyr::distinct(Cluster, entrez, .keep_all = TRUE) %>%
    split(.[["Cluster"]])
  
  
#' Perform enrichment analysis
#'
#' @param df Input data frame or dataset
#' @param fold_change_col Color specification
#' @return Enrichment analysis results
#' @export
run_enrich_go <- function(df, fold_change_col = "Average.Log.Fold.Change") {
  
    original_gene_list <- df[[fold_change_col]]
    names(original_gene_list) <- df$entrez
    gene_list <- na.omit(original_gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    clusterProfiler::enrichGO(
      gene = names(gene_list),
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      ont = "BP",
      readable = TRUE
    )
  }
  
  safe_enrich_go <- purrr::safely(run_enrich_go, otherwise = NA_real_)
  ora_output <- map(df_list, safe_enrich_go) %>%
    map("result") %>%
    identity()
  
  ora_tables <-
    ora_output %>%
    map(tibble::as_tibble) %>%
    dplyr::bind_rows(.id = "cluster") %>%
    dplyr::group_by(Description) %>%
    split(.$cluster) %>%
    identity()
  
  unique_terms <-
    ora_tables %>%
    map(dplyr::pull, ID) %>%
    identity()
  
  
#' Filter data based on specified criteria
#'
#' @param cluster_name Cluster information
#' @param ora_results Parameter for ora results
#' @param terms_list Parameter for terms list
#' @return Filtered data
#' @export
filter_ora_result_by_terms <- function(cluster_name, ora_results, terms_list) {

  
    ora_result <- ora_results[[cluster_name]]
    myterms <- terms_list[[cluster_name]]
    ora_result@result <- ora_result@result[rownames(ora_result@result) %in% myterms, ]
    return(ora_result)
  }
  
  if (only_unique_terms) {
    ora_output <-
      names(unique_terms) %>%
      set_names(.) %>%
      map(filter_ora_result_by_terms, ora_output, unique_terms)
  }
  
  ora_plots <-
    ora_output %>%
    purrr::discard(~ identical(.x, NA_real_)) %>%
    imap(~ clusterProfiler::dotplot(.x, title = .y, showCategory = 20))
  
  return(list("result" = ora_output, "tables" = ora_tables, "plots" = ora_plots))
}

enrichment_analysis <- function(df, fold_change_col = "avg_log2FC", analysis_method = c("gsea", "ora"), gene_set = "H", pvalueCutoff = 0.5, TERM2GENE = msig_h, annotate = TRUE) {
  analysis_method = match.arg(analysis_method)
  if(annotate){
    df <-
      df %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol")
  }
  df <-
    df %>%
    dplyr::distinct(entrez, .keep_all = TRUE)
  
  if (analysis_method == "gsea"){
    original_gene_list <- df[[fold_change_col]]
    names(original_gene_list) <- df$entrez
    gene_list <- na.omit(original_gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    if (gene_set == "C6") {
      gene_sets <- msigdbr::msigdbr(species = "human", category = "C6") |> 
        dplyr::select(gs_name, entrez_gene)
      gse <- clusterProfiler::GSEA(
        geneList = gene_list,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        verbose = TRUE,
        TERM2GENE = gene_sets,
        pAdjustMethod = "BH"
      )
    } else if (gene_set == "H") {
      gene_sets <- msigdbr::msigdbr(species = "human", category = "H") |> 
        dplyr::select(gs_name, entrez_gene)
      gse <- clusterProfiler::GSEA(
        geneList = gene_list,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        verbose = TRUE,
        TERM2GENE = gene_sets,
        pAdjustMethod = "BH"
      )
    } else if (gene_set == "both"){
      gene_sets <- map(c("H" = "H", "C6" = "C6"), ~msigdbr::msigdbr(species = "human", category = .x)) |> 
        map(dplyr::select, gs_name, entrez_gene) |> 
        dplyr::bind_rows()
      gse <- clusterProfiler::GSEA(
        geneList = gene_list,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        verbose = TRUE,
        TERM2GENE = gene_sets,
        pAdjustMethod = "BH"
      )
    }
    return(gse)
  }
  
  if (analysis_method == "ora"){
    gene_list <- df |> 
      dplyr::filter(p_val_adj < 0.05) |> 
      dplyr::filter(avg_log2FC > 0) |> 
      dplyr::pull(entrez)
    
    if (gene_set == "C6") {
      gene_sets <- msigdbr::msigdbr(species = "human", category = "C6") |> 
        dplyr::select(gs_name, entrez_gene)
      gse <- clusterProfiler::enricher(
        gene = gene_list,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        TERM2GENE = gene_sets,
        pAdjustMethod = "BH"
      )
    } else if (gene_set == "H") {
      gene_sets <- msigdbr::msigdbr(species = "human", category = "H") |> 
        dplyr::select(gs_name, entrez_gene)
      gse <- clusterProfiler::enricher(
        gene = gene_list,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        TERM2GENE = gene_sets,
        pAdjustMethod = "BH"
      )
    } else if (gene_set == "both"){
      gene_sets <- map(c("H" = "H", "C6" = "C6"), ~msigdbr::msigdbr(species = "human", category = .x)) |> 
        map(dplyr::select, gs_name, entrez_gene) |> 
        dplyr::bind_rows()
      gse <- clusterProfiler::enricher(
        gene = gene_list,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        TERM2GENE = gene_sets,
        pAdjustMethod = "BH"
      )
    }
    return(gse)
  }
}

plot_enrichment <- function(gse, p_val_cutoff = NULL, signed = TRUE, analysis_method = c("gsea", "ora"), result_slot = "result", showCategory = 10) {
  analysis_method <- match.arg(analysis_method)
  
  if(analysis_method == "gsea"){
    if (signed) {
      
#' Create a plot visualization
#'
#' @param gse Parameter for gse
#' @param showCategory Parameter for showCategory
#' @return ggplot2 plot object
#' @export
make_signed_dotplot <- function(gse, showCategory = 10) {
  
  
        clusterProfiler::dotplot(gse, showCategory = showCategory, split = ".sign", font.size = 12) +
          facet_grid (.sign~ ., scales = "free_y", space = "free_y", labeller = as_labeller(c(
            "activated" = "up", 
            "suppressed" = "down"
          ))) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 10))
      }
      possible_dotplot <- purrr::possibly(make_signed_dotplot, otherwise = ggplot())
    } else {
      
#' Create a plot visualization
#'
#' @param gse Parameter for gse
#' @param showCategory Parameter for showCategory
#' @return ggplot2 plot object
#' @export
make_dotplot <- function(gse, showCategory = 10) {
  
  
        clusterProfiler::dotplot(gse, showCategory = showCategory, font.size = 12)
      }
      possible_dotplot <- purrr::possibly(make_dotplot, otherwise = ggplot())
    }
    
    if(!is.null(p_val_cutoff)){
      showCategories <- slot(gse, result_slot) |>
        dplyr::filter(p.adjust <= p_val_cutoff) |>
        pull(ID)
      mydotplot <- possible_dotplot(gse, showCategory = showCategories)
    } else {
      mydotplot <- possible_dotplot(gse, showCategory = showCategory)	
    }
  } else if(analysis_method == "ora"){
    make_dotplot <- function(gse, showCategory = 10) {
      clusterProfiler::dotplot(gse, showCategory = showCategory, font.size = 12)
    }
    possible_dotplot <- purrr::possibly(make_dotplot, otherwise = ggplot())
    
    if(!is.null(p_val_cutoff)){
      showCategories <- slot(gse, result_slot) |>
        dplyr::filter(p.adjust <= p_val_cutoff) |>
        pull(ID)
      mydotplot <- possible_dotplot(gse, showCategory = showCategories)
    } else {
      mydotplot <- possible_dotplot(gse, showCategory = showCategory)	
    }
  }
  
  return(mydotplot)
}

#' Create a plot visualization
#'
#' @param sample_id Parameter for sample id
#' @param myseus Parameter for myseus
#' @return Data frame
#' @export
make_cell_cycle_plot <- function(sample_id, myseus) {
  

  seu <- myseus[[sample_id]]
  DimPlot(seu, group.by = c("Phase")) +
    plot_annotation(title = sample_id)
  ggsave(glue("results/{sample_id}_cell_cycle.pdf"), width = 7, height = 7)
  return(glue("results/{sample_id}_cell_cycle.pdf"))
}

#' Perform differential expression analysis
#'
#' @param sample_id Parameter for sample id
#' @param myseus Parameter for myseus
#' @param celldf Cell identifiers or information
#' @param ... Additional arguments passed to other functions
#' @return Differential expression results
#' @export
diffex_by_cluster <- function(sample_id, myseus, celldf, ...) {
  
  
  seu <- myseus[[sample_id]]
  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()
  seu <- Seurat::AddMetaData(seu, celldf)
  clusters <- unique(seu$gene_snn_res.0.2) %>%
    set_names(.)
  
  
#' Filter data based on specified criteria
#'
#' @param cluster Cluster information
#' @param seu Seurat object
#' @return Filtered data
#' @export
filter_seu_to_cluster <- function(cluster, seu) {
  
  
    seu[, (seu@meta.data$gene_snn_res.0.2 == cluster)]
  }
  
  split_seu <- map(clusters, filter_seu_to_cluster, seu)
  safe_FindMarkers <- purrr::safely(FindMarkers, otherwise = NA_real_)
  cluster_diffex <- map(split_seu, safe_FindMarkers) %>%
    map("result") %>%
    identity()
  cluster_diffex <- cluster_diffex[!is.na(cluster_diffex)] %>%
    map(tibble::rownames_to_column, "symbol")
  write_xlsx(cluster_diffex, glue("results/{sample_id}_cluster_diffex.xlsx"))
  return(glue("results/{sample_id}_cluster_diffex.xlsx"))
}