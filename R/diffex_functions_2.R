# Diffex Functions (2)

#' Perform differential expression analysis
#'
#' @param per_cluster_clone_diffex Cluster information
#' @param total_clone_diffex Parameter for total clone diffex
#' @return Differential expression results
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - multiple_joins: Combine multiple joins into single join operation where possible

compare_per_cluster_and_total_clone_diffex <- function(per_cluster_clone_diffex, total_clone_diffex) {
  per_cluster_clone_diffex <-
    per_cluster_clone_diffex %>%
    readxl()

  total_clone_diffex <-
    total_clone_diffex %>%
    readxl()
}

#' Perform make clone comparison operation
#'
#' @param mysegs Parameter for mysegs
#' @param comparison Parameter for comparison
#' @param seu Seurat object
#' @param mynb Numbat object
#' @param location Character string (default: "cis")
#' @return Function result
#' @export
make_clone_comparison <- function(mysegs, comparison, seu, mynb, location = "cis") {
  #

  idents <-
    comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split(pattern = "_v_") %>%
    unlist()

  segments <- mynb$clone_post %>%
    dplyr::left_join(mynb$joint_post, by = "cell") %>%
    dplyr::filter(clone_opt %in% idents) %>%
    dplyr::filter(seg %in% mysegs) %>%
    dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
    dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
    dplyr::filter(!cnv_state_map == "neu") %>%
    plyranges::as_granges() %>%
    identity()
  
  # logfc.threshold = 0.1, test.use = "MAST"
  diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt", test.use = "wilcox", assay = "gene", logfc.threshold = 0.1) %>%
  	tibble::rownames_to_column("symbol") %>%
  	dplyr::left_join(annotables::grch38, by = "symbol") %>%
  	dplyr::distinct(ensgene, .keep_all = TRUE) %>%
  	dplyr::mutate(seqnames = chr) %>%
  	dplyr::filter(!is.na(start), !is.na(end)) %>%
  	plyranges::as_granges()
  
  # cis 
  
  cis_diffex <- 
  	diffex |> 
  plyranges::join_overlap_intersect(segments) %>%
  	as_tibble() %>%
  	dplyr::mutate(log2_sign = dplyr::case_when(
  		cnv_state_map == "amp" ~ -1,
  		cnv_state_map == "del" ~ 1
  	)) %>%
  	# dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
  	dplyr::select(-c(
  		"CHROM", "seg_start",
  		"seg_end", "cnv_state_map", "log2_sign"
  	)) %>%
  	dplyr::filter(!str_detect(chr, "CHR_")) %>%
  	dplyr::distinct(symbol, .keep_all = TRUE)
  
  # genes_in_segments <-
  # 	annotables::grch38 %>%
  # 	dplyr::rename(seqnames = chr) %>%
  # 	plyranges::as_granges() %>%
  # 	plyranges::join_overlap_intersect(segments) %>%
  # 	as_tibble() %>%
  # 	dplyr::select(seg, symbol) %>%
  # 	dplyr::group_by(seg) %>%
  # 	dplyr::summarize(genes_in_segment = paste(symbol, collapse = ", ")) %>%
  # 	identity()
  
  cis_diffex <-
  	cis_diffex %>%
  	# dplyr::left_join(genes_in_segments, by = "seg") |> 
  	append_clone_nums(comparison, seu)
  
  # trans
  
  trans_ranges <-
  	diffex %>%
  	plyranges::setdiff_ranges(segments)
  
  trans_diffex <-
  	diffex %>%
  	plyranges::join_overlap_intersect(trans_ranges) %>%
  	as_tibble() %>%
  	dplyr::filter(!str_detect(chr, "CHR_")) %>%
  	dplyr::distinct(symbol, .keep_all = TRUE) |> 
  	append_clone_nums(comparison, seu)
  
  # all 
  all_diffex <- dplyr::bind_rows(list("cis" = cis_diffex, "trans" = trans_diffex), .id = "location")
  
  if(location == "cis"){
  	return(cis_diffex)
  } else if (location == "trans"){
  	return(trans_diffex)
  } else if (location == "all"){
  	return(all_diffex)
  }


}

#' Perform clustering analysis
#'
#' @param mysegs Parameter for mysegs
#' @param clone_comparison Parameter for clone comparison
#' @param seu Seurat object
#' @param mynb Numbat object
#' @param w_scna Parameter for w scna
#' @param wo_scna Parameter for wo scna
#' @return List object
#' @export
make_cluster_comparison <- function(mysegs, clone_comparison, seu, mynb, w_scna = NULL, wo_scna = NULL) {
  #
  message(glue("{w_scna} v. {wo_scna}"))

  w_scna <- str_split_1(w_scna, "-")
  wo_scna <- str_split_1(wo_scna, "-")

  seu$scna_status <-
    seu@meta.data |>
    tibble::rownames_to_column("cell") |>
    dplyr::mutate(scna_status = dplyr::case_when(
      seu$clusters %in% w_scna ~ "w_scna",
      seu$clusters %in% wo_scna ~ "wo_scna"
    )) |>
    dplyr::pull(scna_status)
  
  
#' Extract or pull specific data elements
#'
#' @param mynb Numbat object
#' @param idents Cell identities or groups
#' @param mysegs Parameter for mysegs
#' @return Extracted data elements
#' @export
pull_segments <- function(mynb, idents, mysegs){
  	segments <- mynb$clone_post %>%
  		dplyr::left_join(mynb$joint_post, by = "cell") %>%
  		dplyr::filter(clone_opt %in% idents) %>%
  		dplyr::filter(seg %in% mysegs) %>%
  		dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
  		dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
  		dplyr::filter(!cnv_state_map == "neu") %>%
  		plyranges::as_granges() %>%
  		identity()
  }
  
  if(is.list(mynb)){
  	idents <-
  		map(clone_comparison, ~str_extract(.x, "[0-9]_v_[0-9]")) %>%
  		map(str_split, pattern = "_v_") %>%
  		map(unlist) |>
  		identity()
  	
  	segments <- pmap(list(mynb, idents, mysegs), pull_segments)
  	
  } else {
  	idents <-
  		clone_comparison %>%
  		str_extract("[0-9]_v_[0-9]") %>%
  		str_split(pattern = "_v_") %>%
  		unlist()
  	
  	segments <- pull_segments(mynb, idents, mysegs)
  	
  }

  diffex_all <- FindMarkers(seu, ident.1 = "w_scna", ident.2 = "wo_scna", group.by = "scna_status", logfc.threshold = 0.1, test.use = "MAST", assay = "gene") %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::distinct(ensgene, .keep_all = TRUE) %>%
    dplyr::mutate(seqnames = chr) %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::filter(!str_detect(chr, "CHR_")) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    # append_clone_nums(clone_comparison, seu) %>%
    # dplyr::mutate(to_SCT_snn_res. := {{to_SCT_snn_res.}},
    # 			  to_clust := paste({{to_clust}}, collapse = "_")) %>%
    dplyr::select(-strand)

  if(is.list(segments)){
  	segments <- unlist(as(segments, "GRangesList"))
  }
  
  diffex_cis <-
    diffex_all %>%
    plyranges::as_granges() %>%
    plyranges::join_overlap_intersect(segments) %>%
    as_tibble() %>%
    dplyr::mutate(log2_sign = dplyr::case_when(
      cnv_state_map == "amp" ~ -1,
      cnv_state_map == "del" ~ 1
    )) %>%
    # dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
    dplyr::select(-c(
      "CHROM", "seg_start",
      "seg_end", "cnv_state_map", "log2_sign"
    )) %>%
    dplyr::filter(!str_detect(chr, "CHR_")) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(-strand)

  trans_ranges <-
    diffex_all %>%
    plyranges::as_granges() %>%
    plyranges::setdiff_ranges(segments)

  diffex_trans <-
    diffex_all %>%
    plyranges::as_granges() %>%
    plyranges::join_overlap_intersect(trans_ranges) %>%
    as_tibble() %>%
    dplyr::filter(!str_detect(chr, "CHR_")) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(-strand)

  return(list("all" = diffex_all, "cis" = diffex_cis, "trans" = diffex_trans))
}
#' Perform make clone comparison integrated operation
#'
#' @param mysegs Parameter for mysegs
#' @param comparison Parameter for comparison
#' @param seu Seurat object
#' @param mynbs Parameter for mynbs
#' @param location Character string (default: "cis")
#' @return Function result
#' @export
make_clone_comparison_integrated <- function(mysegs, comparison, seu, mynbs, location = "cis") {
  #

  idents <-
    comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split(pattern = "_v_") %>%
    unlist()

  
#' Extract or pull specific data elements
#'
#' @param mynb Numbat object
#' @param idents Cell identities or groups
#' @return Extracted data elements
#' @export
pull_segment_from_nb <- function(mynb, idents) {
    segments <- mynb$clone_post %>%
      dplyr::left_join(mynb$joint_post, by = "cell") %>%
      dplyr::filter(clone_opt %in% idents) %>%
      dplyr::filter(seg %in% mysegs) %>%
      dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
      dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
      dplyr::filter(!cnv_state_map == "neu") %>%
      plyranges::as_granges() %>%
      identity()
  }

  segments <- map(mynbs, pull_segment_from_nb, idents)

  segments <- unlist(as(segments, "GRangesList"))

  if (location == "cis") {
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt", logfc.threshold = 0.1) %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      plyranges::as_granges() %>%
      plyranges::join_overlap_intersect(segments) %>%
      as_tibble() %>%
      dplyr::mutate(log2_sign = dplyr::case_when(
        cnv_state_map == "amp" ~ -1,
        cnv_state_map == "del" ~ 1
      )) %>%
      # dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
      dplyr::select(-c(
        "CHROM", "seg_start",
        "seg_end", "cnv_state_map", "log2_sign"
      )) %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)
  } else if (location == "trans") {
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt", logfc.threshold = 0.1) %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      plyranges::as_granges()

    out_of_segment_ranges <-
      diffex %>%
      plyranges::setdiff_ranges(segments)

    diffex <-
      diffex %>%
      plyranges::join_overlap_intersect(out_of_segment_ranges) %>%
      as_tibble() %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)
  }
}

