# Automated cluster-dictionary annotation.
#
# Replaces the hand-curated data/cluster_dictionary.tsv: each gene_snn_res.0.2
# cluster is labelled from its real marker genes (top-20 stored in the
# cluster_markers table of batch_hashes.sqlite, or computed from the seu when a
# sample is missing from the DB).
#
# Rule precedence (high -> low) in classify_cluster_abbreviation():
#   1. low_qual : rank-1 marker is mitochondrial (MT-*)
#   2. G2M      : >= 3 of top-`top_n` markers are cell-cycle genes
#   3. HSP      : >= 2 of top-`top_n` markers are heat-shock genes
#   4. cone     : >= 2 of top-`top_n` markers are cone-specific genes
#   5. rod      : >= 2 of top-`top_n` markers are rod-specific genes
#   6. MALAT1   : MALAT1 is among the top-`malat1_rank_max` markers
#   7. <gene>   : otherwise the rank-1 marker gene verbatim

# ---- curated program gene sets ------------------------------------------------
.cluster_program_gene_sets <- list(
  cell_cycle = c("CENPF","TOP2A","UBE2C","MKI67","ASPM","HMGB2","CKS1B","PTTG1",
    "NUF2","SMC4","BIRC5","KPNA2","NUSAP1","CDK1","CCNB1","CCNB2","CDC20","CENPA",
    "CENPE","CKAP2","CKAP2L","AURKA","AURKB","BUB1","BUB1B","CDCA3","CDCA8","ANLN",
    "TPX2","HMMR","GTSE1","ECT2","KIF11","KIF23","KIF2C","NCAPD2","NCAPG","PLK1",
    "RRM2","TYMS","CCNA2","CDCA2","DLGAP5","KIF20A","PRC1","RACGAP1","TROAP",
    "NDC80","HJURP","G2E3","CDC25C","TUBB4B","HMGB3"),
  heat_shock = c("HSPA1A","HSPA1B","HSPA6","HSPA8","HSPA4","HSPB1","HSPH1",
    "HSP90AA1","HSP90AB1","DNAJB1","DNAJA1","DNAJB4","BAG3","HSPE1","HSPD1",
    "AHSA1","CHORDC1","FKBP4","SERPINH1","ZFAND2A"),
  cone = c("ARR3","GNGT2","OPN1SW","OPN1MW","OPN1LW","PDE6H","PDE6C","GNAT2",
    "GUCA1C","GNB3","PCP2","RXRG","CNGB3","CNGA3"),
  rod = c("RHO","NRL","NR2E3","PDE6A","PDE6B","GNAT1","CNGA1","CNGB1","SAG",
    "RP1","RDH12","GNGT1")
)

#' Classify a single cluster's abbreviation from its ranked marker genes
#'
#' @param top_genes Character vector of marker gene symbols, ordered by
#'   marker rank (rank 1 first).
#' @param malat1_rank_max Integer; MALAT1 must appear at or above this rank to
#'   label the cluster "MALAT1" (default 5).
#' @param top_n Integer window used for the program (G2M/HSP/cone/rod) rules
#'   (default 10).
#' @return A list with `abbrev`, `rule`, and `matched` (comma-separated
#'   matched genes).
#' @export
classify_cluster_abbreviation <- function(top_genes, malat1_rank_max = 5, top_n = 10) {
  top_genes <- top_genes[!is.na(top_genes)]
  if (length(top_genes) == 0) {
    return(list(abbrev = NA_character_, rule = "no_markers", matched = ""))
  }
  win  <- utils::head(top_genes, top_n)
  sets <- .cluster_program_gene_sets

  if (grepl("^MT-", top_genes[1])) {
    return(list(abbrev = "low_qual", rule = "mito_top1", matched = top_genes[1]))
  }
  hit <- function(set) win[win %in% set]
  if (length(hit(sets$cell_cycle)) >= 3) {
    return(list(abbrev = "G2M", rule = "cell_cycle",
                matched = paste(hit(sets$cell_cycle), collapse = ",")))
  }
  if (length(hit(sets$heat_shock)) >= 2) {
    return(list(abbrev = "HSP", rule = "heat_shock",
                matched = paste(hit(sets$heat_shock), collapse = ",")))
  }
  if (length(hit(sets$cone)) >= 2) {
    return(list(abbrev = "cone", rule = "cone",
                matched = paste(hit(sets$cone), collapse = ",")))
  }
  if (length(hit(sets$rod)) >= 2) {
    return(list(abbrev = "rod", rule = "rod",
                matched = paste(hit(sets$rod), collapse = ",")))
  }
  if ("MALAT1" %in% utils::head(top_genes, malat1_rank_max)) {
    return(list(abbrev = "MALAT1", rule = "malat1_top_rank", matched = "MALAT1"))
  }
  list(abbrev = top_genes[1], rule = "top1_gene", matched = top_genes[1])
}

# Pull gene_snn_res.0.2 top markers per cluster for one sample from the DB.
.dict_markers_from_db <- function(con, filepath, n_keep = 20) {
  q <- DBI::dbGetQuery(con,
    "SELECT cluster, marker_rank, gene_name
       FROM cluster_markers
      WHERE filepath = ? AND marker_rank <= ?
      ORDER BY cluster, marker_rank",
    params = list(filepath, n_keep))
  if (nrow(q) == 0) return(NULL)
  q
}

# Compute gene_snn_res.0.2 top markers for one sample from its seu (fallback).
.dict_markers_from_seu <- function(seu_path, n_keep = 20) {
  if (!file.exists(seu_path)) return(NULL)
  seu <- readRDS(seu_path)
  if (!"gene_snn_res.0.2" %in% colnames(seu@meta.data)) { rm(seu); gc(); return(NULL) }
  Seurat::DefaultAssay(seu) <- "gene"
  if (inherits(seu[["gene"]], "Assay5")) {
    seu[["gene"]] <- SeuratObject::JoinLayers(seu[["gene"]])
  }
  dat <- tryCatch(SeuratObject::LayerData(seu, assay = "gene", layer = "data"),
                  error = function(e) NULL)
  if (is.null(dat) || length(dat@x) == 0) {
    seu <- Seurat::NormalizeData(seu, assay = "gene", verbose = FALSE)
  }
  m <- tryCatch(presto::wilcoxauc(seu, "gene_snn_res.0.2", seurat_assay = "gene"),
                error = function(e) NULL)
  rm(seu); gc()
  if (is.null(m)) return(NULL)
  m |>
    dplyr::group_by(.data$group) |>
    dplyr::filter(.data$padj < 0.5) |>
    dplyr::arrange(.data$group, dplyr::desc(.data$logFC)) |>
    dplyr::slice_head(n = n_keep) |>
    dplyr::mutate(marker_rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::transmute(cluster = as.character(.data$group),
                     marker_rank = .data$marker_rank,
                     gene_name = .data$feature)
}

#' Build the cluster dictionary from marker genes
#'
#' Drop-in replacement for `read_cluster_dictionary()`: returns a named list of
#' tibbles (one per sample_id) with columns `sample_id`, `abbreviation`,
#' `gene_snn_res.0.2`, `remove`.
#'
#' @param seu_files Named (by sample id) character vector of base unfiltered
#'   seurat .rds paths (e.g. the `all_seu_files` target).
#' @param sqlite_path Path to batch_hashes.sqlite.
#' @param malat1_rank_max,top_n Passed to [classify_cluster_abbreviation()].
#' @param remove_programs Abbreviations whose clusters get `remove = "1"`.
#' @param write_audit_dir If non-NULL, write `cluster_dictionary_auto.tsv` and
#'   a discrepancy report (vs `compare_to`) into this directory.
#' @param compare_to Optional path to the previous hand dictionary tsv for the
#'   discrepancy report.
#' @param .db_ready Ignored; declare a dependency on the metadata-DB target so
#'   targets runs DB population first.
#' @return Named list of per-sample dictionary tibbles.
#' @export
build_cluster_dictionary <- function(seu_files,
                                     sqlite_path = "batch_hashes.sqlite",
                                     malat1_rank_max = 5,
                                     top_n = 10,
                                     remove_programs = c("rod", "cone", "APOE", "low_qual", "MALAT1"),
                                     write_audit_dir = NULL,
                                     compare_to = NULL,
                                     .db_ready = NULL) {
  if (is.null(names(seu_files))) {
    names(seu_files) <- stringr::str_extract(seu_files, "SR[RX][0-9]+")
  }

  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  rows <- list()
  for (sid in names(seu_files)) {
    fp <- unname(seu_files[[sid]])
    mk <- .dict_markers_from_db(con, fp)
    if (is.null(mk)) {
      message("cluster_dictionary: computing markers from seu for ", sid)
      mk <- .dict_markers_from_seu(fp)
    }
    if (is.null(mk)) { message("cluster_dictionary: no markers for ", sid); next }

    clusters <- mk |>
      dplyr::distinct(.data$cluster) |>
      dplyr::mutate(cl_int = suppressWarnings(as.integer(.data$cluster))) |>
      dplyr::arrange(.data$cl_int) |>
      dplyr::pull(.data$cluster)

    for (cl in clusters) {
      g <- mk |>
        dplyr::filter(.data$cluster == cl) |>
        dplyr::arrange(.data$marker_rank) |>
        dplyr::pull(.data$gene_name)
      cls <- classify_cluster_abbreviation(g, malat1_rank_max = malat1_rank_max, top_n = top_n)
      rows[[length(rows) + 1]] <- tibble::tibble(
        sample_id        = sid,
        gene_snn_res.0.2 = suppressWarnings(as.integer(cl)),
        abbreviation     = cls$abbrev,
        remove           = if (!is.na(cls$abbrev) && cls$abbrev %in% remove_programs) "1" else NA_character_,
        rule             = cls$rule,
        matched          = cls$matched,
        top_genes        = paste(utils::head(g, top_n), collapse = ", ")
      )
    }
  }
  auto <- dplyr::bind_rows(rows) |>
    dplyr::arrange(.data$sample_id, .data$gene_snn_res.0.2)

  if (!is.null(write_audit_dir)) {
    fs::dir_create(write_audit_dir)
    auto |>
      dplyr::transmute(.data$sample_id, .data$abbreviation, .data$gene_snn_res.0.2, .data$remove) |>
      readr::write_tsv(file.path(write_audit_dir, "cluster_dictionary_auto.tsv"), na = "NA")
    if (!is.null(compare_to) && file.exists(compare_to)) {
      old <- readr::read_tsv(compare_to, show_col_types = FALSE) |>
        dplyr::mutate(gene_snn_res.0.2 = as.integer(.data$gene_snn_res.0.2)) |>
        dplyr::select(.data$sample_id, .data$gene_snn_res.0.2,
                      old_abbreviation = .data$abbreviation, old_remove = .data$remove)
      rep <- auto |>
        dplyr::left_join(old, by = c("sample_id", "gene_snn_res.0.2")) |>
        dplyr::mutate(changed = is.na(.data$old_abbreviation) | .data$old_abbreviation != .data$abbreviation) |>
        dplyr::select(.data$sample_id, .data$gene_snn_res.0.2, .data$old_abbreviation,
                      .data$abbreviation, .data$changed, .data$rule, .data$remove,
                      .data$matched, .data$top_genes)
      readr::write_tsv(rep, file.path(write_audit_dir, "cluster_dictionary_discrepancies.tsv"), na = "NA")
    }
  }

  out <- auto |>
    dplyr::select(.data$sample_id, .data$abbreviation, .data$gene_snn_res.0.2, .data$remove)
  split(out, out$sample_id)
}
