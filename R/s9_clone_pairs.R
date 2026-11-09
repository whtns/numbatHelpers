# Fig. S9 clone pairs: within-cluster diffex and cross-tumor clone-pair integrations.
#
# An S9 pair is one row of results/fig_s9_t1e5/fig_s9_pages.csv: the clone that acquired an
# RB SCNA (`clone_with`) vs its direct parent (`clone_without`), at one numbat t=1e-5
# iteration. Membership comes from output/numbat_t1e5/<sample>/clone_post_<iteration>.tsv,
# not from the objects' clone_opt (that is the final round only), as in
# src/fig_s9_clone_diffex.R, whose helpers are ported here.
#
# Targets: R/pipeline_targets_s9.R. Method notes: docs/fig_s9_within_cluster_and_integrations.md.

.S9_NONE_TITLE <- "No eligible RB SCNA clone comparison"

.S9_RB_OF <- c("1q+" = "1q+", "1+" = "1q+", "2p+" = "2p+", "2+" = "2p+", "6p+" = "6p+", "6+" = "6p+",
               "16q-" = "16q-", "16-" = "16q-")

# hg38 centromere midpoints (bp), as in src/plot_clone_tree_collage_t1e5.R
.S9_CEN <- c("1" = 123.5, "2" = 93.1, "3" = 92.2, "4" = 50.7, "5" = 48.3, "6" = 59.2, "7" = 59.8,
             "8" = 45.0, "9" = 44.4, "10" = 40.6, "11" = 52.7, "12" = 35.6, "13" = 17.0, "14" = 17.1,
             "15" = 18.4, "16" = 37.3, "17" = 24.9, "18" = 18.2, "19" = 25.9, "20" = 28.2,
             "21" = 11.9, "22" = 14.0) * 1e6

.S9_SCNA_TOKEN <- c("1q+" = "1q", "2p+" = "2p", "6p+" = "6p", "16q-" = "16q")

#' File-name slug of an S9 page title (same as the S9 builder, so names line up)
#' @param x Page title(s)
#' @return Character
#' @export
s9_slug <- function(x) {
  x |>
    stringr::str_replace_all("[+−-]", "") |>
    stringr::str_replace_all("[(),]", "") |>
    stringr::str_squish() |>
    stringr::str_replace_all(" ", "_") |>
    stringr::str_replace_all("co_acquired|coacquired", "coacq") |>
    tolower()
}

#' Read the Fig. S9 clone pairs
#'
#' @param pages_csv results/fig_s9_t1e5/fig_s9_pages.csv
#' @param seu_dir Directory of the <sample>_filtered_seu.rds objects
#' @param nb_dir numbat t=1e-5 output directory
#' @return One row per S9 panel, with clone_with, clone_without, key, n_scna, seu_path,
#'   clone_post and joint_post paths
#' @export
read_s9_pairs <- function(pages_csv, seu_dir = "output/seurat", nb_dir = "output/numbat_t1e5") {
  readr::read_csv(pages_csv, show_col_types = FALSE) |>
    dplyr::filter(page != .S9_NONE_TITLE, !is.na(pair)) |>
    dplyr::mutate(
      iteration = as.integer(iteration),
      clone_with = as.integer(stringr::str_extract(pair, "^[0-9]+")),
      clone_without = as.integer(stringr::str_extract(pair, "[0-9]+$")),
      key = sprintf("%s_it%d_%s", sample, iteration, pair),
      n_scna = lengths(stringr::str_split(acquired, ";")),
      seu_path = file.path(seu_dir, paste0(sample, "_filtered_seu.rds")),
      clone_post = file.path(nb_dir, sample, sprintf("clone_post_%d.tsv", iteration)),
      joint_post = file.path(nb_dir, sample, sprintf("joint_post_%d.tsv", iteration))
    ) |>
    (\(d) {
      no_seu <- !file.exists(d$seu_path)
      if (any(no_seu)) message("no filtered object, left out: ", paste(unique(d$sample[no_seu]), collapse = ", "))
      d[!no_seu, ]
    })()
}

#' Input files of the S9 pairs (for file tracking)
#' @param pairs Output of read_s9_pairs()
#' @return Unique paths: filtered objects, clone_post and joint_post files
#' @export
s9_input_paths <- function(pairs) {
  p <- unique(c(pairs$seu_path, pairs$clone_post, pairs$joint_post))
  missing <- p[!file.exists(p)]
  if (length(missing)) stop("missing S9 inputs: ", paste(missing, collapse = ", "))
  p
}

# A filtered object ready for S9 tests: assay gene, normalized, filter_keep FALSE dropped
.s9_read_seu <- function(seu_path) {
  sid <- stringr::str_extract(basename(seu_path), "SR[RX][0-9]+")
  seu <- readRDS(seu_path)
  Seurat::DefaultAssay(seu) <- "gene"
  if (inherits(seu[["gene"]], "Assay5") && length(SeuratObject::Layers(seu[["gene"]], search = "counts")) > 1)
    seu[["gene"]] <- SeuratObject::JoinLayers(seu[["gene"]])
  if (!"data" %in% SeuratObject::Layers(seu[["gene"]]) ||
      sum(SeuratObject::LayerData(seu, assay = "gene", layer = "data")[, 1:min(50, ncol(seu))]) == 0) {
    message(sid, ": gene assay has no normalized data; running NormalizeData")
    seu <- Seurat::NormalizeData(seu, assay = "gene", verbose = FALSE)
  }
  if ("filter_keep" %in% colnames(seu@meta.data) && is.logical(seu$filter_keep) &&
      any(!seu$filter_keep, na.rm = TRUE)) {
    message(sid, ": dropping ", sum(!seu$filter_keep, na.rm = TRUE), " cells with filter_keep FALSE")
    seu <- seu[, which(seu$filter_keep)]
  }
  seu
}

#' Assign cells to the two clones of an S9 pair
#'
#' @param cells Cell barcodes (colnames of the filtered object)
#' @param clone_post clone_post_<iteration>.tsv of the pair's iteration
#' @param with,without Clone numbers
#' @return Character vector along `cells`: "with", "without" or NA
#' @export
s9_cell_groups <- function(cells, clone_post, with, without) {
  cp <- readr::read_tsv(clone_post, show_col_types = FALSE, col_select = c("cell", "clone_opt"))
  grp <- cp$clone_opt[match(cells, cp$cell)]
  dplyr::case_when(grp == with ~ "with", grp == without ~ "without")
}

#' Segments that differ between the two clones of an S9 pair, with arm-level labels
#'
#' Port of differing_segments() in src/fig_s9_clone_diffex.R.
#' @param clone_post,joint_post numbat files of the pair's iteration
#' @param with,without Clone numbers
#' @return tibble seg, CHROM, seg_start, seg_end, cnv_state, label ("1q+", "7+"), rb
#' @export
s9_differing_segments <- function(clone_post, joint_post, with, without) {
  cp <- readr::read_tsv(clone_post, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  gt <- function(k) {
    g <- unique(cp$GT_opt[cp$clone_opt == as.character(k)])
    stopifnot(length(g) == 1)
    s <- stringr::str_split(tidyr::replace_na(g, ""), ",")[[1]]
    s[s != ""]
  }
  segs <- union(setdiff(gt(with), gt(without)), setdiff(gt(without), gt(with)))
  jp <- data.table::fread(joint_post, select = c("seg", "CHROM", "seg_start", "seg_end", "cnv_state")) |>
    tibble::as_tibble() |>
    dplyr::distinct(seg, .keep_all = TRUE)
  jp |>
    dplyr::filter(seg %in% segs) |>
    dplyr::mutate(
      cen = .S9_CEN[as.character(CHROM)],
      arm = dplyr::case_when(seg_start < cen - 5e6 & seg_end > cen + 5e6 ~ "",
                             (seg_start + seg_end) / 2 < cen ~ "p", TRUE ~ "q"),
      sign = dplyr::case_when(cnv_state %in% c("amp", "bamp") ~ "+",
                              cnv_state %in% c("del", "bdel") ~ "-", TRUE ~ paste0(" ", cnv_state)),
      label = paste0(CHROM, arm, sign),
      rb = unname(.S9_RB_OF[label])
    ) |>
    dplyr::select(seg, CHROM, seg_start, seg_end, cnv_state, label, rb)
}

#' Locate genes relative to an S9 pair's differing segments
#'
#' Port of locate_genes() in src/fig_s9_clone_diffex.R.
#' @param symbols Gene symbols
#' @param segs Output of s9_differing_segments()
#' @param acquired Acquired RB SCNAs, e.g. c("1q+", "16q-")
#' @return tibble symbol, chr, start, end, location (cis / edge_other / trans)
#' @export
s9_locate_genes <- function(symbols, segs, acquired) {
  genes <- annotables::grch38 |>
    dplyr::filter(symbol %in% symbols, chr %in% c(1:22, "X", "Y")) |>
    dplyr::distinct(symbol, .keep_all = TRUE) |>
    dplyr::select(symbol, chr, start, end)
  hit <- function(s) {
    if (nrow(s) == 0) return(character())
    genes |>
      dplyr::inner_join(s, by = c(chr = "CHROM"), relationship = "many-to-many") |>
      dplyr::filter(start <= seg_end, end >= seg_start) |>
      dplyr::pull(symbol) |>
      unique()
  }
  segs <- dplyr::mutate(segs, CHROM = as.character(CHROM))
  cis <- hit(dplyr::filter(segs, rb %in% acquired))
  other <- setdiff(hit(dplyr::filter(segs, !rb %in% acquired | is.na(rb))), cis)
  tibble::tibble(symbol = symbols) |>
    dplyr::left_join(genes, by = "symbol") |>
    dplyr::mutate(location = dplyr::case_when(symbol %in% cis ~ "cis", symbol %in% other ~ "edge_other",
                                              TRUE ~ "trans"))
}

# Wilcoxon with vs without on the cells of `seu` (assay gene), as the whole-sample S9 test
.s9_wilcox <- function(seu, group_col = "s9_group") {
  Seurat::FindMarkers(seu, ident.1 = "with", ident.2 = "without", group.by = group_col,
                      assay = "gene", test.use = "wilcox", logfc.threshold = 0, min.pct = 0.05,
                      verbose = FALSE) |>
    tibble::rownames_to_column("symbol") |>
    tibble::as_tibble()
}

# Cells per cluster x group, and whether the cluster is tested
.s9_counts <- function(md, by, min_cells) {
  md |>
    dplyr::filter(!is.na(s9_group)) |>
    dplyr::count(dplyr::across(dplyr::all_of(by)), s9_group) |>
    tidyr::pivot_wider(names_from = s9_group, values_from = n, values_fill = 0) |>
    (\(d) { for (g in c("with", "without")) if (!g %in% names(d)) d[[g]] <- 0L; d })() |>
    dplyr::rename(cells_with = with, cells_without = without) |>
    dplyr::mutate(tested = cells_with >= min_cells & cells_without >= min_cells)
}

#' Within-cluster diffex of every S9 pair of one tumor
#'
#' Clusters are `resolution` of the filtered object, phase-labelled with
#' assign_auto_phase_clusters(). In each cluster where both clones have at least
#' `min_cells` cells, the clone with the SCNA is tested against its comparison clone
#' (Wilcoxon, assay gene, logfc.threshold 0, min.pct 0.05, as the whole-sample S9 diffex).
#'
#' @param pairs Rows of read_s9_pairs() for one sample
#' @param resolution Cluster column (default SCT_snn_res.0.4)
#' @param min_cells Minimum cells per clone in a cluster
#' @param out_dir Output directory
#' @return Paths: per pair, the diffex CSV (<key>.csv) and the cell-count CSV (<key>_counts.csv)
#' @export
s9_within_cluster_diffex <- function(pairs, resolution = "SCT_snn_res.0.4", min_cells = 20,
                                     out_dir = "results/fig_s9_t1e5/cluster_diffex") {
  sid <- unique(pairs$sample)
  stopifnot(length(sid) == 1)
  pairs <- dplyr::distinct(pairs, key, .keep_all = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seu <- .s9_read_seu(pairs$seu_path[1])
  if (!resolution %in% colnames(seu@meta.data)) stop(sid, ": no ", resolution, " in ", pairs$seu_path[1])
  # the filtered objects carry no Phase / S.Score (filtered_seus_with_phase is a separate
  # target); without them every cluster is labelled "other"
  if (!all(c("Phase", "S.Score", "G2M.Score") %in% colnames(seu@meta.data)) || all(is.na(seu$Phase))) {
    message(sid, ": scoring cell cycle on the gene assay")
    seu <- Seurat::CellCycleScoring(seu, s.features = Seurat::cc.genes.updated.2019$s.genes,
                                    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
  }
  seu <- assign_auto_phase_clusters(seu, resolution)
  paths <- character()
  for (i in seq_len(nrow(pairs))) {
    r <- pairs[i, ]
    seu$s9_group <- s9_cell_groups(colnames(seu), r$clone_post, r$clone_with, r$clone_without)
    counts <- .s9_counts(seu@meta.data, c("clusters", "phase_level"), min_cells) |>
      dplyr::mutate(sample = sid, key = r$key, page = r$page, resolution = resolution,
                    error = NA_character_)
    message(sprintf("%s: %d with, %d without; %d of %d clusters tested", r$key,
                    sum(counts$cells_with), sum(counts$cells_without), sum(counts$tested), nrow(counts)))
    segs <- s9_differing_segments(r$clone_post, r$joint_post, r$clone_with, r$clone_without)
    acq <- stringr::str_split(r$acquired, ";")[[1]]
    if (!any(segs$rb %in% acq)) message("  WARNING: no differing segment maps to ", r$acquired)
    de <- list()
    for (cl in as.character(counts$clusters[counts$tested])) {
      res <- tryCatch(.s9_wilcox(seu[, which(seu$clusters == cl & !is.na(seu$s9_group))]),
                      error = function(e) conditionMessage(e))
      if (is.character(res)) {
        message("  ", cl, ": FindMarkers failed: ", res)
        counts$error[as.character(counts$clusters) == cl] <- res
        next
      }
      de[[cl]] <- res
    }
    de <- dplyr::bind_rows(de, .id = "clusters")
    if (nrow(de) > 0) {
      de <- de |>
        dplyr::left_join(s9_locate_genes(unique(de$symbol), segs, acq), by = "symbol") |>
        dplyr::left_join(dplyr::mutate(counts, clusters = as.character(clusters)) |>
                           dplyr::select(clusters, phase_level, cells_with, cells_without), by = "clusters")
    }
    de <- dplyr::mutate(de, sample = sid, iteration = r$iteration, pair = r$pair, acquired = r$acquired,
                        page = r$page, key = r$key, resolution = resolution,
                        cis_segments = paste(segs$seg[segs$rb %in% acq], collapse = ","),
                        edge_other = paste(segs$label[!segs$rb %in% acq | is.na(segs$rb)], collapse = ","))
    p_de <- file.path(out_dir, paste0(r$key, ".csv"))
    p_n <- file.path(out_dir, paste0(r$key, "_counts.csv"))
    readr::write_csv(de, p_de)
    readr::write_csv(counts, p_n)
    paths <- c(paths, p_de, p_n)
  }
  paths
}

#' Integration groups of S9 pairs
#'
#' Two kinds, one pair per tumor in each:
#' - `page_<slug>`: an S9 page with at least `min_tumors` tumors;
#' - `alone_<scna>`: every pair that acquires that RB SCNA without another RB SCNA,
#'   whatever the background.
#' A tumor with several pairs in a group keeps the one with adjusted == FALSE, then the
#' later iteration. `expect_iteration` asserts that choice for named samples.
#'
#' @param pairs Output of read_s9_pairs()
#' @param min_tumors Minimum tumors per group
#' @param expect_iteration Named integer: sample -> iteration it must resolve to
#' @return Pair rows with group_id, group_label, group_kind
#' @export
make_s9_integration_groups <- function(pairs, min_tumors = 2, expect_iteration = c(SRX11133592 = 2L)) {
  pick_one <- function(d) {
    d |>
      dplyr::group_by(group_id, sample) |>
      dplyr::arrange(adjusted, dplyr::desc(iteration), .by_group = TRUE) |>
      dplyr::slice(1) |>
      dplyr::ungroup()
  }
  pages <- pairs |>
    dplyr::mutate(group_id = paste0("page_", s9_slug(page)), group_label = page, group_kind = "page")
  alone <- purrr::imap_dfr(.S9_SCNA_TOKEN, function(tok, scna) {
    pairs |>
      dplyr::filter(acquired == scna) |>
      dplyr::mutate(group_id = paste0("alone_", tok),
                    group_label = paste0(sub("-$", "−", scna), " alone (not co-acquired with another RB SCNA)"),
                    group_kind = "alone")
  })
  groups <- dplyr::bind_rows(pages, alone) |>
    pick_one() |>
    dplyr::group_by(group_id) |>
    dplyr::filter(dplyr::n_distinct(sample) >= min_tumors) |>
    dplyr::ungroup() |>
    dplyr::arrange(group_kind != "page", group_id, sample)
  stopifnot(!anyDuplicated(groups[c("group_id", "sample")]))
  for (s in names(expect_iteration)) {
    got <- groups$iteration[groups$sample == s]
    if (length(got) && any(got != expect_iteration[[s]]))
      stop(s, " resolved to iteration ", paste(unique(got), collapse = "/"), ", expected ", expect_iteration[[s]])
  }
  groups
}

# CCA integration as seuratTools::seurat_integrate(): per-batch preprocessing, small batches
# merged, FindIntegrationAnchors(dims 1:30, k.filter 50), IntegrateData at k.weight =
# min(smallest batch - 15, 100). One difference: when IntegrateData fails (too few anchor cells
# in a small tumor; e.g. page_1q_over_2p_6p, 113 anchors, SRX10031194 375 cells) it retries at
# a lower k.weight instead of falling back to harmony (not installed). The first attempt is
# identical, so integrations that succeed in seuratTools come out the same.
.s9_integrate_cca <- function(seus, k_weights = c(50, 30, 20)) {
  seus <- purrr::imap(seus, function(s, nm) {
    s[["gene"]] <- seuratTools::seurat_preprocess(s[["gene"]], scale = TRUE)
    s$batch <- nm
    s
  })
  seus <- seuratTools::merge_small_seus(seus)
  if (length(seus) < 2) stop("fewer than 2 batches left after merging small ones")
  anchors <- Seurat::FindIntegrationAnchors(object.list = seus, dims = 1:30, k.filter = 50)
  kw0 <- min(min(sapply(seus, ncol)) - 15, 100)
  res <- NULL
  for (kw in unique(c(kw0, k_weights[k_weights < kw0]))) {
    res <- tryCatch(Seurat::IntegrateData(anchorset = anchors, dims = 1:30, k.weight = kw),
                    error = function(e) e)
    if (!inherits(res, "error")) break
    message("IntegrateData failed at k.weight ", kw, ": ", conditionMessage(res))
  }
  if (inherits(res, "error")) stop("IntegrateData failed at every k.weight: ", conditionMessage(res))
  if (kw != kw0) message("integrated at k.weight ", kw, " (default ", kw0, " failed)")
  res@misc$s9_k_weight <- kw
  Seurat::Idents(res) <- "batch"
  res[["batch"]] <- Seurat::Idents(res)
  Seurat::DefaultAssay(res) <- "integrated"
  res
}

#' Integrate the S9 clone pairs of one group
#'
#' Each tumor contributes only the cells of its pair's two clones (s9_group = with / without).
#' Same engine as integration_by_scna_clones(): seuratTools::integration_workflow(), then
#' ScaleData, dimension reduction and clustering at 0.2-2.0 on the integrated assay. The gene
#' assay is joined and cell cycle rescored (as annotate_integrated_states()), and `clusters` /
#' `phase_level` come from assign_auto_phase_clusters() at `resolution`.
#'
#' @param group_tbl Rows of make_s9_integration_groups() for one group
#' @param resolution Integrated cluster column to phase-label
#' @param out_dir Output directory
#' @return Path to <group_id>_filtered_seu.rds
#' @export
integrate_s9_group <- function(group_tbl, resolution = "integrated_snn_res.0.4", min_cells = 20,
                               out_dir = "output/seurat/integrated_s9") {
  gid <- unique(group_tbl$group_id)
  stopifnot(length(gid) == 1, !anyDuplicated(group_tbl$sample))
  message(gid, ": ", paste(sprintf("%s %s", group_tbl$sample, group_tbl$key), collapse = ", "))
  seus <- purrr::map(seq_len(nrow(group_tbl)), function(i) {
    r <- group_tbl[i, ]
    seu <- .s9_read_seu(r$seu_path)
    grp <- s9_cell_groups(colnames(seu), r$clone_post, r$clone_with, r$clone_without)
    seu <- seu[, which(!is.na(grp))]
    seu$s9_group <- grp[!is.na(grp)]
    seu$s9_tumor <- r$sample
    seu$s9_key <- r$key
    seu$s9_cell <- colnames(seu)
    drop <- grepl("_snn_res\\.", colnames(seu@meta.data)) |
      colnames(seu@meta.data) %in% c("clusters", "phase_level", "seurat_clusters")
    seu@meta.data <- seu@meta.data[, !drop, drop = FALSE]
    n <- c(sum(seu$s9_group == "with"), sum(seu$s9_group == "without"))
    message(sprintf("  %s: %d with, %d without", r$key, n[1], n[2]))
    gc()
    # a pair with almost no cells in one clone breaks IntegrateData (and adds nothing)
    if (any(n < min_cells)) { message("  dropped: fewer than ", min_cells, " cells in a clone"); return(NULL) }
    seu
  }) |>
    rlang::set_names(group_tbl$sample)
  dropped <- names(seus)[purrr::map_lgl(seus, is.null)]
  seus <- purrr::compact(seus)
  if (length(seus) < 2) stop(gid, ": fewer than 2 tumors with >= ", min_cells, " cells per clone")

  seu <- .s9_integrate_cca(seus)
  rm(seus); gc()
  seu <- Seurat::ScaleData(seu)
  seu <- seuratTools::seurat_reduce_dimensions(seu)
  seu <- seuratTools::seurat_cluster(seu, resolution = seq(0.2, 2.0, by = 0.2), seurat_assay = "integrated")

  # v5 assays come out of integration split one layer per batch; a v3 Assay (e.g.
  # page_2p_vs_diploid's inputs) has no layers to join
  if (inherits(seu[["gene"]], "Assay5")) seu[["gene"]] <- SeuratObject::JoinLayers(seu[["gene"]])
  default_assay <- Seurat::DefaultAssay(seu)
  Seurat::DefaultAssay(seu) <- "gene"
  seu <- Seurat::CellCycleScoring(seu, s.features = Seurat::cc.genes.updated.2019$s.genes,
                                  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
  Seurat::DefaultAssay(seu) <- default_assay
  if (!resolution %in% colnames(seu@meta.data)) stop(gid, ": no ", resolution, " after clustering")
  seu <- assign_auto_phase_clusters(seu, resolution, seurat_assay = "gene")
  seu$s9_group <- factor(seu$s9_group, levels = c("without", "with"))
  seu$s9_integration <- gid
  seu@misc$s9_dropped_tumors <- dropped

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(gid, "_filtered_seu.rds"))
  saveRDS(seu, path)
  path
}

# Paired pseudobulk DESeq2 (~ tumor + group) on the cells of `seu`, tumors given
.s9_pseudobulk <- function(seu, tumors, min_count = 10) {
  seu <- seu[, which(seu$s9_tumor %in% tumors)]
  cnt <- SeuratObject::LayerData(seu, assay = "gene", layer = "counts")
  sample_of <- paste(seu$s9_tumor, seu$s9_group, sep = "|")
  mm <- Matrix::sparse.model.matrix(~ 0 + f, data.frame(f = factor(sample_of)))
  colnames(mm) <- sub("^f", "", colnames(mm))
  pb <- as.matrix(cnt %*% mm)
  pb <- round(pb[rowSums(pb) >= min_count, , drop = FALSE])
  cd <- data.frame(row.names = colnames(pb),
                   tumor = factor(sub("\\|.*", "", colnames(pb))),
                   group = factor(sub(".*\\|", "", colnames(pb)), levels = c("without", "with")))
  dds <- DESeq2::DESeqDataSetFromMatrix(pb, cd, design = ~ tumor + group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  DESeq2::results(dds, name = "group_with_vs_without") |>
    as.data.frame() |>
    tibble::rownames_to_column("symbol") |>
    tibble::as_tibble()
}

#' With-vs-without diffex inside each cluster of an S9 integration
#'
#' (a) per tumor: Wilcoxon in each (cluster, tumor) with at least `min_cells` per clone;
#'     gene locations from that tumor's own differing segments.
#' (b) pseudobulk: per cluster with at least `min_tumors` qualifying tumors, gene counts
#'     summed per tumor x clone, DESeq2 ~ tumor + group (Wald, with vs without). Also
#'     cluster "all" (every cell). n_tumors_cis = tumors in which the gene is cis.
#' Also writes the cell table (UMAP, clusters, group, tumor) for the plots.
#'
#' @param seu_path Output of integrate_s9_group()
#' @param group_tbl Rows of make_s9_integration_groups() for that group
#' @param min_cells,min_tumors Thresholds
#' @param out_dir Output directory
#' @return Paths: per_tumor, pseudobulk, counts and cells files
#' @export
s9_integrated_cluster_diffex <- function(seu_path, group_tbl, min_cells = 20, min_tumors = 2,
                                         out_dir = "results/fig_s9_t1e5/integrated") {
  gid <- unique(group_tbl$group_id)
  stopifnot(length(gid) == 1)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seu <- readRDS(seu_path)
  Seurat::DefaultAssay(seu) <- "gene"
  seu$s9_group <- as.character(seu$s9_group)
  group_tbl <- group_tbl[group_tbl$sample %in% unique(seu$s9_tumor), ]   # tumors left out by integrate_s9_group()
  md <- seu@meta.data

  umap <- Seurat::Embeddings(seu, "umap")
  cells_path <- file.path(out_dir, paste0(gid, "_cells.csv.gz"))
  tibble::tibble(cell = colnames(seu), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
                 clusters = as.character(md$clusters), phase_level = as.character(md$phase_level),
                 Phase = md$Phase, s9_group = md$s9_group, s9_tumor = md$s9_tumor, s9_key = md$s9_key,
                 group_id = gid) |>
    readr::write_csv(cells_path)

  loc <- purrr::map(rlang::set_names(seq_len(nrow(group_tbl)), group_tbl$sample), function(i) {
    r <- group_tbl[i, ]
    list(segs = s9_differing_segments(r$clone_post, r$joint_post, r$clone_with, r$clone_without),
         acq = stringr::str_split(r$acquired, ";")[[1]])
  })

  counts <- .s9_counts(md, c("clusters", "phase_level", "s9_tumor"), min_cells) |>
    dplyr::mutate(group_id = gid, error = NA_character_)

  # (a) per tumor
  per_tumor <- list()
  for (j in which(counts$tested)) {
    cl <- as.character(counts$clusters[j]); tu <- counts$s9_tumor[j]
    res <- tryCatch(.s9_wilcox(seu[, which(seu$clusters == cl & seu$s9_tumor == tu)]),
                    error = function(e) conditionMessage(e))
    if (is.character(res)) { counts$error[j] <- res; message(cl, " ", tu, ": ", res); next }
    per_tumor[[length(per_tumor) + 1]] <- res |>
      dplyr::left_join(s9_locate_genes(res$symbol, loc[[tu]]$segs, loc[[tu]]$acq), by = "symbol") |>
      dplyr::mutate(clusters = cl, phase_level = as.character(counts$phase_level[j]), s9_tumor = tu,
                    cells_with = counts$cells_with[j], cells_without = counts$cells_without[j])
  }
  per_tumor <- dplyr::bind_rows(per_tumor) |> dplyr::mutate(group_id = gid)

  # (b) pseudobulk
  cis_by_tumor <- purrr::imap_dfr(loc, function(l, tu) {
    s9_locate_genes(rownames(seu), l$segs, l$acq) |> dplyr::transmute(symbol, s9_tumor = tu, location)
  })
  n_cis <- cis_by_tumor |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(n_tumors_cis = sum(location == "cis"), .groups = "drop")
  pb_sets <- counts |>
    dplyr::filter(tested) |>
    dplyr::group_by(clusters, phase_level) |>
    dplyr::summarise(tumors = list(s9_tumor), .groups = "drop") |>
    dplyr::filter(lengths(tumors) >= min_tumors) |>
    dplyr::mutate(clusters = as.character(clusters), phase_level = as.character(phase_level))
  all_ok <- counts |>
    dplyr::group_by(s9_tumor) |>
    dplyr::summarise(ok = sum(cells_with) >= min_cells & sum(cells_without) >= min_cells) |>
    dplyr::filter(ok) |>
    dplyr::pull(s9_tumor)
  if (length(all_ok) >= min_tumors)
    pb_sets <- dplyr::bind_rows(tibble::tibble(clusters = "all", phase_level = "all", tumors = list(all_ok)), pb_sets)
  pseudobulk <- purrr::pmap_dfr(pb_sets, function(clusters, phase_level, tumors) {
    cells <- if (clusters == "all") colnames(seu) else colnames(seu)[seu$clusters == clusters]
    res <- tryCatch(.s9_pseudobulk(seu[, cells], tumors), error = function(e) conditionMessage(e))
    if (is.character(res)) { message("pseudobulk ", clusters, ": ", res); return(NULL) }
    res |> dplyr::mutate(clusters = clusters, phase_level = phase_level, n_tumors = length(tumors),
                         tumors = paste(tumors, collapse = ";"))
  })
  if (nrow(pseudobulk) > 0)
    pseudobulk <- pseudobulk |>
      dplyr::left_join(n_cis, by = "symbol") |>
      dplyr::left_join(dplyr::distinct(annotables::grch38, symbol, .keep_all = TRUE) |>
                         dplyr::select(symbol, chr), by = "symbol") |>
      dplyr::mutate(group_id = gid)

  paths <- c(per_tumor = file.path(out_dir, paste0(gid, "_per_tumor.csv")),
             pseudobulk = file.path(out_dir, paste0(gid, "_pseudobulk.csv")),
             counts = file.path(out_dir, paste0(gid, "_counts.csv")),
             cells = cells_path)
  readr::write_csv(per_tumor, paths[["per_tumor"]])
  readr::write_csv(pseudobulk, paths[["pseudobulk"]])
  readr::write_csv(counts, paths[["counts"]])
  unname(paths)
}

# ---- plots -------------------------------------------------------------------

.S9_P_ADJ_MAX <- 0.1
.S9_MIN_LOG2FC <- 0.5

# Genes significant in the most columns, then by mean |FC|; half up, half down if split
.s9_pick_genes <- function(d, n, split_sign, fc = "avg_log2FC", p = "p_val_adj") {
  sig <- d[!is.na(d[[p]]) & d[[p]] <= .S9_P_ADJ_MAX & abs(d[[fc]]) >= .S9_MIN_LOG2FC, ]
  if (nrow(sig) == 0) return(character())
  ranked <- sig |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(rec = dplyr::n_distinct(column), m = mean(.data[[fc]]), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(rec), dplyr::desc(abs(m)))
  if (!split_sign) return(utils::head(ranked$symbol, n))
  c(utils::head(ranked$symbol[ranked$m > 0], n / 2), utils::head(ranked$symbol[ranked$m < 0], n / 2))
}

.s9_bubble <- function(d, genes, title, fc = "avg_log2FC", p = "p_val_adj", fc_cap = 3, nlp_cap = 50) {
  if (length(genes) == 0)
    return(ggplot() + theme_void() + labs(title = title, subtitle = "no significant genes"))
  d <- d[d$symbol %in% genes, ]
  ord <- d |> dplyr::group_by(symbol) |> dplyr::summarise(m = mean(.data[[fc]])) |> dplyr::arrange(m) |> dplyr::pull(symbol)
  d |>
    dplyr::mutate(symbol = factor(symbol, levels = ord),
                  sig = factor(!is.na(.data[[p]]) & .data[[p]] <= .S9_P_ADJ_MAX & abs(.data[[fc]]) >= .S9_MIN_LOG2FC,
                               c(TRUE, FALSE), c("yes", "no")),
                  nlp = pmin(-log10(.data[[p]]), nlp_cap), lfc = .data[[fc]]) |>
    ggplot(aes(x = column, y = symbol)) +
    geom_point(aes(size = nlp, fill = lfc, alpha = sig), shape = 21, stroke = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "log2FC",
                         limits = c(-fc_cap, fc_cap), oob = scales::squish) +
    scale_alpha_manual(values = c(yes = 1, no = 0.3), name = "significant", drop = FALSE) +
    scale_size_area(max_size = 4, limits = c(0, nlp_cap), breaks = c(1, 5, 10, 25, 50),
                    name = "-log10 p adj\n(capped at 50)") +
    scale_x_discrete(drop = FALSE) +
    guides(alpha = guide_legend(override.aes = list(size = 3, fill = "grey50"))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid.major = element_line(linewidth = 0.2))
}

.s9_save_pages <- function(plots, sizes, out_prefix) {
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  pdfs <- character()
  for (nm in names(plots)) {
    ggsave(sprintf("%s_%s.png", out_prefix, nm), plots[[nm]], width = sizes[[nm]][1],
           height = sizes[[nm]][2], dpi = 150, bg = "white", limitsize = FALSE)
    pdfs[nm] <- tempfile(tmpdir = Sys.getenv("TMPDIR", tempdir()), fileext = ".pdf")
    ggsave(pdfs[nm], plots[[nm]], width = sizes[[nm]][1], height = sizes[[nm]][2],
           device = grDevices::cairo_pdf, limitsize = FALSE)
  }
  qpdf::pdf_combine(pdfs, paste0(out_prefix, ".pdf"))
  unlink(pdfs)
  c(paste0(out_prefix, ".pdf"), sprintf("%s_%s.png", out_prefix, names(plots)))
}

# Read and bind CSVs, skipping zero-row tables (they come back all-character)
.s9_read_many <- function(paths, ...) {
  purrr::map(paths, readr::read_csv, show_col_types = FALSE, ...) |>
    purrr::keep(~ nrow(.x) > 0) |>
    dplyr::bind_rows()
}

#' Bubble pages and tables for the per-tumor within-cluster S9 diffex
#'
#' One page per S9 page; columns = pair x cluster; cis and outside-cis panels.
#' @param diffex_paths Output of s9_within_cluster_diffex() (all branches)
#' @param pairs Output of read_s9_pairs()
#' @param out_prefix PDF/PNG prefix
#' @param xlsx Summary workbook
#' @return Paths written
#' @export
plot_s9_cluster_diffex <- function(diffex_paths, pairs,
                                   out_prefix = "doc/external_figures_and_tables/Fig_S9_cluster_diffex_t1e5",
                                   xlsx = "results/fig_s9_t1e5/cluster_diffex/fig_s9_cluster_diffex.xlsx",
                                   n_genes = 40) {
  diffex_paths <- unlist(diffex_paths)
  de <- .s9_read_many(stringr::str_subset(diffex_paths, "_counts\\.csv$", negate = TRUE),
                       col_types = readr::cols(chr = "c", edge_other = "c", cis_segments = "c", clusters = "c"))
  counts <- .s9_read_many(stringr::str_subset(diffex_paths, "_counts\\.csv$"),
                           col_types = readr::cols(clusters = "c", error = "c"))
  if (nrow(de) == 0) stop("no cluster was tested in any S9 pair")
  de <- de |> dplyr::mutate(column = sprintf("%s it%d %s\n%s", sample, iteration, pair, clusters))
  plots <- list(); sizes <- list()
  for (pg in unique(pairs$page)) {
    d <- dplyr::filter(de, page == pg)
    if (nrow(d) == 0) { message("page skipped (nothing tested): ", pg); next }
    lev <- d |> dplyr::distinct(column, sample, iteration, pair, phase_level, clusters) |>
      dplyr::arrange(sample, iteration, pair, phase_level, clusters) |> dplyr::pull(column)
    d$column <- factor(d$column, levels = unique(lev))
    gain <- !all(stringr::str_detect(stringr::str_split(d$acquired[1], ";")[[1]], "-$"))
    cis <- dplyr::filter(d, location == "cis"); trans <- dplyr::filter(d, location != "cis")
    g_cis <- .s9_pick_genes(cis, n_genes, split_sign = !gain)
    g_trans <- .s9_pick_genes(trans, n_genes, split_sign = TRUE)
    nm <- s9_slug(pg)
    plots[[nm]] <- (.s9_bubble(cis, g_cis, "cis (in the acquired RB SCNA segments)") |
                      .s9_bubble(trans, g_trans, "outside them")) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = paste("Fig. S9 within-cluster clone diffex:", pg),
        caption = stringr::str_wrap(paste0(
          "Clone with the SCNA vs its comparison clone inside each cluster (", d$resolution[1],
          ", phase-labelled) with >= 20 cells per clone; Wilcoxon, assay gene. Genes significant (p adj <= 0.1, ",
          "|log2FC| >= 0.5) in the most columns, then by mean |log2FC|; faded = not significant there."), 160),
        theme = theme(plot.caption = element_text(hjust = 0)))
    sizes[[nm]] <- c(max(9, 2 * (3 + 0.28 * nlevels(d$column))),
                     max(5, 2.8 + 0.12 * max(length(g_cis), length(g_trans))))
  }
  out <- .s9_save_pages(plots, sizes, out_prefix)
  sig <- de |> dplyr::filter(p_val_adj <= .S9_P_ADJ_MAX, abs(avg_log2FC) >= .S9_MIN_LOG2FC)
  summary <- counts |>
    dplyr::left_join(
      sig |> dplyr::mutate(where = ifelse(location == "cis", "cis", "trans"),
                           dir = ifelse(avg_log2FC > 0, "up", "down")) |>
        dplyr::count(key, clusters, where, dir) |>
        tidyr::pivot_wider(names_from = c(where, dir), values_from = n, names_prefix = "n_sig_", values_fill = 0),
      by = c("key", "clusters")) |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("n_sig_"), ~ tidyr::replace_na(.x, 0)))
  writexl::write_xlsx(list(
    readme = tibble::tibble(note = c(
      "Fig. S9 within-cluster clone diffex (numbatHelpers s9_within_cluster_diffex; targets s9_cluster_diffex).",
      "Each S9 pair (clone with the SCNA vs its comparison clone, t=1e-5 clone_post at the panel's iteration) tested inside each cluster of the filtered object.",
      "Clusters: SCT_snn_res.0.4, phase-labelled by assign_auto_phase_clusters(). Tested when both clones have >= 20 cells there.",
      "Wilcoxon (FindMarkers, assay gene, logfc.threshold 0, min.pct 0.05). Significant: p adj <= 0.1 and |log2FC| >= 0.5; positive = higher with the SCNA.",
      "Full tables: results/fig_s9_t1e5/cluster_diffex/<sample>_it<r>_<pair>.csv; cells per cluster: <key>_counts.csv.")),
    clusters = summary,
    all_significant = sig |> dplyr::select(-column) |> dplyr::arrange(page, key, clusters, location != "cis", p_val_adj)),
    xlsx)
  c(out, xlsx)
}

#' Pages and tables for the S9 clone-pair integrations
#'
#' One page per group: UMAPs (cluster, clone, tumor), per-tumor cluster composition by
#' clone, pseudobulk cis / outside-cis bubbles by cluster, and per-tumor recurrence.
#' @param diffex_paths Output of s9_integrated_cluster_diffex() (all branches)
#' @param groups Output of make_s9_integration_groups()
#' @param out_prefix PDF/PNG prefix
#' @param xlsx Summary workbook
#' @return Paths written
#' @export
plot_s9_integrations <- function(diffex_paths, groups,
                                 out_prefix = "doc/external_figures_and_tables/Fig_S9_integrations_t1e5",
                                 xlsx = "results/fig_s9_t1e5/integrated/fig_s9_integrations.xlsx",
                                 n_genes = 40) {
  diffex_paths <- unlist(diffex_paths)
  rd <- function(suffix, ...) .s9_read_many(stringr::str_subset(diffex_paths, paste0(suffix, "$")), ...)
  pt <- rd("_per_tumor\\.csv", col_types = readr::cols(chr = "c", clusters = "c"))
  pb <- rd("_pseudobulk\\.csv", col_types = readr::cols(chr = "c", clusters = "c"))
  ct <- rd("_counts\\.csv", col_types = readr::cols(clusters = "c", error = "c"))
  cells <- rd("_cells\\.csv\\.gz", col_types = readr::cols(clusters = "c"))
  if (nrow(pt) == 0) pt <- tibble::tibble(group_id = character(), symbol = character(), clusters = character(),
                                          location = character(), s9_tumor = character(),
                                          avg_log2FC = double(), p_val_adj = double())
  if (nrow(pb) == 0) pb <- tibble::tibble(group_id = character(), symbol = character(), clusters = character(),
                                          n_tumors = integer(), n_tumors_cis = integer(),
                                          log2FoldChange = double(), padj = double())
  plots <- list(); sizes <- list(); summ <- list()
  for (gid in unique(groups$group_id)) {
    g <- dplyr::filter(groups, group_id == gid)
    ce <- cells[cells$group_id == gid, ]
    if (nrow(ce) == 0) { message("group skipped (no cells table): ", gid); next }
    lev <- ce |> dplyr::distinct(phase_level, clusters) |> dplyr::arrange(phase_level, clusters) |> dplyr::pull(clusters)
    lev <- unique(lev)
    um <- function(col, title) {
      ggplot(ce[sample(nrow(ce)), ], aes(UMAP_1, UMAP_2, color = .data[[col]])) +
        geom_point(size = 0.2, alpha = 0.6) + labs(title = title, color = NULL) +
        theme_void(base_size = 8) + guides(color = guide_legend(override.aes = list(size = 2)))
    }
    comp <- ce |>
      dplyr::count(s9_tumor, s9_group, clusters) |>
      dplyr::group_by(s9_tumor, s9_group) |> dplyr::mutate(frac = n / sum(n)) |> dplyr::ungroup() |>
      dplyr::mutate(clusters = factor(clusters, levels = lev),
                    s9_group = factor(s9_group, levels = c("without", "with"))) |>
      ggplot(aes(x = s9_group, y = frac, fill = clusters)) + geom_col(width = 0.8) +
      facet_wrap(~s9_tumor, nrow = 1) + labs(title = "cluster composition by clone", x = NULL, y = "fraction of cells") +
      theme_minimal(base_size = 8)
    pbg <- dplyr::filter(pb, group_id == gid) |>
      dplyr::mutate(column = factor(sprintf("%s (n=%d)", clusters, n_tumors),
                                    levels = unique(sprintf("%s (n=%d)", clusters, n_tumors)[order(clusters != "all", match(clusters, lev))])))
    pb_cis <- dplyr::filter(pbg, n_tumors_cis > 0); pb_trans <- dplyr::filter(pbg, n_tumors_cis == 0 | is.na(n_tumors_cis))
    gain <- !all(stringr::str_detect(stringr::str_split(g$acquired[1], ";")[[1]], "-$"))
    gc_ <- .s9_pick_genes(pb_cis, n_genes, !gain, fc = "log2FoldChange", p = "padj")
    gt_ <- .s9_pick_genes(pb_trans, n_genes, TRUE, fc = "log2FoldChange", p = "padj")
    rec <- dplyr::filter(pt, group_id == gid, p_val_adj <= .S9_P_ADJ_MAX, abs(avg_log2FC) >= .S9_MIN_LOG2FC) |>
      dplyr::group_by(symbol, clusters, location) |>
      dplyr::summarise(n_tumors = dplyr::n_distinct(s9_tumor), m = mean(avg_log2FC), .groups = "drop")
    top_rec <- rec |> dplyr::group_by(symbol) |> dplyr::summarise(r = max(n_tumors), a = mean(abs(m))) |>
      dplyr::arrange(dplyr::desc(r), dplyr::desc(a)) |> utils::head(n_genes) |> dplyr::pull(symbol)
    rec_plot <- if (length(top_rec) == 0) ggplot() + theme_void() + labs(title = "per-tumor recurrence", subtitle = "no significant genes") else
      rec |> dplyr::filter(symbol %in% top_rec) |>
      dplyr::mutate(symbol = factor(symbol, levels = rev(top_rec)), clusters = factor(clusters, levels = lev)) |>
      ggplot(aes(clusters, symbol, size = n_tumors, fill = m, shape = location == "cis")) +
      geom_point(stroke = 0.3) + scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 22), name = "cis") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "mean log2FC", limits = c(-3, 3), oob = scales::squish) +
      scale_size_area(max_size = 4, name = "tumors significant") +
      labs(title = "per-tumor Wilcoxon: tumors significant per cluster", x = NULL, y = NULL) +
      theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    top <- um("clusters", "clusters") | um("s9_group", "clone") | um("s9_tumor", "tumor")
    bottom <- .s9_bubble(pb_cis, gc_, "pseudobulk: cis in >= 1 tumor", "log2FoldChange", "padj") |
      .s9_bubble(pb_trans, gt_, "pseudobulk: outside cis", "log2FoldChange", "padj") | rec_plot
    plots[[gid]] <- (top / comp / bottom) + plot_layout(heights = c(1, 0.6, 1.6)) +
      plot_annotation(
        title = paste("S9 clone-pair integration:", g$group_label[1]),
        subtitle = paste0(
          paste(sprintf("%s it%d %s", g$sample, g$iteration, g$pair)[g$sample %in% ce$s9_tumor], collapse = "; "),
          if (any(!g$sample %in% ce$s9_tumor))
            paste0("   (left out, < 20 cells in a clone: ", paste(g$sample[!g$sample %in% ce$s9_tumor], collapse = ", "), ")")),
        caption = stringr::str_wrap(paste0(
          "Each tumor contributes the two clones of its S9 pair (CCA integration, seuratTools). Clusters integrated_snn_res.0.4, ",
          "phase-labelled. Pseudobulk: per cluster with >= 2 tumors having >= 20 cells per clone, counts summed per tumor x clone, ",
          "DESeq2 ~ tumor + group; 'all' = every cell. Significant: p adj <= 0.1, |log2FC| >= 0.5; positive = higher with the SCNA."), 180),
        theme = theme(plot.caption = element_text(hjust = 0)))
    sizes[[gid]] <- c(16, 7 + 0.12 * max(length(gc_), length(gt_), length(top_rec)))
    summ[[gid]] <- pbg |>
      dplyr::filter(!is.na(padj), padj <= .S9_P_ADJ_MAX, abs(log2FoldChange) >= .S9_MIN_LOG2FC) |>
      dplyr::mutate(where = ifelse(n_tumors_cis > 0, "cis", "trans"), dir = ifelse(log2FoldChange > 0, "up", "down")) |>
      dplyr::count(group_id, clusters, n_tumors, where, dir) |>
      tidyr::pivot_wider(names_from = c(where, dir), values_from = n, names_prefix = "n_sig_pb_", values_fill = 0)
  }
  out <- .s9_save_pages(plots, sizes, out_prefix)
  writexl::write_xlsx(list(
    readme = tibble::tibble(note = c(
      "S9 clone-pair integrations (numbatHelpers integrate_s9_group / s9_integrated_cluster_diffex; targets s9_integrated_*).",
      "Groups: each S9 page with >= 2 tumors (page_*), and each RB SCNA acquired without another RB SCNA (alone_*). One pair per tumor.",
      "Clusters: integrated_snn_res.0.4, phase-labelled. per_tumor: Wilcoxon per (cluster, tumor) with >= 20 cells per clone.",
      "pseudobulk: DESeq2 ~ tumor + group per cluster with >= 2 such tumors; 'all' = all cells. n_tumors_cis = tumors where the gene is cis.",
      "Full tables: results/fig_s9_t1e5/integrated/<group_id>_{per_tumor,pseudobulk,counts}.csv; objects: output/seurat/integrated_s9/.")),
    groups = groups |> dplyr::select(group_id, group_label, group_kind, sample, iteration, pair, acquired, page, background, key),
    cells_per_cluster = ct,
    pseudobulk_sig_counts = dplyr::bind_rows(summ),
    pseudobulk_significant = pb |> dplyr::filter(!is.na(padj), padj <= .S9_P_ADJ_MAX, abs(log2FoldChange) >= .S9_MIN_LOG2FC) |>
      dplyr::arrange(group_id, clusters, padj)),
    xlsx)
  c(out, xlsx)
}
