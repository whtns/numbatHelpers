# Per-cell cell-type calls and cell-type exclusion (issue #54).
#
# A cell is excluded at the filtering stage only when BOTH
#   1. its per-cell call against the plae fetal-retina pseudobulk (30 types; "M2") is a
#      non-photoreceptor-lineage type, and
#   2. the matching marker panel is positive, calibrated per sample against depth-matched
#      confident tumor cells ("M5").
# M2 alone is not enough: on tumor cells its RBC / microglia calls are near-ties carried by
# housekeeping and ribosomal genes, with no marker expression (docs/tumor_celltype_label_audit.md).
# Ported from src/audit_composite_diploid.R (label_group, cor_call, log_pb, per_cell_call, M5)
# in external_rb_scrnaseq_proj.

#' Harmonise cell-type names onto shared groups
#'
#' Maps plae / Sridhar / clustifyr vocabularies onto: cone, pr_precursor, rod, bipolar, rbc,
#' microglia, other_immune, other_nonretinal, rpe, neurogenic, rpc, other_retinal.
#' @param x Character vector of cell-type names.
#' @return Character vector of group names (NA where x is NA).
#' @export
celltype_group <- function(x) {
  y <- tolower(x)
  out <- rep("other_retinal", length(y))
  out[grepl("^cone", y)] <- "cone"
  out[grepl("photoreceptor precursor", y)] <- "pr_precursor"
  out[grepl("^rods?$", y)] <- "rod"
  out[grepl("bipolar", y)] <- "bipolar"
  out[grepl("red blood", y)] <- "rbc"
  out[grepl("microglia", y)] <- "microglia"
  out[grepl("macrophage|monocyte|t-cell|t/nk|b-cell|natural killer|mast", y)] <- "other_immune"
  out[grepl("endothelial|vein|blood vessel|pericyte|fibroblast|schwann|melanocyte|oligodendro|astrocyte|muller", y)] <-
    "other_nonretinal"
  out[grepl("^rpe", y)] <- "rpe"
  out[grepl("neurogenic", y)] <- "neurogenic"
  out[grepl("rpc", y) & !grepl("rpe", y)] <- "rpc"
  out[is.na(x)] <- NA_character_
  out
}

#' Cell-type groups excluded at the filtering stage (decided 2026-11-13, #54)
#'
#' Every non-photoreceptor-lineage group with a marker panel. Exclusion also needs panel support
#' (see [celltype_marker_support()]).
#' @return Character vector of groups.
#' @export
celltype_m2_excluded_groups <- function() {
  c("rod", "bipolar", "rbc", "microglia", "other_immune", "other_nonretinal", "rpe")
}

#' Marker panels used to confirm an excluded call
#'
#' Panels from src/audit_composite_diploid.R (rod, rbc, microglia, bipolar) and
#' src/audit_diploid_identity_v2.R (glia, vascular, stromal, rpe), plus an immune panel that
#' covers lymphoid cells.
#' @return Named list of gene vectors.
#' @export
celltype_marker_panels <- function() {
  list(
    rod       = c("RHO", "NRL", "NR2E3", "GNAT1"),
    bipolar   = c("VSX2", "VSX1", "PRKCA", "CABP5", "GRM6", "TRPM1"),
    rbc       = c("HBB", "HBA1", "HBA2", "HBG1", "HBG2", "HBD"),
    microglia = c("PTPRC", "C1QA", "C1QB", "TYROBP", "CX3CR1", "P2RY12", "AIF1"),
    immune    = c("PTPRC", "CD74", "HLA-DRA", "CD3E", "NKG7", "CD68"),
    glia      = c("RLBP1", "GLUL", "SLC1A3", "GFAP", "CLU"),
    vascular  = c("PECAM1", "VWF", "CLDN5"),
    stromal   = c("COL1A2", "COL3A1", "DCN"),
    rpe       = c("RPE65", "BEST1", "PMEL"))
}

#' Which panels confirm each excluded group (any one suffices)
#' @return Named list: group -> panel names.
#' @export
celltype_group_panels <- function() {
  list(rod = "rod", bipolar = "bipolar", rbc = "rbc", microglia = "microglia",
       other_immune = c("immune", "microglia"),
       other_nonretinal = c("glia", "vascular", "stromal"),
       rpe = "rpe")
}

# log1p(CP10k) of an assay's counts (JoinLayers first for Assay5).
.lognorm_counts <- function(seu, assay = "gene") {
  a <- seu[[assay]]
  if (inherits(a, "Assay5")) a <- SeuratObject::JoinLayers(a)
  counts <- SeuratObject::LayerData(a, layer = "counts")
  list(m = log1p(Matrix::t(Matrix::t(counts) / pmax(Matrix::colSums(counts), 1)) * 1e4),
       umi = Matrix::colSums(counts))
}

#' Calibrated marker-panel support
#'
#' For each panel: score = mean log1p(CP10k) over its genes; cells are binned by log10 depth
#' (quantiles of the calibration cells); a cell is positive when its score exceeds the
#' calibration cells' q-quantile in its depth bin and is > 0.
#'
#' @param seu Seurat object.
#' @param calib Logical vector over cells: the calibration population (confident tumor cells).
#' @param q Quantile (default 0.99).
#' @param n_bins Depth bins (default 5).
#' @param assay Assay whose counts are used.
#' @return data.frame with `celltype_pos_<panel>` columns, rownames = cells. All FALSE when
#'   there are fewer than 20 calibration cells.
#' @export
celltype_marker_support <- function(seu, calib, q = 0.99, n_bins = 5L, assay = "gene", lognorm = NULL) {
  ln <- if (is.null(lognorm)) .lognorm_counts(seu, assay) else lognorm
  m <- ln$m; umi <- ln$umi
  panels <- celltype_marker_panels()
  out <- data.frame(row.names = colnames(m))
  calib <- !is.na(calib) & calib
  if (sum(calib) < 20) {
    message("celltype_marker_support: only ", sum(calib), " calibration cells; no cell is marker-positive")
    for (p in names(panels)) out[[paste0("celltype_pos_", p)]] <- FALSE
    return(out)
  }
  lu <- log10(umi + 1)
  brks <- unique(stats::quantile(lu[calib], seq(0, 1, length.out = n_bins + 1)))
  bin <- if (length(brks) < 2) rep(1L, length(lu)) else
    pmin(pmax(findInterval(lu, brks, all.inside = TRUE), 1L), length(brks) - 1L)
  for (p in names(panels)) {
    g <- intersect(panels[[p]], rownames(m))
    sc <- if (length(g)) Matrix::colMeans(m[g, , drop = FALSE]) else rep(0, ncol(m))
    thr <- tapply(sc[calib], bin[calib], stats::quantile, probs = q, na.rm = TRUE)
    t_cell <- thr[as.character(bin)]
    t_cell[is.na(t_cell)] <- max(thr, na.rm = TRUE)   # a bin with no calibration cells
    out[[paste0("celltype_pos_", p)]] <- sc > t_cell & sc > 0
  }
  out
}

#' Per-cell M2 cell-type call and cell-type exclusion
#'
#' Correlates each cell's log1p(CP10k) expression with the log1p(CP10k) plae pseudobulk (all 30
#' types) over the object's variable features and takes the best-correlated type. With `calib`,
#' also computes calibrated marker support and `celltype_exclude` = the call's group is in
#' [celltype_m2_excluded_groups()] AND one of its panels ([celltype_group_panels()]) is positive.
#'
#' @param seu Seurat object.
#' @param calib Optional logical vector over cells (confident tumor cells) for marker calibration.
#'   NULL: no marker columns and `celltype_exclude` is FALSE everywhere.
#' @param ref_path Path to the plae pseudobulk (genes x cell types, counts).
#' @param assay Assay whose counts are used (default "gene").
#' @param genes Optional gene set; default the object's variable features.
#' @param chunk Cells per correlation chunk.
#' @return `seu` with metadata `celltype_m2`, `celltype_m2_r`, `celltype_m2_margin`,
#'   `celltype_m2_group`, `celltype_pos_<panel>` and `celltype_exclude`.
#' @export
annotate_celltype_m2 <- function(seu, calib = NULL,
                                 ref_path = "/project2/cobrinik_1090/Homo_sapiens/numbat/plae_ref.rds",
                                 assay = "gene", genes = NULL, chunk = 2000L) {
  ln <- .lognorm_counts(seu, assay)
  calls <- celltype_m2_calls(seu, ref_path = ref_path, assay = assay, genes = genes, chunk = chunk, lognorm = ln)
  calls$celltype_exclude <- FALSE
  if (!is.null(calib)) {
    stopifnot(length(calib) == ncol(seu))
    pos <- celltype_marker_support(seu, calib, assay = assay, lognorm = ln)
    gp <- celltype_group_panels()
    pos <- pos[rownames(calls), , drop = FALSE]
    supported <- rep(FALSE, nrow(calls))
    for (g in celltype_m2_excluded_groups()) {
      rows <- which(calls$celltype_m2_group %in% g)
      if (!length(rows)) next
      cols <- paste0("celltype_pos_", gp[[g]])
      supported[rows] <- Reduce(`|`, lapply(cols, function(cc) pos[[cc]][rows]))
    }
    calls$celltype_exclude <- supported
    calls <- cbind(calls, pos)
  }
  Seurat::AddMetaData(seu, calls)
}

#' @rdname annotate_celltype_m2
#' @param lognorm Optional precomputed `.lognorm_counts()` result.
#' @return For `celltype_m2_calls()`, a data.frame with `celltype_m2`, `celltype_m2_r`,
#'   `celltype_m2_margin` and `celltype_m2_group`, rownames = cells.
#' @export
celltype_m2_calls <- function(seu,
                              ref_path = "/project2/cobrinik_1090/Homo_sapiens/numbat/plae_ref.rds",
                              assay = "gene", genes = NULL, chunk = 2000L, lognorm = NULL) {
  m <- (if (is.null(lognorm)) .lognorm_counts(seu, assay) else lognorm)$m

  ref <- readRDS(ref_path)
  ref <- as.matrix(ref[, setdiff(colnames(ref), "NA"), drop = FALSE])
  ref <- log1p(sweep(ref, 2, pmax(colSums(ref), 1), "/") * 1e4)

  if (is.null(genes)) {
    genes <- Seurat::VariableFeatures(seu)
    if (length(genes) < 500) {
      genes <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(
        seu, assay = assay, nfeatures = 2000, verbose = FALSE), assay = assay)
    }
  }
  genes <- intersect(genes, intersect(rownames(m), rownames(ref)))
  if (length(genes) < 100) stop("only ", length(genes), " genes shared with the plae reference")

  rg <- ref[genes, , drop = FALSE]
  idx <- split(seq_len(ncol(m)), ceiling(seq_len(ncol(m)) / chunk))
  res <- lapply(idx, function(j) {
    r <- suppressWarnings(cor(as.matrix(m[genes, j, drop = FALSE]), rg))
    r[!is.finite(r)] <- NA
    best <- apply(r, 1, function(v) if (all(is.na(v))) NA_integer_ else which.max(v))
    data.frame(
      celltype_m2 = colnames(rg)[best],
      celltype_m2_r = r[cbind(seq_len(nrow(r)), best)],
      celltype_m2_margin = apply(r, 1, function(v) {
        v <- sort(v, decreasing = TRUE); if (length(v) < 2) NA_real_ else v[1] - v[2]
      }),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  rownames(out) <- colnames(m)
  out$celltype_m2_group <- celltype_group(out$celltype_m2)
  out
}

#' Metadata columns written by [annotate_celltype_m2()] (refreshed on cached objects)
#' @return Character vector.
#' @export
celltype_meta_cols <- function() {
  c("celltype_m2", "celltype_m2_r", "celltype_m2_margin", "celltype_m2_group", "celltype_exclude",
    paste0("celltype_pos_", names(celltype_marker_panels())))
}

#' Fixed, high-contrast colours for the plae 30-type M2 calls
#'
#' Glasbey colours (pals::glasbey) with near-black entries dropped, assigned in lineage order so a
#' type has the same colour in every panel and every sample. "other (<0.2%)" and "no call" are
#' grey.
#' @return Named character vector of colours.
#' @export
celltype_m2_palette <- function() {
  types <- c(
    "Cones", "Photoreceptor Precursors", "Early RPCs", "RPCs", "Late RPCs", "Neurogenic Cells",
    "Horizontal Cells", "AC/HC_Precurs", "Amacrine Cells", "Retinal Ganglion Cells",
    "Rods", "Bipolar Cells", "Rod Bipolar Cells", "Red Blood Cells", "Microglia",
    "Macrophage", "Monocyte", "T-Cell", "B-Cell", "Natural Killer", "Mast",
    "Muller Glia", "Astrocytes", "Endothelial", "Pericytes", "Fibroblasts", "Vein",
    "Schwann", "Melanocytes", "RPE")
  g <- pals::glasbey(32)
  g <- g[!g %in% c("#000033", "#201A01")]
  c(stats::setNames(g[seq_along(types)], types),
    "other (<0.2%)" = "grey80", "no call" = "grey80")
}
