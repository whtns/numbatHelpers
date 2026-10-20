#' Differential expression of one tumor clone against the aggregated diploid object
#'
#' For tumors with no usable within-sample contrast -- every tumor clone carries
#' the RB SCNA, or the diploid clone is too small -- the only available
#' comparator is the cohort-level diploid object built by
#' [assemble_diploid_seu()]. This runs that comparison with the same test and
#' annotation as [find_diffex_clones()] / `make_clone_comparison()`:
#' `FindMarkers(test.use = "wilcox", assay = "gene")`, then genes split into
#' `cis` (inside the clone's RB SCNA segments) and `trans` (everything else).
#'
#' Integration is deliberately NOT run. It is needed for the UMAP figures, not
#' for a Wilcoxon test on log-normalised counts, and it is the expensive step.
#' The two cell sets are merged at the counts level and normalised together.
#'
#' What this comparison is: the SCNA effect confounded with donor and batch.
#' The tumor cells come from one sample and the reference from others, so there
#' is no design in which those separate. Read it as tumor-versus-normal with the
#' named SCNAs among the differences.
#'
#' @param seu_path Path to the tumor's low-hypoxia Seurat object.
#' @param diploid_seu_path Path to the aggregated diploid object.
#' @param numbat_rds_file Path to the sample's numbat rds, for segment
#'   coordinates. If NULL or unreadable, every gene is reported as `trans`.
#' @param clone Clone number (`clone_opt`) standing in for the tumor.
#' @param rb_segs Segment ids of the clone's RB SCNAs, used for the cis split.
#'   Matched against `joint_post$seg`, falling back to `segs_consensus$seg_cons`.
#' @param exclude_own_cells Drop the tumor's own cells from the reference
#'   (leave-one-out), so every tumor is compared against other donors only.
#' @param min_cells Refuse to run when either side has fewer than this.
#' @param logfc_threshold Passed to [Seurat::FindMarkers()].
#' @param out_dir Directory for the result CSV.
#' @return Path to the CSV, or NULL when the comparison could not be run.
#' @export
find_diffex_clone_vs_diploid <- function(seu_path,
                                         diploid_seu_path,
                                         numbat_rds_file = NULL,
                                         clone,
                                         rb_segs = character(0),
                                         exclude_own_cells = TRUE,
                                         min_cells = 50L,
                                         logfc_threshold = 0.1,
                                         out_dir = "results/diploid_reference_diffex") {
  sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  old_assay_version <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")
  on.exit(options(Seurat.object.assay.version = old_assay_version), add = TRUE)

  # keep the raw gene counts only: the merged object is re-normalised below, and
  # SCT/integrated assays from the two sources are not comparable anyway
  prep_v3 <- function(s) {
    Seurat::DefaultAssay(s) <- "gene"
    if (inherits(s[["gene"]], "Assay5")) {
      counts_mat <- SeuratObject::GetAssayData(s[["gene"]], layer = "counts")
      s[["gene"]] <- SeuratObject::CreateAssayObject(counts = counts_mat)
    }
    for (a in setdiff(SeuratObject::Assays(s), "gene")) s[[a]] <- NULL
    s
  }

  tumor <- prep_v3(readRDS(seu_path))
  if (!"clone_opt" %in% colnames(tumor@meta.data))
    stop(sample_id, ": no clone_opt column; cannot select a clone")
  keep <- !is.na(tumor$clone_opt) & as.character(tumor$clone_opt) == as.character(clone)
  if (sum(keep) < min_cells) {
    warning(sample_id, ": clone ", clone, " has ", sum(keep), " cells (< ", min_cells, "); skipping.")
    return(NULL)
  }
  tumor <- tumor[, keep]

  ref <- prep_v3(readRDS(diploid_seu_path))
  if (exclude_own_cells && "sample_source" %in% colnames(ref@meta.data)) {
    own <- !is.na(ref$sample_source) & ref$sample_source == sample_id
    if (any(own)) ref <- ref[, !own]
  }
  # barcodes repeat across samples, so a shared name is not a shared cell; drop
  # them anyway rather than let merge() silently rename and pair them up
  shared <- intersect(colnames(tumor), colnames(ref))
  if (length(shared)) ref <- ref[, !colnames(ref) %in% shared]
  if (ncol(ref) < min_cells) {
    warning(sample_id, ": diploid reference has ", ncol(ref), " cells (< ", min_cells, "); skipping.")
    return(NULL)
  }

  tumor$comparison_group <- "tumor_clone"
  ref$comparison_group   <- "diploid_ref"
  n_tumor <- ncol(tumor); n_ref <- ncol(ref)

  seu <- merge(tumor, ref)
  if (inherits(seu[["gene"]], "Assay5")) seu[["gene"]] <- SeuratObject::JoinLayers(seu[["gene"]])
  seu <- Seurat::NormalizeData(seu, assay = "gene", verbose = FALSE)

  diffex <- Seurat::FindMarkers(seu, ident.1 = "tumor_clone", ident.2 = "diploid_ref",
                                group.by = "comparison_group", test.use = "wilcox",
                                assay = "gene", logfc.threshold = logfc_threshold) |>
    tibble::rownames_to_column("symbol") |>
    dplyr::left_join(annotables::grch38, by = "symbol") |>
    dplyr::distinct(ensgene, .keep_all = TRUE) |>
    # Canonical chromosomes only. make_clone_comparison()'s `!str_detect(chr,
    # "CHR_")` does not match annotables' alt/patch names (HSCHR3_1_CTG3,
    # HG405_PATCH, ...), which carry their own ensgene ids and so survive the
    # distinct() above -- MUC20-OT1 came back 8 times in the SRX10264520 test.
    dplyr::filter(!is.na(start), !is.na(end), chr %in% c(as.character(1:22), "X", "Y")) |>
    dplyr::distinct(symbol, .keep_all = TRUE)

  # cis = inside the clone's RB SCNA segments. Same source as
  # make_clone_comparison(): the segment coordinates numbat inferred for this
  # sample. Without them nothing can be called cis, so say so rather than
  # quietly reporting every gene as trans.
  segs <- .rb_segment_ranges(numbat_rds_file, clone, rb_segs)
  if (is.null(segs)) {
    warning(sample_id, ": no segment coordinates matched; every gene reported as trans.")
    diffex$location <- "trans"
  } else {
    gr <- diffex |> dplyr::mutate(seqnames = chr) |> plyranges::as_granges()
    cis_sym <- gr |> plyranges::join_overlap_intersect(segs) |> tibble::as_tibble() |>
      dplyr::pull(symbol) |> unique()
    diffex$location <- ifelse(diffex$symbol %in% cis_sym, "cis", "trans")
  }

  out <- diffex |>
    dplyr::mutate(sample_id = sample_id, clone = clone,
                  comparison = paste0("clone", clone, "_v_aggdiploid"),
                  n_cells_clone = n_tumor, n_cells_reference = n_ref,
                  reference_excludes_own_cells = exclude_own_cells) |>
    dplyr::relocate(sample_id, clone, comparison, location, symbol)

  out_path <- file.path(out_dir, paste0(sample_id, "_clone", clone, "_v_aggdiploid_diffex.csv"))
  readr::write_csv(out, out_path)
  message(sample_id, " clone ", clone, ": ", n_tumor, " tumor vs ", n_ref, " reference cells; ",
          nrow(out), " genes (", sum(out$location == "cis"), " cis)")
  out_path
}

# Segment coordinates for the clone's RB SCNAs, as a GRanges.
#
# A clone's GT_opt tokens are `seg_cons` ids, and segs_consensus carries their
# coordinates directly -- so look there first. joint_post$seg is a per-cell
# segment id covering only the part of the genome numbat evaluated per cell (7
# of 42 segments in SRX10264520), so matching there alone would drop most cis
# genes. It stays as a fallback for objects with no usable segs_consensus.
.rb_segment_ranges <- function(numbat_rds_file, clone, rb_segs) {
  if (is.null(numbat_rds_file) || !length(rb_segs) || !file.exists(numbat_rds_file)) return(NULL)
  nb <- tryCatch(readRDS(numbat_rds_file), error = function(e) NULL)
  if (is.null(nb)) return(NULL)

  sc <- nb$segs_consensus
  if (!is.null(sc) && all(c("seg_cons", "CHROM", "seg_start", "seg_end") %in% names(sc))) {
    hit <- sc[sc$seg_cons %in% rb_segs, , drop = FALSE]
    if ("cnv_state" %in% names(hit)) hit <- hit[!is.na(hit$cnv_state) & hit$cnv_state != "neu", , drop = FALSE]
    if (nrow(hit)) {
      return(hit |>
        dplyr::distinct(CHROM, seg_cons, seg_start, seg_end) |>
        dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) |>
        plyranges::as_granges())
    }
  }

  jp <- nb$joint_post
  if (is.null(jp)) return(NULL)
  if (!is.null(nb$clone_post)) {
    cells <- nb$clone_post$cell[as.character(nb$clone_post$clone_opt) == as.character(clone)]
    if (length(cells)) jp <- jp[jp$cell %in% cells, , drop = FALSE]
  }
  hit <- jp[jp$seg %in% rb_segs, , drop = FALSE]
  if (!nrow(hit)) return(NULL)

  hit |>
    dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) |>
    dplyr::filter(cnv_state_map != "neu") |>
    dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) |>
    plyranges::as_granges()
}
