# Round-level QC for the rebuilt *_numbat.rds objects.
#
# WHY THIS EXISTS
#
# scripts/process_numbat_rds.R now builds each object at a per-sample consensus
# round chosen by src/select_numbat_round.R, instead of numbat's package default
# i = 2. Nothing in the RDS filename or its location records which round was
# used, so the only way to confirm what a given object actually holds is to open
# it and compare its segmentation against the segs_consensus_k.tsv files on disk.
#
# These functions do that from the RDS alone -- no Seurat object, no cluster
# dictionary, no clone simplifications. That independence is the point: if the
# heavy downstream targets look wrong, this tells you whether the numbat inputs
# are wrong or whether the problem is further down.
#
# The PDF has three pages: an SCNA map (per-chromosome, real Mb axis, arms,
# states, LLR and cell fractions, canonical RB windows shaded), the called-segment
# table behind it, and numbat's phylo heatmap for the cell-level view.
# numbat's own plot_consensus is deliberately NOT used -- it labels segments
# "a"/"b"/"c'" with no arm, no coordinates and no state name, which is not
# enough to inspect a call.

# Arm coverage required before a call counts as the canonical arm-level event.
#
# WHY THIS EXISTS. The original rule counted any del/loh distal to the chr16
# centromere as "16q_loss", and any amp in the 6p window as "6p_gain", with no
# size floor. That let focal specks masquerade as arm-level RB SCNAs:
# SRX10264524 was credited with 16q_loss on the strength of a 1.6 Mb deletion at
# 83.9-85.5 Mb -- 3% of the arm, 11 genes, LLR 8.7, invisible on the heatmap.
#
# WHERE 0.15 COMES FROM. Across the 124 non-zero (sample, event) coverages in
# the 39-sample SRX cohort at their selected rounds, the values fall into a
# cluster of 18 running 0.009-0.121, then a gap, then a dense continuum from
# 0.228 upward. 0.15 sits in that gap: it removes exactly the speck cluster and
# nothing else. Any value in (0.121, 0.170) gives the same answer.
#
# 0.25 was tried first and rejected. It was chosen from the 20 events gained by
# round selection, where the apparent gap runs 0.052-0.297 -- but that subset is
# biased, and against the full cohort 0.25 cuts the dense region, separating
# 0.233 from 0.257 for no reason. GISTIC's broad/focal convention of 0.50 is
# more aggressive still: it would discard SRX11133587's 23 Mb 16q deletion at
# LLR 9283, plainly an arm-level event at arm_frac 0.43.
RB_MIN_ARM_FRAC <- 0.15

#' Canonical RB SCNAs present in a numbat consensus segmentation
#'
#' An event counts only when the segments supporting it cover at least
#' `min_arm_frac` of the target arm. Without that floor a focal deletion of a
#' few Mb is indistinguishable from arm-level loss in the output, which is not
#' what "16q-" means in retinoblastoma. See [numbat_rb_arm_frac()].
#'
#' @param segs A `segs_consensus` table (data.frame/data.table) from a numbat
#'   object or a `segs_consensus_k.tsv` file.
#' @param min_arm_frac Minimum fraction of the target arm that supporting
#'   segments must cover. Defaults to `RB_MIN_ARM_FRAC` (0.15). Pass 0 to
#'   reproduce the old size-agnostic behaviour.
#' @return Character vector of event names, e.g. `c("1q_gain", "16q_loss")`,
#'   carrying an `arm_frac` attribute with the coverage of all five events.
#' @export
numbat_rb_events <- function(segs, min_arm_frac = RB_MIN_ARM_FRAC) {
  fr <- numbat_rb_arm_frac(segs)
  # fr > 0 as well as the threshold: an uncalled arm has coverage exactly 0, and
  # `0 >= 0` would otherwise credit every sample with all five events whenever
  # min_arm_frac is 0.
  hit <- names(fr)[fr > 0 & fr >= min_arm_frac]
  structure(hit, arm_frac = fr)
}

#' Which consensus round on disk a numbat object's segmentation came from
#'
#' Compares the object's `segs_consensus` against every `segs_consensus_k.tsv`
#' in the sample's numbat directory. Matching on the sorted segment breakpoints
#' rather than on the whole table is deliberate: the object's copy has been
#' through `relevel_chrom()` and column classes may differ, but breakpoints are
#' what define a round's segmentation.
#'
#' @param segs The object's `segs_consensus`.
#' @param sample_dir Directory holding the sample's per-round numbat output.
#' @param max_round Highest round index to look for.
#' @return Integer round index, or `NA_integer_` if nothing on disk matches.
#' @export
numbat_match_round <- function(segs, sample_dir, max_round = 12L) {
  if (is.null(segs) || !is.data.frame(segs) || nrow(segs) == 0) return(NA_integer_)
  if (!dir.exists(sample_dir)) return(NA_integer_)

  obs <- sort(as.numeric(segs$seg_start))
  hit <- NA_integer_
  for (k in seq_len(max_round)) {
    p <- file.path(sample_dir, sprintf("segs_consensus_%d.tsv", k))
    if (!file.exists(p)) next
    d <- tryCatch(data.table::fread(p, showProgress = FALSE), error = function(e) NULL)
    if (is.null(d) || nrow(d) != length(obs)) next
    if (isTRUE(all.equal(sort(as.numeric(d$seg_start)), obs))) hit <- k
  }
  hit
}

#' Inspect one rebuilt numbat RDS: summary row plus a diagnostic PDF
#'
#' Opens the RDS exactly once and reports what it actually contains -- which
#' consensus round, whether that is the round the selection manifest asked for,
#' which components survived the build, and which canonical RB events the
#' segmentation carries. Optionally writes a three-page PDF -- SCNA map, called
#' -segment table, phylo heatmap -- so the object can be eyeballed. A panel that
#' fails to render is replaced by a page naming the error rather than aborting.
#'
#' Nothing here writes into `output/numbat_sridhar/`; the RDS is read-only.
#'
#' @param numbat_rds_file Path to a `*_numbat.rds` file.
#' @param out_dir Directory for the diagnostic PDF.
#' @param manifest Round-selection manifest written by `src/select_numbat_round.R`.
#' @param make_pdf Whether to render the PDF. FALSE gives the table only.
#' @param width,height PDF dimensions in inches.
#' @return A one-row data.frame. `pdf` holds the PDF path (NA if not written).
#' @export
qc_numbat_rds <- function(numbat_rds_file,
                          out_dir  = "results/numbat_round_qc",
                          manifest = "results/numbat_selected_round.csv",
                          make_pdf = TRUE,
                          width    = 14,
                          height   = 8) {

  numbat_rds_file <- unname(numbat_rds_file[[1]])
  sample_id  <- stringr::str_extract(basename(numbat_rds_file), "SR[RX][0-9]+")
  sample_dir <- sub("_numbat\\.rds$", "", numbat_rds_file)
  finfo      <- file.info(numbat_rds_file)

  # Field access has to tolerate both the R6 Numbat object and the plain-list
  # "numbat_recovered" fallback that process_numbat_rds.R writes when
  # Numbat$new() cannot be constructed at all.
  fld <- function(nb, name) tryCatch(nb[[name]], error = function(e) NULL)

  nb <- tryCatch(readRDS(numbat_rds_file), error = function(e) NULL)

  na_row <- function(msg) data.frame(
    sample_id = sample_id, readable = FALSE, object_class = NA_character_,
    stamped_round = NA_integer_, manifest_round = NA_integer_,
    disk_round = NA_integer_, rebuilt = NA, round_ok = NA,
    n_cells = NA_integer_, n_clones = NA_integer_, n_segs_nonneutral = NA_integer_,
    n_rb_events = NA_integer_, rb_events = NA_character_,
    null_components = NA_character_, size_mb = round(finfo$size / 1e6, 1),
    mtime = as.character(finfo$mtime), pdf = NA_character_,
    note = msg, stringsAsFactors = FALSE
  )

  if (is.null(nb)) return(na_row("readRDS failed"))

  segs  <- fld(nb, "segs_consensus")
  cpost <- fld(nb, "clone_post")

  # Components the downstream pipeline reaches for. A NULL here is the usual
  # cause of a blank clone tree or an empty heatmap panel, so name them.
  wanted <- c("joint_post", "segs_consensus", "clone_post", "gtree",
              "bulk_clones", "mut_graph", "treeML")
  nulls  <- wanted[vapply(wanted, function(f) is.null(fld(nb, f)), logical(1))]

  state <- if (is.data.frame(segs)) {
    if ("cnv_state_post" %in% names(segs)) segs$cnv_state_post else segs$cnv_state
  } else NULL

  events     <- numbat_rb_events(segs)
  disk_round <- numbat_match_round(segs, sample_dir)

  stamped <- fld(nb, "selected_round")
  stamped <- if (is.null(stamped)) NA_integer_ else as.integer(stamped)

  man_round <- NA_integer_
  if (!is.null(manifest) && file.exists(manifest)) {
    # as.data.frame, not the data.table itself: inside `[.data.table` the `i`
    # expression is evaluated in the table's own frame, so a bare `sample_id`
    # on the right of `==` resolves to the COLUMN, not to this function's
    # argument, and every row matches.
    m <- tryCatch(as.data.frame(data.table::fread(manifest)), error = function(e) NULL)
    if (!is.null(m) && all(c("sample_id", "selected_round") %in% names(m))) {
      hit <- m[m$sample_id == sample_id, , drop = FALSE]
      if (nrow(hit) == 1L) man_round <- as.integer(hit$selected_round[1])
    }
  }

  clone_col <- if (is.data.frame(cpost) && "clone_opt" %in% names(cpost)) cpost$clone_opt else NULL

  # -------------------------------------------------------------------------
  # Diagnostic PDF: numbat's own two views, nothing added.
  # -------------------------------------------------------------------------
  pdf_path <- NA_character_
  pdf_note <- ""
  if (isTRUE(make_pdf)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path <- file.path(out_dir, paste0(sample_id, "_round_qc.pdf"))

    ttl <- sprintf("%s  |  round: stamped %s / manifest %s / on-disk %s  |  RB events: %s",
                   sample_id,
                   ifelse(is.na(stamped), "-", stamped),
                   ifelse(is.na(man_round), "-", man_round),
                   ifelse(is.na(disk_round), "no match", disk_round),
                   if (length(events)) paste(events, collapse = ", ") else "none")

    heat <- tryCatch({
      p <- nb$plot_phylo_heatmap(raster = FALSE, show_phylo = TRUE)
      p + patchwork::plot_annotation(title = ttl)
    }, error = function(e) { pdf_note <<- paste0("heatmap: ", conditionMessage(e)); NULL })

    # The SCNA pages. plot_consensus is deliberately not used: it labels segments
    # "a"/"b"/"c'" with no arm, no coordinates and no state, which is not enough
    # to inspect a call.
    scna_map <- tryCatch(
      plot_numbat_scna_map(nb, title = ttl),
      error = function(e) { pdf_note <<- paste0(pdf_note, " scna_map: ", conditionMessage(e)); NULL }
    )
    scna_tab <- tryCatch(
      plot_numbat_scna_table(numbat_scna_table(nb),
                             title = paste0(sample_id, " - called segments")),
      error = function(e) { pdf_note <<- paste0(pdf_note, " scna_table: ", conditionMessage(e)); NULL }
    )

    panels <- Filter(Negate(is.null), list(scna_map, scna_tab, heat))
    names(panels) <- c("scna_map", "scna_table", "phylo_heatmap")[
      c(!is.null(scna_map), !is.null(scna_tab), !is.null(heat))]

    if (length(panels) == 0) {
      pdf_path <- NA_character_
    } else {
      grDevices::pdf(pdf_path, width = width, height = height)
      on.exit(grDevices::dev.off(), add = TRUE)

      # ggplot is lazy: aesthetics are evaluated at PRINT time, not when the
      # object is built. A tryCatch around construction therefore catches
      # nothing, which is how a single degenerate sample (SRX10031191, whose
      # numbat run never finished, leaving a gtree with no `node` column) took
      # down a whole pipeline run. Guard the draw itself, and give the failed
      # panel a page saying why so the sample is not silently short a page.
      for (nm in names(panels)) {
        pn  <- panels[[nm]]
        err <- tryCatch({
          if (inherits(pn, "gtable") || inherits(pn, "grob")) {
            grid::grid.newpage(); grid::grid.draw(pn)
          } else {
            print(pn)
          }
          NULL
        }, error = function(e) conditionMessage(e))

        if (!is.null(err)) {
          # cli puts U+2139 and friends in its messages; the pdf device cannot
          # encode them and warns per glyph. Flatten to ASCII for both the page
          # and the CSV cell.
          err <- gsub("[^\x20-\x7E]", " ", gsub("[\r\n]+", " ", err))
          err <- trimws(gsub(" +", " ", err))
          pdf_note <- paste0(pdf_note, " ", nm, ": ", err, ";")
          grid::grid.newpage()
          grid::grid.text(
            paste0(sample_id, "\n\n", nm, " could not be rendered\n\n",
                   paste(strwrap(err, width = 70), collapse = "\n")),
            gp = grid::gpar(fontsize = 11, col = "firebrick"))
        }
      }
    }
  }

  data.frame(
    sample_id         = sample_id,
    readable          = TRUE,
    object_class      = paste(class(nb), collapse = ","),
    stamped_round     = stamped,
    manifest_round    = man_round,
    disk_round        = disk_round,
    # Only objects written by the round-selecting build carry selected_round, so
    # that stamp is what separates "rebuilt" from "left as it was". round_ok is
    # NA rather than FALSE for the latter: an object built before the change is
    # not wrong, it just predates the question.
    rebuilt           = !is.na(stamped),
    round_ok          = if (is.na(stamped)) NA else isTRUE(disk_round == man_round) && isTRUE(stamped == man_round),
    n_cells           = if (is.data.frame(cpost)) nrow(cpost) else NA_integer_,
    n_clones          = if (!is.null(clone_col)) length(unique(stats::na.omit(clone_col))) else NA_integer_,
    n_segs_nonneutral = if (!is.null(state)) sum(state != "neu", na.rm = TRUE) else NA_integer_,
    n_rb_events       = length(events),
    rb_events         = if (length(events)) paste(events, collapse = ";") else "",
    null_components   = if (length(nulls)) paste(nulls, collapse = ";") else "",
    size_mb           = round(finfo$size / 1e6, 1),
    mtime             = as.character(finfo$mtime),
    pdf               = pdf_path,
    note              = trimws(pdf_note),
    stringsAsFactors  = FALSE
  )
}
