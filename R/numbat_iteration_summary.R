# Per-iteration (consensus round) numbat summaries
#
# The pipeline collapses each numbat run to ONE consensus round: the RDS is built
# with round_mode = "final" (pipeline/scripts/process_numbat_rds.R) and every
# downstream product -- clone trees, heatmaps, clone diffex, the RB SCNA triage --
# reads only that round. results/numbat_active_round.csv records the result.
#
# That round was chosen by rule and, for most samples, never looked at. Because
# clone numbering is round-specific, a wrong round silently invalidates every
# clone key downstream -- which has already happened once (SRX11133593, see
# docs/rb_scna_triage.md). These functions lay the rounds side by side so the
# choice can be checked by eye instead of assumed.
#
# Layout: one COLUMN per round, three stacked panels per column
#   1. the fig_s03a-style numbat heatmap  (plot_numbat)
#   2. the waterfall                      (plot_variability_at_SCNA)
#   3. numbat's own bulk-clone panel      (bulk_clones_<k>.pdf, already on disk)
# plus one full-width header row carrying the smoothed-expression heatmap, which
# numbat writes ONCE per run rather than per round -- it is labelled as such so
# the reader does not take it for a per-round panel.
#
# Nothing here writes into a numbat output directory. SRX11133592/93/94 are held
# out by project policy (AGENTS.md) and their run dirs are read-only to us; the
# other samples' runs must not be silently mutated either. Note in particular
# that make_all_numbat_plots() regenerates PNGs back into the run dir -- it is
# deliberately NOT used.


#' Consensus rounds a numbat output directory can be loaded at
#'
#' A round is usable only when all nine artefacts `Numbat$new(i = k)` reads are
#' present. Round counts vary across the cohort -- four for most samples, two for
#' SRX11133592 and SRX11133594 -- so callers must discover them rather than
#' assume `1:4`.
#'
#' The same nine-file rule is also implemented in
#' `pipeline/scripts/process_numbat_rds.R` and `src/resolve_numbat_round.R`,
#' which run standalone under sbatch. This is the first packaged copy; the
#' others are deliberately left alone.
#'
#' @param sample_dir Directory holding one sample's per-round numbat output.
#' @param max_round Highest round index to look for.
#' @return Integer vector of complete rounds, ascending; `integer(0)` if none.
#' @export
numbat_complete_rounds <- function(sample_dir, max_round = 12L) {
  if (length(sample_dir) != 1L || is.na(sample_dir) || !dir.exists(sample_dir)) {
    return(integer(0))
  }
  artefacts <- function(k) sprintf(
    c("segs_consensus_%d.tsv", "clone_post_%d.tsv", "joint_post_%d.tsv",
      "exp_post_%d.tsv", "allele_post_%d.tsv", "geno_%d.tsv",
      "treeML_%d.rds", "mut_graph_%d.rds", "tree_final_%d.rds"), k)
  ks <- seq_len(as.integer(max_round))
  ok <- vapply(ks, function(k) {
    all(file.exists(file.path(sample_dir, artefacts(k))))
  }, logical(1L))
  ks[ok]
}


# --- private magick helpers -------------------------------------------------
#
# collate_sample_summary() has equivalents of these, but they are closures inside
# that function and so are not callable from here. Copies are deliberate: lifting
# them out would mean editing the function whose 33 outputs were just rebuilt and
# verified, for no benefit to this file. Two differences from the originals:
# these read PNG as well as PDF (so a raw numbat PNG works when no converted PDF
# exists), and every one tolerates NULL so a round that fails to render leaves a
# labelled gap instead of aborting the sample.

.iter_read_panel <- function(path, density = 300, page = 1L) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch({
    if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
      pages <- magick::image_read_pdf(path, density = density)
      if (page > length(pages)) return(NULL)
      pages[page]
    } else {
      magick::image_read(path)
    }
  }, error = function(e) NULL)
}

.iter_annotate <- function(img, label, size = 34L) {
  if (is.null(img)) return(NULL)
  info <- magick::image_info(img)[1L, ]
  # Size the strip by line count. A single long line overruns the column width
  # and is clipped mid-word -- which is how "ROUND 1" first rendered as "JND 1".
  n_lines <- length(strsplit(label, "\n", fixed = TRUE)[[1L]])
  strip_h <- max(as.integer(round(size * 1.7 * n_lines)) + 16L,
                 as.integer(round(info$height * 0.05)))
  strip <- magick::image_blank(info$width, strip_h, color = "white")
  strip <- magick::image_annotate(strip, label, size = size, color = "black",
                                  gravity = "center", weight = 700)
  magick::image_append(c(strip, img), stack = TRUE)
}

# Scale to the target width BEFORE annotating. Annotating first and scaling after
# shrinks the label text with the image, which at a 1200 px column renders a
# 34 pt strip at about 14 pt -- unreadable, and these labels carry the round
# number and RB event set that the whole panel is there to be judged by.
.iter_panel <- function(path, label, width = NULL, density = 300, page = 1L, size = 34L) {
  img <- .iter_read_panel(path, density = density, page = page)
  if (is.null(img)) return(NULL)
  if (!is.null(width)) img <- magick::image_scale(img, paste0(as.integer(width), "x"))
  .iter_annotate(img, label, size = size)
}

# A visible, labelled gap. Silence would read as "this round has no such panel",
# which is a different claim from "this panel could not be rendered".
.iter_missing <- function(width, height, label) {
  img <- magick::image_blank(max(200L, as.integer(width)), max(120L, as.integer(height)),
                             color = "grey94")
  magick::image_annotate(img, label, size = 34, color = "grey35",
                         gravity = "center", weight = 700)
}

.iter_stack <- function(img_list) {
  img_list <- purrr::compact(img_list)
  if (length(img_list) == 0) return(NULL)
  w <- max(vapply(img_list, function(i) magick::image_info(i)$width[1L], integer(1L)))
  padded <- lapply(img_list, function(img) {
    info <- magick::image_info(img)[1L, ]
    if (info$width < w) {
      magick::image_extent(img, glue::glue("{w}x{info$height}"),
                           gravity = "West", color = "white")
    } else img
  })
  do.call(c, padded) |> magick::image_append(stack = TRUE)
}

.iter_append <- function(cols) {
  cols <- purrr::compact(cols)
  if (length(cols) == 0) return(NULL)
  max_h <- max(vapply(cols, function(cc) magick::image_info(cc)$height[1L], integer(1L)))
  padded <- lapply(cols, function(col) {
    info <- magick::image_info(col)[1L, ]
    if (info$height < max_h) {
      magick::image_extent(col, glue::glue("{info$width}x{max_h}"),
                           gravity = "North", color = "white")
    } else col
  })
  do.call(c, padded) |> magick::image_append(stack = FALSE)
}

# Render a ggplot to a temporary PDF, returning NULL rather than erroring.
.iter_ggsave <- function(plot_obj, width, height) {
  if (is.null(plot_obj) || identical(plot_obj, NA_real_)) return(NULL)
  path <- tempfile(fileext = ".pdf")
  ok <- tryCatch({
    ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, limitsize = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(path)) NULL else path
}


#' Collate one sample's numbat consensus rounds side by side
#'
#' One column per complete consensus round, each stacking the fig_s03a-style
#' numbat heatmap, the waterfall (fig_s03a's second page -- P(SCNA) per cell,
#' faceted by segment) and numbat's per-round bulk-clone panel, with the
#' round-invariant smoothed-expression heatmap as a full-width header. The round
#' the sample's `*_numbat.rds` actually holds is marked `[ACTIVE]`, so what the
#' pipeline currently uses sits next to the alternatives it was chosen over.
#'
#' The heatmap needs neither a Seurat object nor the pipeline's clone assignment:
#' `plot_numbat()` never touches its `myseu` argument, and the only thing
#' `make_numbat_heatmaps()` takes from a Seurat object is `cell`/`clone_opt` --
#' which each round already publishes as `clone_post_<k>.tsv`.
#'
#' @param numbat_rds_file One sample's `*_numbat.rds` path. Only its name and
#'   location are used; the object itself is never read, since each round is
#'   loaded from the run directory instead.
#' @param numbat_plot_pdfs Converted numbat PDFs for this sample, i.e. one branch
#'   of `large_numbat_pdfs` (`convert_numbat_pngs()` output). Supplies the
#'   per-round `bulk_clones_<k>.pdf` and the `exp_roll_clust.pdf` header. When
#'   absent, the raw PNGs in the run directory are used instead.
#' @param eligible_samples_csv Triage table; samples whose `best_priority` is not
#'   in `eligible_priority` return `NULL`. Read from the table rather than
#'   hardcoded so the scope follows the triage when it is rerun. Pass `NULL` to
#'   collate every sample.
#' @param eligible_priority Which `best_priority` values count as eligible.
#' @param active_round_csv Table of the round each object holds, used only for the
#'   `[ACTIVE]` mark.
#' @param out_dir Where the PDF is written. Must not be under `output/`.
#' @param p_min,line_width Passed to `plot_numbat()`; the defaults match the
#'   `numbat_heatmap_plots_*` targets so the panels are comparable.
#' @param max_round Highest round index to look for.
#' @param col_width Pixel width each round's column is scaled to. This is the
#'   file-size control: composited at native 300-DPI size a four-round collage
#'   reaches 230 MB, which is too slow to open for the eyeballing this exists to
#'   support.
#' @param density Rasterisation DPI for reading panels.
#' @return Path to the written PDF, or `NULL` when the sample is not eligible or
#'   has no complete round.
#' @export
collate_iteration_summary <- function(numbat_rds_file,
                                      numbat_plot_pdfs = NULL,
                                      eligible_samples_csv = "results/diploid_audit/rb_scna_triage_samples.csv",
                                      eligible_priority = "P1_sufficient",
                                      active_round_csv = "results/numbat_active_round.csv",
                                      out_dir = "results/iteration_summaries",
                                      p_min = 0.9,
                                      line_width = 0.1,
                                      max_round = 12L,
                                      col_width = 1200L,
                                      density = 150) {

  rds <- unlist(numbat_rds_file, use.names = FALSE)
  rds <- rds[!is.na(rds)]
  if (length(rds) == 0) return(NULL)
  rds <- rds[[1L]]

  sample_id <- stringr::str_extract(rds, "SR[RX][0-9]+")
  if (is.na(sample_id)) return(NULL)

  # Hard stop rather than a warning: a path under output/ would write into a
  # numbat run directory, and three samples are held out by policy (AGENTS.md).
  if (grepl("(^|/)output/", out_dir)) {
    stop("collate_iteration_summary() must not write under output/: ", out_dir)
  }

  # --- eligibility ----------------------------------------------------------
  if (!is.null(eligible_samples_csv) && !is.na(eligible_samples_csv) &&
      nzchar(eligible_samples_csv) && file.exists(eligible_samples_csv)) {
    tri <- tryCatch(readr::read_csv(eligible_samples_csv, show_col_types = FALSE),
                    error = function(e) NULL)
    if (!is.null(tri) && all(c("sample_id", "best_priority") %in% names(tri))) {
      if (!sample_id %in% tri$sample_id[tri$best_priority %in% eligible_priority]) {
        return(NULL)
      }
    }
  }

  sample_dir <- sub("_numbat\\.rds$", "", rds)
  ks <- numbat_complete_rounds(sample_dir, max_round = max_round)
  if (length(ks) == 0) {
    warning("No complete consensus round for ", sample_id, " in ", sample_dir)
    return(NULL)
  }

  # --- which round the object actually holds --------------------------------
  active_k <- NA_integer_
  if (!is.null(active_round_csv) && !is.na(active_round_csv) &&
      nzchar(active_round_csv) && file.exists(active_round_csv)) {
    ar <- tryCatch(readr::read_csv(active_round_csv, show_col_types = FALSE),
                   error = function(e) NULL)
    if (!is.null(ar) && all(c("sample_id", "round") %in% names(ar))) {
      hit <- ar$round[ar$sample_id == sample_id]
      if (length(hit) > 0) active_k <- suppressWarnings(as.integer(hit[[1L]]))
    }
  }

  # --- this sample's converted panels, with a raw-PNG fallback --------------
  pdfs <- unlist(numbat_plot_pdfs, use.names = FALSE)
  pdfs <- pdfs[!is.na(pdfs) & nzchar(pdfs)]
  pdfs <- pdfs[stringr::str_detect(pdfs, sample_id)]

  pick_panel <- function(basename_pdf, basename_png = basename_pdf) {
    hit <- pdfs[basename(pdfs) == basename_pdf]
    if (length(hit) > 0 && file.exists(hit[[1L]])) return(hit[[1L]])
    fallback <- file.path(sample_dir, sub("\\.pdf$", ".png", basename_png))
    if (file.exists(fallback)) fallback else NA_character_
  }

  # --- one column per round -------------------------------------------------
  cols <- list()
  for (k in ks) {
    segs_path <- file.path(sample_dir, sprintf("segs_consensus_%d.tsv", k))
    cp_path   <- file.path(sample_dir, sprintf("clone_post_%d.tsv", k))

    segs <- tryCatch(as.data.frame(data.table::fread(segs_path, showProgress = FALSE)),
                     error = function(e) NULL)
    cp   <- tryCatch(as.data.frame(data.table::fread(cp_path, showProgress = FALSE)),
                     error = function(e) NULL)

    # numbat_rb_arm_frac() -- which numbat_rb_events() calls -- requires a
    # data.frame and returns all-zero coverage for anything else. Handing it the
    # TSV path would therefore report "no RB event" for every round, silently.
    ev <- tryCatch(as.character(numbat_rb_events(segs)), error = function(e) character(0))
    n_clones <- if (!is.null(cp) && "clone_opt" %in% names(cp)) {
      length(unique(stats::na.omit(cp$clone_opt)))
    } else NA_integer_

    lbl <- sprintf(
      "ROUND %d%s   |   %s clones\n%s",
      k,
      if (!is.na(active_k) && k == active_k) "   [ACTIVE]" else "",
      if (is.na(n_clones)) "?" else as.character(n_clones),
      if (length(ev) > 0) paste(ev, collapse = ", ") else "no canonical RB event"
    )

    # Build the round's numbat object once and use it for both rendered panels.
    # Numbat$new() reads bulk_clones_final.tsv.gz whatever i is, so the bulk
    # panel below deliberately comes from bulk_clones_<k> on disk instead.
    # gtf must be passed explicitly. Numbat$new()'s default for it is the bare
    # symbol `gtf_hg38`, which resolves only when numbat is attached with
    # library(); called namespaced -- as it is here and on every crew worker --
    # the default fails with "object 'gtf_hg38' not found" and every rendered
    # panel silently degrades to a grey "not rendered" box.
    nb_k <- tryCatch(
      numbat::Numbat$new(out_dir = sample_dir, i = k, gtf = numbat::gtf_hg38),
      error = function(e) {
        warning("Numbat$new() failed for ", sample_id, " round ", k, ": ",
                conditionMessage(e), call. = FALSE)
        NULL
      })

    hm_path <- NULL
    wf_path <- NULL
    if (!is.null(nb_k)) {
      clone_annot <- if (!is.null(cp) && all(c("cell", "clone_opt") %in% names(cp))) {
        ca <- cp[, c("cell", "clone_opt")]
        ca$clone_opt <- as.character(ca$clone_opt)
        ca
      } else NULL

      if (!is.null(clone_annot)) {
        hm <- safe_plot_numbat(
          nb_k,
          myseu = NULL,                 # unused by plot_numbat(); see file header
          myannot = clone_annot,
          mytitle = glue::glue("{sample_id} - round {k}"),
          clone_bar = FALSE,
          p_min = p_min,
          line_width = line_width,
          show_segment_names_on_x = TRUE
        )
        if (!is.null(hm[["error"]])) {
          warning("plot_numbat() failed for ", sample_id, " round ", k, ": ",
                  conditionMessage(hm[["error"]]), call. = FALSE)
        }
        # Taller than make_numbat_heatmaps' 10x5: this is the panel the round
        # judgement actually rests on, and at 2:1 it is the shortest row on the
        # page, dwarfed by the SCNA map beneath it.
        hm_res  <- hm[["result"]]
        hm_path <- .iter_ggsave(hm_res, width = 10, height = 7)

        # The waterfall is fig_s03a's SECOND page: p_cnv for every cell, faceted
        # by segment, with a clone tile beneath. It is derived from panel 3 of the
        # heatmap just rendered -- the same recipe make_numbat_heatmaps() uses for
        # its *_scna_var.pdf -- rather than drawn independently, so it is
        # guaranteed to describe the same cells and segments as the heatmap it
        # sits under. That also means it cannot exist when the heatmap fails,
        # which is why both live in this one branch.
        if (!is.null(hm_res) && !identical(hm_res, NA_real_) && length(hm_res) >= 3) {
          wf <- tryCatch(
            plot_variability_at_SCNA(
              dplyr::left_join(hm_res[[3]][["data"]], clone_annot, by = "cell"),
              p_min = p_min),
            error = function(e) {
              warning("plot_variability_at_SCNA() failed for ", sample_id, " round ",
                      k, ": ", conditionMessage(e), call. = FALSE)
              NULL
            })
          wf_path <- .iter_ggsave(wf, width = 12, height = 9)
        }
      }
    }

    # Free before the next round: joint_post alone runs past 100 MB per round, so
    # holding four of them is how this OOMs.
    rm(nb_k, segs, cp)
    gc(verbose = FALSE)

    bulk_path <- pick_panel(sprintf("bulk_clones_%d.pdf", k))

    hm_img   <- .iter_panel(hm_path,   lbl,            width = col_width, density = density, size = 30L)
    wf_img   <- .iter_panel(wf_path,   "Waterfall: P(SCNA) per cell, by segment",
                            width = col_width, density = density, size = 26L)
    bulk_img <- .iter_panel(bulk_path, "Bulk clones",  width = col_width, density = density, size = 26L)

    if (is.null(hm_img))   hm_img   <- .iter_missing(col_width, 500L, paste0(lbl, " -- heatmap not rendered"))
    if (is.null(wf_img))   wf_img   <- .iter_missing(col_width, 500L, "waterfall not rendered")
    if (is.null(bulk_img)) bulk_img <- .iter_missing(col_width, 500L, sprintf("bulk_clones_%d not found", k))

    cols[[length(cols) + 1L]] <- .iter_stack(list(hm_img, wf_img, bulk_img))
  }

  body <- .iter_append(cols)
  if (is.null(body)) {
    warning("No round rendered for ", sample_id)
    return(NULL)
  }

  rows <- list()

  # Smoothed expression: numbat computes this once, before the consensus loop, so
  # it is genuinely round-invariant. Shown once and labelled, rather than repeated
  # in every column where it would read as a per-round panel.
  exp_img <- .iter_panel(
    pick_panel("exp_roll_clust.pdf"),
    "Smoothed expression -- numbat writes this ONCE per run, not per round",
    width = magick::image_info(body)$width[1L],
    density = density,
    size = 40L
  )
  if (!is.null(exp_img)) rows[[length(rows) + 1L]] <- exp_img
  rows[[length(rows) + 1L]] <- body

  max_width <- max(vapply(rows, function(r) magick::image_info(r)$width[1L], integer(1L)))
  rows <- lapply(rows, function(r) {
    info <- magick::image_info(r)[1L, ]
    if (info$width < max_width) {
      magick::image_extent(r, glue::glue("{max_width}x{info$height}"),
                           gravity = "West", color = "white")
    } else r
  })
  out_img <- do.call(c, rows) |> magick::image_append(stack = TRUE)

  title <- sprintf("%s  --  numbat consensus rounds %s%s",
                   sample_id, paste(ks, collapse = ", "),
                   if (is.na(active_k)) "" else sprintf("  (object holds round %d)", active_k))
  info <- magick::image_info(out_img)[1L, ]
  title_h <- max(80L, as.integer(round(info$height * 0.02)))
  title_strip <- magick::image_blank(info$width, title_h, color = "white")
  title_strip <- magick::image_annotate(title_strip, title, size = 58, color = "black",
                                        gravity = "north", weight = 700, location = "+0+14")
  out_img <- magick::image_append(c(title_strip, out_img), stack = TRUE)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir, sprintf("%s_iterations.pdf", sample_id))
  magick::image_write(out_img, format = "pdf", path = out_path)
  out_path
}
