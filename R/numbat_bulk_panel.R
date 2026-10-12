# Round-correct pseudobulk clone panels, with the allele/haplotype track.
#
# Why this file exists (github issues #41 and #42)
# ------------------------------------------------
# convert_numbat_pngs() builds the bulk-clone panel by converting
# bulk_clones_final.png out of the numbat output directory. That PNG is numbat's
# LAST consensus round. Since the majority-of-rounds change, <sample>_numbat.rds
# is rebuilt at a SELECTED round, which for 34 of 39 samples is not the last one
# -- and for 20 of them carries a different canonical RB event set. So the bulk
# panel and the heatmap beside it have been describing different rounds.
#
# The fix is to stop reading numbat's PNG and render from nb$bulk_clones, which
# is the selected round's own pseudobulk table (present in 39/39 objects). The
# same table carries the phased-haplotype columns -- pBAF, pAD, haplo_post,
# major_count/minor_count, theta_mle -- so numbat::plot_psbulk(allele_only=TRUE)
# gives us the haplotype view (#41) from exactly the same data, with no second
# source to keep in sync.
#
# The round is stamped into the panel title. A mismatch of this kind should
# never again be invisible on the figure itself.
#
# Consensus colouring (segments = "consensus", the default)
# ---------------------------------------------------------
# plot_psbulk() has no segment ribbon -- it colours the logFC and pHF points by
# `state_post`, which in a normal run follows `cnv_state_post`, the PER-CLONE
# RETEST. That retest skips every consensus-neutral segment outright and floors
# at min_LLR, so the panel shows a strict subset of the consensus segments and a
# different computation from the heatmap beside it: 43 of 252 clone-segment
# pairs silenced on SRX10264524 round 4.
#
# numbat_consensus_bulk() recolours the tracks by the consensus states instead,
# so the panel and the heatmap show the same segment set. The cost: the
# consensus call is asserted on every clone, including the root/normal clone,
# so colour marks "this segment is an event in this tumor", not "this clone
# carries it" -- per-clone carriage is read from the logFC offset and the pHF
# band separation. Pass segments = "retest" for numbat's own colouring.

#' Render a numbat object's pseudobulk clone profiles at its selected round
#'
#' Reads `bulk_clones` from a rebuilt `*_numbat.rds` and draws numbat's own
#' pseudobulk panel, optionally followed by an allele-only page showing the
#' phased-haplotype track. Both pages are titled with the consensus round the
#' object actually holds.
#'
#' Replaces `retrieve_numbat_plot_type(convert_numbat_pngs(...),
#' "bulk_clones_final.pdf")`, which returns numbat's final round regardless of
#' which round the object was built at.
#'
#' @param numbat_rds_file Path to a `*_numbat.rds` file.
#' @param out_dir Directory for the PDF.
#' @param manifest Round table (`sample_id` plus `selected_round` or `round`);
#'   used only to report the requested round alongside the one the object holds.
#'   Missing file is not an error.
#' @param allele_only_page Whether to append the allele/haplotype page.
#' @param segments Which colouring to use. `"clone_geno"` colours each clone
#'   only by the segments its phylogeny genotype (`GT_opt`) carries, so the
#'   panel agrees with the clone tree; `"consensus"` colours every clone by the
#'   consensus state; `"retest"` is numbat's own per-clone retest. See
#'   `R/numbat_bulk_colouring.R`.
#' @param depth_histograms Whether to draw a per-clone read-count histogram at
#'   the right of each clone's row. Needs round-correct clone membership; it is
#'   skipped (with a warning) when that is unavailable.
#' @param sqlite_path Metadata SQLite database holding `cell_qc_values`, read by
#'   [get_sample_cell_depth()] for the per-clone histograms.
#' @param width,height PDF dimensions in inches.
#' @return Path to the PDF, or `NA_character_` if nothing could be rendered.
#' @export
plot_numbat_bulk_clones <- function(numbat_rds_file,
                                    out_dir  = "results/numbat_bulk_clones",
                                    manifest = "results/numbat_selected_round.csv",
                                    allele_only_page = TRUE,
                                    segments = c("clone_geno", "consensus", "retest"),
                                    depth_histograms = TRUE,
                                    sqlite_path = "batch_hashes.sqlite",
                                    width  = 14,
                                    height = 8) {

  segments <- match.arg(segments)

  numbat_rds_file <- unname(unlist(numbat_rds_file)[[1]])
  sample_id  <- stringr::str_extract(basename(numbat_rds_file), "SR[RX][0-9]+")
  sample_dir <- sub("_numbat\\.rds$", "", numbat_rds_file)

  nb <- tryCatch(readRDS(numbat_rds_file), error = function(e) NULL)
  if (is.null(nb)) return(NA_character_)

  # Tolerate both the R6 Numbat object and the plain-list "numbat_recovered"
  # fallback that process_numbat_rds.R writes, as qc_numbat_rds() does.
  fld <- function(name) tryCatch(nb[[name]], error = function(e) NULL)

  # Which round the object holds.
  stamped <- fld("selected_round")
  held <- if (!is.null(stamped) && !is.na(stamped)) as.integer(stamped)
          else numbat_match_round(fld("segs_consensus"), sample_dir)

  # nb$bulk_clones is NOT the right table to colour per clone.
  # process_numbat_rds.R repoints it to bulk_clones_<held>, whose clone
  # MEMBERSHIP is round held-1's -- verified on SRX10264524, where
  # bulk_clones_4 carries 209/3850/2045/3169/938/1432 (= clone_post_3) while
  # nb$clone_post is 249/3804/1999/3226/906/1459 (= clone_post_4). Colouring
  # those rows with round held's GT_opt would attach each genotype to a
  # different clone. The table whose clones ARE round held's is
  # bulk_clones_<held+1>, or bulk_clones_final when held is the last round.
  bulk <- numbat_membership_bulk(sample_dir, held)
  membership_ok <- !is.null(bulk)
  if (!membership_ok) bulk <- fld("bulk_clones")
  if (!is.data.frame(bulk) || nrow(bulk) == 0) return(NA_character_)

  # Segment letters are round-specific, so nothing here joins on the letter
  # across rounds: numbat_clone_geno_bulk() resolves GT_opt to genomic
  # intervals from segs_consensus and matches them against the bulk table's own
  # segment extents.
  if (segments == "clone_geno" && membership_ok) {
    bulk <- numbat_clone_geno_bulk(bulk, fld("segs_consensus"), fld("clone_post"))
  } else if (segments == "consensus") {
    bulk <- numbat_consensus_bulk(bulk, fld("segs_consensus"))
  } else if (segments == "clone_geno") {
    # Falling back silently to numbat's retest colouring would look like a
    # per-clone genotype and would not be one.
    warning("no round-correct bulk_clones for ", sample_id,
            "; drawing numbat's retest colouring instead", call. = FALSE)
    segments <- "retest"
  }

  round_label <- numbat_round_label(nb, sample_dir, sample_id, manifest, fld)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(sample_id, "_bulk_clones.pdf"))

  events <- numbat_rb_events(fld("segs_consensus"))
  ttl <- sprintf("%s  |  %s  |  %s segments  |  RB events: %s",
                 sample_id, round_label, segments,
                 if (length(events)) paste(events, collapse = ", ") else "none")

  # plot_psbulk() resolves `gaps_hg38` as a BARE symbol with no argument to
  # override it. numbat ships that as LazyData, which lives in the package's
  # lazy-load DB and is reachable only once numbat is ATTACHED -- it is not in
  # the namespace, so importing is not enough and the panel dies with
  # "object 'gaps_hg38' not found". Attach once, idempotently.
  if (!"package:numbat" %in% search()) {
    suppressPackageStartupMessages(
      library(numbat, quietly = TRUE, warn.conflicts = FALSE))
  }

  # min_LLR = 0 only under consensus colouring: plot_psbulk() otherwise
  # re-applies its floor to the PER-CLONE LLR and neutralises most of what
  # numbat_consensus_bulk() just coloured. Under "retest" the floor is numbat's
  # own behaviour and is left alone.
  min_llr <- if (segments == "retest") 5 else 0

  # Per-clone depth. Only when the bulk rows ARE round held's clones: the
  # fallback nb$bulk_clones carries round held-1's membership, so joining it to
  # nb$clone_post would attach each clone's histogram to a different cell set.
  clone_depth <- NULL
  if (isTRUE(depth_histograms)) {
    if (membership_ok) {
      clone_depth <- numbat_clone_depth(fld("clone_post"),
        get_sample_cell_depth(sample_id, sqlite_path = sqlite_path))
    } else {
      warning("no round-correct bulk_clones for ", sample_id,
              "; per-clone depth histograms skipped", call. = FALSE)
    }
  }

  pages <- list(expression_and_allele =
                  function() plot_bulks_with_depth(bulk, clone_depth,
                                                   min_LLR = min_llr))
  if (isTRUE(allele_only_page)) {
    pages$allele_only <-
      function() plot_bulks_with_depth(bulk, clone_depth,
                                       min_LLR = min_llr, allele_only = TRUE)
  }

  grDevices::pdf(pdf_path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  drawn <- 0L
  for (nm in names(pages)) {
    # ggplot is lazy -- aesthetics evaluate at print time, so guarding
    # construction catches nothing. Guard the draw, and give a failed page a
    # page saying why rather than losing it silently.
    err <- tryCatch({
      p <- pages[[nm]]()
      sub <- if (nm == "allele_only") paste(ttl, " |  allele / haplotype") else ttl
      print(p + patchwork::plot_annotation(title = sub))
      drawn <- drawn + 1L
      NULL
    }, error = function(e) conditionMessage(e))

    if (!is.null(err)) {
      err <- trimws(gsub(" +", " ", gsub("[^\x20-\x7E]", " ",
                                         gsub("[\r\n]+", " ", err))))
      grid::grid.newpage()
      grid::grid.text(
        paste0(sample_id, "\n\n", nm, " could not be rendered\n\n",
               paste(strwrap(err, width = 70), collapse = "\n")),
        gp = grid::gpar(fontsize = 11, col = "firebrick"))
    }
  }

  if (drawn == 0L) return(NA_character_)
  pdf_path
}

#' Describe which consensus round a numbat object holds
#'
#' Prefers the round stamped onto the object at build time, falls back to
#' matching `segs_consensus` against the per-round tables on disk, and reports
#' the manifest's requested round alongside so a disagreement is visible.
#'
#' @param nb The numbat object.
#' @param sample_dir Directory holding the sample's per-round numbat output.
#' @param sample_id Sample accession.
#' @param manifest Path to the round-selection manifest, or NULL.
#' @param fld Field accessor tolerant of R6 and list objects.
#' @return A one-line character label.
#' @keywords internal
numbat_round_label <- function(nb, sample_dir, sample_id, manifest, fld) {
  stamped <- fld("selected_round")
  stamped <- if (is.null(stamped)) NA_integer_ else as.integer(stamped)

  on_disk <- numbat_match_round(fld("segs_consensus"), sample_dir)

  requested <- NA_integer_
  if (!is.null(manifest) && file.exists(manifest)) {
    m <- tryCatch(as.data.frame(data.table::fread(manifest)),
                  error = function(e) NULL)
    # numbat_selected_round*.csv calls it selected_round; numbat_active_round.csv
    # (the t=1e-5 cohort's table) calls it round.
    rcol <- intersect(c("selected_round", "round"), names(m))[1]
    if (!is.null(m) && "sample_id" %in% names(m) && !is.na(rcol)) {
      hit <- m[m$sample_id == sample_id, , drop = FALSE]
      if (nrow(hit) == 1L) requested <- as.integer(hit[[rcol]][1])
    }
  }

  held <- if (!is.na(stamped)) stamped else on_disk
  lab  <- paste0("round ", if (is.na(held)) "unknown" else held)

  # Only say anything about the manifest when it disagrees; a matching round is
  # the expected case and does not need column inches on every figure.
  if (!is.na(requested) && !is.na(held) && requested != held) {
    lab <- paste0(lab, " (manifest asked for ", requested, ")")
  }
  lab
}


#' The bulk_clones table whose clone membership is round k's
#'
#' `bulk_clones_<k>` holds round k's segments but round k-1's membership, so the
#' table describing round k's clones is `bulk_clones_<k+1>` -- or
#' `bulk_clones_final`, which is the only table whose membership and segments
#' are both the last round's.
#'
#' @param sample_dir The sample's numbat output directory.
#' @param k The round whose clones are wanted.
#' @return A data.frame, or `NULL` if no such table exists.
#' @keywords internal
numbat_membership_bulk <- function(sample_dir, k) {

  if (length(k) != 1L || is.na(k)) return(NULL)

  cands <- c(file.path(sample_dir, sprintf("bulk_clones_%d.tsv.gz", as.integer(k) + 1L)),
             file.path(sample_dir, "bulk_clones_final.tsv.gz"))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) return(NULL)

  tryCatch(as.data.frame(data.table::fread(hit[[1L]], showProgress = FALSE)),
           error = function(e) NULL)
}


#' Per-cell read counts keyed to numbat clone
#'
#' @param clone_post numbat's `clone_post` table (`cell`, `clone_opt`).
#' @param depth Output of [get_sample_cell_depth()].
#' @return A data frame with `clone` (character) and `nCount_gene`, or `NULL`
#'   when nothing joins.
#' @keywords internal
numbat_clone_depth <- function(clone_post, depth) {
  if (!is.data.frame(clone_post) || !all(c("cell", "clone_opt") %in% names(clone_post)) ||
      !is.data.frame(depth) || nrow(depth) == 0) return(NULL)
  # Barcodes are the same strings in both (AAACCTGAGAAGAAGC-1): numbat is run on
  # the unfiltered *_seu.rds, which is the stage get_sample_cell_depth() reads.
  out <- merge(clone_post[, c("cell", "clone_opt")],
               depth[, c("cell", "nCount_gene")], by = "cell")
  out <- out[!is.na(out$nCount_gene) & out$nCount_gene > 0, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  data.frame(clone = as.character(out$clone_opt), nCount_gene = out$nCount_gene,
             stringsAsFactors = FALSE)
}

#' numbat::plot_bulks() with a read-count histogram beside each clone's row
#'
#' Reproduces `numbat::plot_bulks(ncol = 1, title = TRUE)` row for row, then,
#' when `clone_depth` is given, puts that clone's per-cell `nCount_gene`
#' histogram at the right of its row. The x axis (log10) is shared by every
#' row, so clones can be compared by eye; the y axis is the fraction of that
#' clone's cells, so a 200-cell clone and a 4,000-cell clone draw at the same
#' height. Dashed red is the pipeline's `nCount_gene > 1000` filter; solid blue
#' the clone median; dotted grey the whole-sample median.
#'
#' @param bulk A `bulk_clones` table; `sample` is the clone id.
#' @param clone_depth Output of [numbat_clone_depth()], or `NULL` for plain
#'   numbat rows.
#' @param ... Passed to `numbat::plot_psbulk()` (e.g. `min_LLR`, `allele_only`).
#' @param threshold Count filter to mark.
#' @param title_size Row title size, as in `plot_bulks()`.
#' @param hist_width Relative width of the histogram column.
#' @return A patchwork.
#' @keywords internal
plot_bulks_with_depth <- function(bulk, clone_depth = NULL, ...,
                                  threshold = 1000, title_size = 8,
                                  hist_width = 0.12) {
  if (!"sample" %in% colnames(bulk)) bulk$sample <- 1

  rows <- split(bulk, bulk$sample)
  plot_list <- lapply(rows, function(b) {
    clone   <- unique(b$sample)
    n_cells <- unique(b$n_cells)
    p <- numbat::plot_psbulk(b, ...) +
      ggplot2::theme(title = ggplot2::element_text(size = title_size),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0, 0, 0.25, 0, "cm"))
    p + ggplot2::ggtitle(if (is.null(n_cells)) clone
                         else glue::glue("{clone} (n={n_cells})"))
  })

  if (is.null(clone_depth) || nrow(clone_depth) == 0) {
    return(patchwork::wrap_plots(plot_list, ncol = 1, guides = "collect"))
  }

  x_lim  <- range(clone_depth$nCount_gene)
  med_all <- stats::median(clone_depth$nCount_gene)
  last   <- names(rows)[length(rows)]

  hist_list <- lapply(names(rows), function(clone) {
    v <- clone_depth$nCount_gene[clone_depth$clone == clone]
    h <- ggplot2::ggplot() +
      ggplot2::scale_x_log10(limits = x_lim,
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      ggplot2::theme_classic(base_size = 7) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = title_size),
                     plot.margin = ggplot2::margin(0, 0.1, 0.25, 0.2, "cm"))
    if (length(v) == 0) {
      return(h + ggplot2::ggtitle("no depth data"))
    }
    h <- h +
      ggplot2::geom_histogram(data = data.frame(x = v),
                              ggplot2::aes(x = x, y = ggplot2::after_stat(count / sum(count))),
                              bins = 40, fill = "grey35", colour = NA) +
      ggplot2::geom_vline(xintercept = threshold, colour = "firebrick",
                          linetype = "dashed", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = med_all, colour = "grey50",
                          linetype = "dotted", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = stats::median(v), colour = "#2166AC",
                          linewidth = 0.5) +
      ggplot2::ggtitle(sprintf("median %s  <%s: %.0f%%",
                               format(round(stats::median(v)), big.mark = ","),
                               format(threshold, big.mark = ","),
                               100 * mean(v < threshold)))
    # Axis label once, on the bottom row, as the psbulk rows do with chromosomes.
    if (clone == last) h + ggplot2::labs(x = "nCount_gene")
    else h + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  })

  interleaved <- as.vector(rbind(plot_list, hist_list))
  patchwork::wrap_plots(interleaved, ncol = 2, widths = c(1, hist_width),
                        guides = "collect")
}
