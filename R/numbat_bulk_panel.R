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
#' @param manifest Round-selection manifest; used only to report the requested
#'   round alongside the one the object holds. Missing file is not an error.
#' @param allele_only_page Whether to append the allele/haplotype page.
#' @param width,height PDF dimensions in inches.
#' @return Path to the PDF, or `NA_character_` if nothing could be rendered.
#' @export
plot_numbat_bulk_clones <- function(numbat_rds_file,
                                    out_dir  = "results/numbat_bulk_clones",
                                    manifest = "results/numbat_selected_round.csv",
                                    allele_only_page = TRUE,
                                    width  = 14,
                                    height = 8) {

  numbat_rds_file <- unname(unlist(numbat_rds_file)[[1]])
  sample_id  <- stringr::str_extract(basename(numbat_rds_file), "SR[RX][0-9]+")
  sample_dir <- sub("_numbat\\.rds$", "", numbat_rds_file)

  nb <- tryCatch(readRDS(numbat_rds_file), error = function(e) NULL)
  if (is.null(nb)) return(NA_character_)

  # Tolerate both the R6 Numbat object and the plain-list "numbat_recovered"
  # fallback that process_numbat_rds.R writes, as qc_numbat_rds() does.
  fld <- function(name) tryCatch(nb[[name]], error = function(e) NULL)

  bulk <- fld("bulk_clones")
  if (!is.data.frame(bulk) || nrow(bulk) == 0) return(NA_character_)

  round_label <- numbat_round_label(nb, sample_dir, sample_id, manifest, fld)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(sample_id, "_bulk_clones.pdf"))

  events <- numbat_rb_events(fld("segs_consensus"))
  ttl <- sprintf("%s  |  %s  |  RB events: %s",
                 sample_id, round_label,
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

  pages <- list(expression_and_allele =
                  function() numbat::plot_bulks(bulk, ncol = 1, title = TRUE))
  if (isTRUE(allele_only_page)) {
    pages$allele_only <-
      function() numbat::plot_bulks(bulk, ncol = 1, title = TRUE,
                                    allele_only = TRUE)
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
    if (!is.null(m) && all(c("sample_id", "selected_round") %in% names(m))) {
      hit <- m[m$sample_id == sample_id, , drop = FALSE]
      if (nrow(hit) == 1L) requested <- as.integer(hit$selected_round[1])
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
