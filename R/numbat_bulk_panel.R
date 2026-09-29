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
#' @param manifest Round-selection manifest; used only to report the requested
#'   round alongside the one the object holds. Missing file is not an error.
#' @param allele_only_page Whether to append the allele/haplotype page.
#' @param segments Which colouring to use. `"clone_geno"` colours each clone
#'   only by the segments its phylogeny genotype (`GT_opt`) carries, so the
#'   panel agrees with the clone tree; `"consensus"` colours every clone by the
#'   consensus state; `"retest"` is numbat's own per-clone retest. See
#'   `R/numbat_bulk_colouring.R`.
#' @param width,height PDF dimensions in inches.
#' @return Path to the PDF, or `NA_character_` if nothing could be rendered.
#' @export
plot_numbat_bulk_clones <- function(numbat_rds_file,
                                    out_dir  = "results/numbat_bulk_clones",
                                    manifest = "results/numbat_selected_round.csv",
                                    allele_only_page = TRUE,
                                    segments = c("clone_geno", "consensus", "retest"),
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

  pages <- list(expression_and_allele =
                  function() numbat::plot_bulks(bulk, ncol = 1, title = TRUE,
                                                min_LLR = min_llr))
  if (isTRUE(allele_only_page)) {
    pages$allele_only <-
      function() numbat::plot_bulks(bulk, ncol = 1, title = TRUE,
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
