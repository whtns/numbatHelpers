# How the pseudobulk clone panel gets its colours.
#
# plot_psbulk() has no segment ribbon: it draws logFC per gene and pHF per SNP,
# point-coloured by `state_post`. In a normal run `state_post` follows
# `cnv_state_post`, the PER-CLONE RETEST -- each clone's own pseudobulk tested
# against every consensus segment, floored at min_LLR. Three colourings are
# therefore available, and they are genuinely different claims:
#
#   "retest"      numbat's own. Colours where that clone's pseudobulk passes
#                 LLR >= min_LLR. Silences 43 of 252 clone-segment pairs on
#                 SRX10264524 round 4, and is not what the clone tree shows.
#
#   "consensus"   colours every clone by the consensus state. Makes the panel
#                 agree with the heatmap, but asserts the call on every clone
#                 including the root/normal one, so colour stops meaning
#                 "this clone carries it".
#
#   "clone_geno"  colours each clone only by the segments the PHYLOGENY says it
#                 carries (`GT_opt` in clone_post). This is what the clone tree
#                 panel draws, so the two agree by construction. Trunk events
#                 inherited by a descendant are coloured even where that clone's
#                 own pseudobulk is too weak to call them; segments `geno`
#                 dropped as not tree-informative are never coloured.
#
# See docs/numbat_per_clone_retest.md for the measured differences.
#
# ROUND ALIGNMENT -- the trap in all of this
# ------------------------------------------
# `bulk_clones_<k>` holds round k's SEGMENTS but round k-1's clone MEMBERSHIP
# (numbat pseudobulks the inherited clones at the start of iteration k, before
# the phylogeny that defines round k's). `bulk_clones_final` is the exception:
# its membership AND its segments are both the last round's.
#
# So the bulk table whose clones are round k's is `bulk_clones_<k+1>`, or
# `bulk_clones_final` when k is the last round -- and its `seg_cons` letters
# then belong to round k+1, while `GT_opt` names round k's. Segment letters are
# round-specific, so these functions never join on the letter across rounds:
# carried segments are resolved to GENOMIC INTERVALS from the caller's
# `segs_consensus` and matched against the bulk table's own segment extents.

#' Per-clone carried segments, as genomic intervals
#'
#' @param segs A `segs_consensus` table for the round `clone_post` came from.
#' @param clone_post A `clone_post` table carrying `clone_opt` and `GT_opt`.
#' @return data.frame of clone / CHROM / gs / ge / cons_state, or NULL.
#' @keywords internal
numbat_clone_intervals <- function(segs, clone_post) {

  if (!is.data.frame(segs) || !is.data.frame(clone_post)) return(NULL)
  if (!all(c("clone_opt", "GT_opt") %in% names(clone_post))) return(NULL)
  if (!all(c("seg_cons", "CHROM", "seg_start", "seg_end",
             "cnv_state", "cnv_state_post") %in% names(segs))) return(NULL)

  gt <- clone_post |>
    dplyr::distinct(clone_opt, GT_opt) |>
    dplyr::filter(!is.na(clone_opt))
  if (nrow(gt) == 0) return(NULL)

  carried <- do.call(rbind, lapply(seq_len(nrow(gt)), function(i) {
    g <- gt$GT_opt[i]
    # The root clone's GT_opt is the empty string: it carries nothing, and that
    # is a real answer, not missing data.
    if (is.na(g) || !nzchar(trimws(g))) return(NULL)
    data.frame(clone    = as.character(gt$clone_opt[i]),
               seg_cons = trimws(strsplit(g, ",")[[1]]),
               stringsAsFactors = FALSE)
  }))
  if (is.null(carried) || nrow(carried) == 0) return(NULL)

  seg_iv <- segs |>
    dplyr::distinct(seg_cons, CHROM, seg_start, seg_end, cnv_state, cnv_state_post) |>
    dplyr::transmute(
      seg_cons,
      CHROM = as.character(CHROM),
      gs    = seg_start,
      ge    = seg_end,
      # The state numbat itself uses for the per-cell posteriors (main.R:442).
      cons_state = ifelse(cnv_state == "neu", cnv_state, cnv_state_post))

  out <- dplyr::inner_join(carried, seg_iv, by = "seg_cons")
  if (nrow(out) == 0) return(NULL)
  out[, c("clone", "CHROM", "gs", "ge", "cons_state")]
}

#' Colour a pseudobulk table by each clone's phylogeny genotype
#'
#' Overwrites `cnv_state` / `cnv_state_post` / `state_post` so that a segment is
#' coloured in a clone only when that clone's `GT_opt` carries it, using the
#' consensus state for the colour. Everything else becomes neutral.
#'
#' Pass the result to `plot_bulks(..., min_LLR = 0)`; `plot_psbulk()` would
#' otherwise re-apply its floor to the per-clone `LLR` and neutralise most of
#' what this just coloured.
#'
#' @param bulk A `bulk_clones` table whose `sample` column is the clone id.
#'   Its membership must match `clone_post` -- see the round-alignment note at
#'   the top of this file.
#' @param segs The `segs_consensus` table for `clone_post`'s round.
#' @param clone_post The matching `clone_post`, carrying `clone_opt`/`GT_opt`.
#' @param min_overlap Fraction of a bulk segment that must fall inside a carried
#'   interval for it to be coloured.
#' @return `bulk`, recoloured. Unchanged if the inputs are unusable.
#' @export
numbat_clone_geno_bulk <- function(bulk, segs, clone_post, min_overlap = 0.5) {

  if (!is.data.frame(bulk)) return(bulk)
  if (!all(c("seg_cons", "CHROM", "POS", "sample") %in% names(bulk))) return(bulk)

  iv <- numbat_clone_intervals(segs, clone_post)
  if (is.null(iv)) return(bulk)

  # The bulk table's OWN segment extents, in its own segment space. Derived from
  # POS rather than trusting seg_start/seg_end, which after annot_consensus()
  # may describe the HMM segment rather than the consensus one.
  bulk_iv <- bulk |>
    dplyr::group_by(CHROM, seg_cons) |>
    dplyr::summarise(bs = min(POS, na.rm = TRUE), be = max(POS, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(CHROM = as.character(CHROM))

  clones <- unique(iv$clone)
  grid <- merge(data.frame(clone = clones, stringsAsFactors = FALSE),
                as.data.frame(bulk_iv), by = NULL)

  hits <- grid |>
    dplyr::inner_join(iv, by = c("clone", "CHROM")) |>
    dplyr::mutate(ov = pmax(0, pmin(be, ge) - pmax(bs, gs))) |>
    dplyr::filter(ov > 0, ov >= min_overlap * pmax(be - bs, 1)) |>
    # A bulk segment spanning two carried intervals takes the one it overlaps
    # most, rather than an arbitrary first match.
    dplyr::group_by(clone, CHROM, seg_cons) |>
    dplyr::slice_max(ov, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(clone, CHROM, seg_cons, geno_state = cons_state)

  bulk |>
    dplyr::mutate(.clone = as.character(sample), .chrom = as.character(CHROM)) |>
    dplyr::left_join(hits, by = c(".clone" = "clone", ".chrom" = "CHROM", "seg_cons")) |>
    dplyr::mutate(
      geno_state = ifelse(is.na(geno_state), "neu", geno_state),
      # p_up is produced by classify_alleles(), which runs only on segments the
      # clone's RETEST called -- so it is NA in exactly the rows a trunk event
      # adds, and plot_psbulk's ifelse(p_up > 0.5, ...) would return NA and drop
      # those points entirely (na.translate = FALSE). pBAF above/below 0.5 is
      # the mirrored band the panel draws anyway.
      p_up = ifelse(is.na(p_up),
                    ifelse(is.na(pBAF), 0.5, as.numeric(pBAF > 0.5)),
                    p_up),
      # plot_psbulk only rewrites state_post when cnv_state_post is amp/loh/del,
      # so rows going neutral must be cleared here or they keep their old colour.
      state_post     = ifelse(geno_state == "neu", "neu", state_post),
      cnv_state_post = geno_state,
      cnv_state      = geno_state) |>
    dplyr::select(-geno_state, -.clone, -.chrom)
}

#' Colour a pseudobulk table by the consensus segment states
#'
#' Every clone is coloured identically, by the consensus call. Makes the panel
#' agree with the heatmap, at the cost of asserting each call on every clone --
#' including the root/normal clone. Use `numbat_clone_geno_bulk()` for per-clone
#' carriage.
#'
#' @param bulk A `bulk_clones` table carrying `seg_cons`.
#' @param segs The `segs_consensus` table for the round `bulk` was retested
#'   against: `bulk_clones_<k>` goes with `segs_consensus_<k>`.
#' @return `bulk`, recoloured.
#' @export
numbat_consensus_bulk <- function(bulk, segs) {

  if (!is.data.frame(bulk) || !is.data.frame(segs)) return(bulk)
  if (!"seg_cons" %in% names(bulk) || !"seg_cons" %in% names(segs)) return(bulk)
  if (!all(c("cnv_state", "cnv_state_post") %in% names(segs)))     return(bulk)

  cons <- segs |>
    dplyr::distinct(seg_cons, cnv_state, cnv_state_post) |>
    dplyr::transmute(
      seg_cons,
      cons_state = ifelse(cnv_state == "neu", cnv_state, cnv_state_post))

  bulk |>
    dplyr::left_join(cons, by = "seg_cons") |>
    dplyr::mutate(
      cons_state = ifelse(is.na(cons_state), cnv_state_post, cons_state),
      p_up = ifelse(is.na(p_up),
                    ifelse(is.na(pBAF), 0.5, as.numeric(pBAF > 0.5)),
                    p_up),
      state_post     = ifelse(cons_state == "neu", "neu", state_post),
      cnv_state_post = cons_state,
      cnv_state      = cons_state) |>
    dplyr::select(-cons_state)
}
