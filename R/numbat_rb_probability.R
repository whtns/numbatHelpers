# Per-cell SCNA posterior for the canonical RB events.
#
# plot_variability_at_SCNA() (plot_functions_1.R) already draws probability on y
# against cells on x, faceted by segment, with a clone ribbon underneath. It is
# normally fed numbat_heatmap[[3]][["data"]] from make_numbat_heatmaps(), which
# needs a Seurat object, a cluster dictionary and clone simplifications.
#
# This wrapper feeds it the same panel taken straight off the numbat object, and
# restricts the facets to the segments that actually support a canonical RB
# event -- so the plot answers "in what fraction of cells, and how confidently,
# is this sample's 1q gain present" rather than showing every segment numbat
# called. That keeps it consistent with the QC panels, which apply the same
# arm-coverage floor.

#' Segments supporting each canonical RB event
#'
#' @param segs A `segs_consensus` table.
#' @param min_arm_frac Arm-coverage floor; events below it are not returned.
#' @return data.frame of `seg`, `CHROM` and the `event` each supports. Empty if
#'   the sample carries no canonical event above the floor.
#' @export
numbat_rb_segments <- function(segs, min_arm_frac = RB_MIN_ARM_FRAC) {
  empty <- data.frame(seg = character(0), CHROM = character(0),
                      event = character(0), stringsAsFactors = FALSE)
  if (is.null(segs) || !is.data.frame(segs) || nrow(segs) == 0) return(empty)
  segs <- as.data.frame(segs)
  if (!all(c("CHROM", "seg", "seg_start", "seg_end") %in% names(segs))) return(empty)
  state <- if ("cnv_state_post" %in% names(segs)) segs$cnv_state_post else segs$cnv_state
  if (is.null(state)) return(empty)

  keep_ev <- numbat_rb_events(segs, min_arm_frac = min_arm_frac)
  if (!length(keep_ev)) return(empty)

  out <- lapply(which(.RB_WINDOWS$event %in% keep_ev), function(i) {
    w   <- .RB_WINDOWS[i, ]
    hit <- as.character(segs$CHROM) == w$CHROM &
      grepl(if (w$want == "gain") "amp" else "del|loh", state) &
      segs$seg_end > w$xmin * 1e6 & segs$seg_start < w$xmax * 1e6
    if (!any(hit, na.rm = TRUE)) return(NULL)
    data.frame(seg = as.character(segs$seg[which(hit)]),
               CHROM = as.character(segs$CHROM[which(hit)]),
               event = w$event, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) return(empty)
  unique(out)
}

#' Per-cell SCNA probability at the canonical RB events of one sample
#'
#' Probability on y, cells on x, one facet per supporting segment, coloured by
#' CNV state with the clone assignment as a ribbon beneath -- the same figure
#' [plot_variability_at_SCNA()] produces, but sourced from the numbat object
#' alone and restricted to the canonical RB segments.
#'
#' @param nb A numbat object (R6 `Numbat`).
#' @param sample_id Label for the title; taken from `nb$label` when NULL.
#' @param min_arm_frac Arm-coverage floor passed to [numbat_rb_segments()].
#' @param p_min Where to draw the reference line.
#' @return A ggplot, or NULL when the object carries no canonical event above
#'   the floor or the heatmap panel cannot be built.
#' @export
plot_rb_scna_probability <- function(nb, sample_id = NULL,
                                     min_arm_frac = RB_MIN_ARM_FRAC,
                                     p_min = 0.9) {
  if (is.null(sample_id)) sample_id <- tryCatch(nb$label, error = function(e) "")

  segs <- tryCatch(nb[["segs_consensus"]], error = function(e) NULL)
  rb   <- numbat_rb_segments(segs, min_arm_frac = min_arm_frac)
  if (nrow(rb) == 0) return(NULL)

  # Panel 3 of the phylo heatmap is the per-cell posterior table. Built here
  # rather than via make_numbat_heatmaps() so no Seurat object is required.
  d <- tryCatch(nb$plot_phylo_heatmap(raster = FALSE, show_phylo = TRUE)[[3]]$data,
                error = function(e) NULL)
  if (is.null(d) || !all(c("cell", "seg", "p_cnv", "cnv_state") %in% names(d))) return(NULL)
  d <- as.data.frame(d)

  # plot_variability_at_SCNA draws a clone ribbon from clone_opt, which is not
  # in that panel; it lives in clone_post and joins on cell.
  cp <- tryCatch(as.data.frame(nb[["clone_post"]]), error = function(e) NULL)
  d$clone_opt <- if (!is.null(cp) && all(c("cell", "clone_opt") %in% names(cp))) {
    cp$clone_opt[match(d$cell, cp$cell)]
  } else NA

  # DO NOT match segs_consensus$seg against this panel's seg. numbat relabels
  # segments between the two tables -- for SRX10831280 the consensus calls chr16
  # 16b/16d/16f/16i while the per-cell table calls the same regions
  # 16b/16c/16e/16g -- so a name join silently drops facets. Match on genomic
  # coordinates instead, which are shared.
  if (!all(c("CHROM", "seg_start", "seg_end") %in% names(d))) return(NULL)
  keep_ev <- unique(rb$event)
  d$event <- NA_character_
  for (i in which(.RB_WINDOWS$event %in% keep_ev)) {
    w   <- .RB_WINDOWS[i, ]
    hit <- as.character(d$CHROM) == w$CHROM &
      grepl(if (w$want == "gain") "amp" else "del|loh", d$cnv_state) &
      d$seg_end > w$xmin * 1e6 & d$seg_start < w$xmax * 1e6
    d$event[which(hit)] <- w$event
  }
  d <- d[!is.na(d$event), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  # Facet by event and segment extent, so the arm is named and the reader can
  # see that an event may rest on more than one segment.
  d$seg <- sprintf("%s %s  (%s: %.1f-%.1f Mb)",
                   sub("_.*", "", d$event), sub(".*_", "", d$event),
                   as.character(d$seg), d$seg_start / 1e6, d$seg_end / 1e6)
  d$seg <- factor(d$seg, levels = unique(d$seg[order(d$event, d$seg_start)]))

  ttl <- sprintf("%s - per-cell SCNA posterior at canonical RB events (%d cells)",
                 sample_id, length(unique(d$cell)))
  p <- plot_variability_at_SCNA(d, p_min = p_min)
  p + ggplot2::labs(title = ttl,
                    subtitle = sprintf("dashed line p = %.2f; x is cells ordered as in the phylo heatmap", p_min))
}
