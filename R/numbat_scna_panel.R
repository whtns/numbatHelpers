# A readable view of the SCNAs a numbat object actually called.
#
# numbat's own two plots are not enough to inspect an SCNA. plot_phylo_heatmap
# resolves cells beautifully but its x axis is whole chromosomes, so 1q and 1p
# are indistinguishable; plot_consensus labels segments "a", "b", "c'" with no
# coordinates, no arm, and no state names. Neither answers "which arm, how big,
# how strong, and in what fraction of cells".
#
# This draws one panel per chromosome with a real Mb axis, the centromere
# marked, each segment labelled with its state, LLR and cell fraction, and the
# five canonical RB windows shaded so a missing call is as visible as a present
# one.

# hg38 centromere midpoints (Mb). These are the same boundaries the canonical RB
# event definitions use -- 1q > 125, 2p < 93, 6p < 59, 16q > 36.8 -- so arm calls
# here and event calls in numbat_rb_events() cannot disagree.
.HG38_CEN <- c(
  `1` = 123.4, `2` = 93.9, `3` = 90.9, `4` = 50.0, `5` = 48.8, `6` = 59.8,
  `7` = 60.1, `8` = 45.2, `9` = 43.0, `10` = 39.8, `11` = 53.4, `12` = 35.5,
  `13` = 17.7, `14` = 17.2, `15` = 19.0, `16` = 36.8, `17` = 25.1, `18` = 18.5,
  `19` = 26.2, `20` = 28.1, `21` = 12.0, `22` = 15.0, X = 60.6, Y = 10.4
)

.HG38_LEN <- c(
  `1` = 248.9, `2` = 242.2, `3` = 198.3, `4` = 190.2, `5` = 181.5, `6` = 170.8,
  `7` = 159.3, `8` = 145.1, `9` = 138.4, `10` = 133.8, `11` = 135.1, `12` = 133.3,
  `13` = 114.4, `14` = 107.0, `15` = 102.0, `16` = 90.3, `17` = 83.3, `18` = 80.4,
  `19` = 58.6, `20` = 64.4, `21` = 46.7, `22` = 50.8, X = 156.0, Y = 57.2
)

# Chromosomes carrying a canonical RB event, and the window on each. Always
# drawn, whether or not the sample calls anything there.
.RB_WINDOWS <- data.frame(
  CHROM = c("1", "2", "6", "13", "16"),
  arm   = c("1q", "2p", "6p", "13q", "16q"),
  event = c("1q_gain", "2p_gain", "6p_gain", "13q_loss", "16q_loss"),
  want  = c("gain", "gain", "gain", "loss", "loss"),
  xmin  = c(123.4, 0, 0, 17.7, 36.8),
  xmax  = c(248.9, 93.9, 59.8, 114.4, 90.3),
  stringsAsFactors = FALSE
)

# Total length covered by a union of intervals. segs_consensus carries one row
# per consensus component, so supporting segments routinely overlap or repeat --
# summing their lengths overcounts and can exceed the arm.
.union_len <- function(lo, hi) {
  if (!length(lo)) return(0)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  tot <- 0; cs <- lo[1]; ce <- hi[1]
  for (i in seq_along(lo)[-1]) {
    if (lo[i] > ce) { tot <- tot + (ce - cs); cs <- lo[i]; ce <- hi[i] }
    else ce <- max(ce, hi[i])
  }
  tot + (ce - cs)
}

#' Fraction of each canonical RB arm covered by a supporting call
#'
#' @param segs A `segs_consensus` table.
#' @return Named numeric, one entry per canonical event, each the union extent of
#'   the segments supporting it divided by the length of the target arm.
#' @export
numbat_rb_arm_frac <- function(segs) {
  out <- stats::setNames(rep(0, nrow(.RB_WINDOWS)), .RB_WINDOWS$event)
  if (is.null(segs) || !is.data.frame(segs) || nrow(segs) == 0) return(out)
  segs <- as.data.frame(segs)
  if (!all(c("CHROM", "seg_start", "seg_end") %in% names(segs))) return(out)
  state <- if ("cnv_state_post" %in% names(segs)) segs$cnv_state_post else segs$cnv_state
  if (is.null(state)) return(out)

  d <- unique(data.frame(
    CHROM = as.character(segs$CHROM),
    lo    = as.numeric(segs$seg_start),
    hi    = as.numeric(segs$seg_end),
    state = as.character(state),
    stringsAsFactors = FALSE
  ))

  for (i in seq_len(nrow(.RB_WINDOWS))) {
    w    <- .RB_WINDOWS[i, ]
    alo  <- w$xmin * 1e6; ahi <- w$xmax * 1e6
    keep <- d$CHROM == w$CHROM &
      grepl(if (w$want == "gain") "amp" else "del|loh", d$state) &
      d$hi > alo & d$lo < ahi
    if (!any(keep, na.rm = TRUE)) next
    sub <- d[which(keep), , drop = FALSE]
    out[w$event] <- .union_len(pmax(sub$lo, alo), pmin(sub$hi, ahi)) / (ahi - alo)
  }
  out
}

# numbat's own state palette, so this panel and the heatmap agree on colour.
.CNV_PAL <- c(amp = "darkred", bamp = "salmon", del = "royalblue",
              bdel = "darkblue", loh = "darkgreen", neu = "gray90")

.arm_of <- function(chrom, start_mb, end_mb) {
  cen <- .HG38_CEN[as.character(chrom)]
  ifelse(is.na(cen), "",
    ifelse(end_mb   <= cen, paste0(chrom, "p"),
    ifelse(start_mb >= cen, paste0(chrom, "q"),
           paste0(chrom, "p-q"))))
}

#' Tidy the segments a numbat object called, with arms and cell fractions
#'
#' @param nb A numbat object (R6 `Numbat` or the recovered-list fallback).
#' @param non_neutral_only Drop `neu` segments.
#' @return data.frame with one row per consensus segment: arm, Mb coordinates,
#'   state, LLR, and `cell_frac` (fraction of cells with posterior p_cnv > 0.5).
#' @export
numbat_scna_table <- function(nb, non_neutral_only = TRUE) {
  segs <- tryCatch(nb[["segs_consensus"]], error = function(e) NULL)
  if (is.null(segs) || !is.data.frame(segs) || nrow(segs) == 0) return(NULL)
  segs <- as.data.frame(segs)

  state <- if ("cnv_state_post" %in% names(segs)) segs$cnv_state_post else segs$cnv_state
  d <- data.frame(
    CHROM    = as.character(segs$CHROM),
    seg      = as.character(segs$seg),
    start_Mb = round(segs$seg_start / 1e6, 2),
    end_Mb   = round(segs$seg_end   / 1e6, 2),
    state    = as.character(state),
    LLR      = round(suppressWarnings(as.numeric(segs$LLR)), 1),
    n_genes  = segs$n_genes,
    n_snps   = segs$n_snps,
    stringsAsFactors = FALSE
  )
  # segs_consensus carries one row per consensus component, so the same segment
  # can appear more than once with identical coordinates and state.
  d <- unique(d)
  d$size_Mb <- round(d$end_Mb - d$start_Mb, 1)
  d$arm     <- .arm_of(d$CHROM, d$start_Mb, d$end_Mb)

  # Fraction of cells carrying the segment, from the per-cell posterior.
  jp <- tryCatch(nb[["joint_post"]], error = function(e) NULL)
  d$cell_frac <- NA_real_
  if (is.data.frame(jp) && all(c("seg", "p_cnv") %in% names(jp))) {
    jp <- as.data.frame(jp)
    fr <- tapply(jp$p_cnv > 0.5, jp$seg, mean, na.rm = TRUE)
    d$cell_frac <- round(unname(fr[d$seg]), 3)
  }

  # rb_event is only stamped when the event clears the arm-coverage floor, so a
  # focal speck is still listed as a segment but is not labelled as the canonical
  # arm-level SCNA. arm_frac carries the coverage that decision was made on.
  fr <- numbat_rb_arm_frac(segs)
  d$rb_event <- ""
  d$arm_frac <- NA_real_
  for (i in seq_len(nrow(.RB_WINDOWS))) {
    w    <- .RB_WINDOWS[i, ]
    hit  <- d$CHROM == w$CHROM &
      grepl(if (w$want == "gain") "amp" else "del|loh", d$state) &
      d$end_Mb > w$xmin & d$start_Mb < w$xmax
    if (!any(hit)) next
    d$arm_frac[hit] <- round(fr[[w$event]], 3)
    if (fr[[w$event]] >= RB_MIN_ARM_FRAC) {
      d$rb_event[hit] <- paste0(w$arm, if (w$want == "gain") "+" else "-")
    }
  }

  if (non_neutral_only) d <- d[d$state != "neu", , drop = FALSE]
  d <- d[order(suppressWarnings(as.integer(d$CHROM)), d$start_Mb), , drop = FALSE]
  rownames(d) <- NULL
  d[, c("CHROM", "arm", "seg", "start_Mb", "end_Mb", "size_Mb",
        "state", "LLR", "cell_frac", "n_genes", "n_snps", "arm_frac", "rb_event")]
}

#' Per-chromosome SCNA map with arms, coordinates, states and cell fractions
#'
#' One facet per chromosome, real Mb axis, centromere dashed, segments coloured
#' by state and labelled `state LLR=.. cells=..%`. The five canonical RB windows
#' are shaded on chr1/2/6/13/16 and those chromosomes are always shown, so an
#' absent RB call reads as clearly as a present one.
#'
#' @param nb A numbat object.
#' @param title Panel title.
#' @param ncol Facet columns.
#' @return A ggplot, or NULL if the object has no usable segmentation.
#' @export
plot_numbat_scna_map <- function(nb, title = "", ncol = 4) {
  all_segs <- numbat_scna_table(nb, non_neutral_only = FALSE)
  if (is.null(all_segs) || nrow(all_segs) == 0) return(NULL)

  # Show every chromosome that calls something, plus the RB chromosomes always.
  keep <- union(unique(all_segs$CHROM[all_segs$state != "neu"]), .RB_WINDOWS$CHROM)
  d    <- all_segs[all_segs$CHROM %in% keep, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  chrom_levels <- as.character(sort(unique(suppressWarnings(as.integer(d$CHROM)))))
  chrom_levels <- c(chrom_levels, setdiff(unique(d$CHROM), chrom_levels))
  d$CHROM <- factor(d$CHROM, levels = chrom_levels)

  nn <- d[d$state != "neu", , drop = FALSE]
  nn$lab <- sprintf("%s %s  LLR=%s%s", nn$arm, nn$state, nn$LLR,
                    ifelse(is.na(nn$cell_frac), "",
                           sprintf("  cells=%d%%", round(100 * nn$cell_frac))))

  # Full chromosome extent, so a segment's size is read against the chromosome
  # rather than against whatever happened to be called.
  ext <- data.frame(
    CHROM = factor(chrom_levels, levels = chrom_levels),
    len   = unname(.HG38_LEN[chrom_levels]),
    stringsAsFactors = FALSE
  )
  ext$len[is.na(ext$len)] <- vapply(
    as.character(ext$CHROM[is.na(ext$len)]),
    function(cc) max(d$end_Mb[d$CHROM == cc], na.rm = TRUE), numeric(1))

  cen <- data.frame(
    CHROM = factor(chrom_levels, levels = chrom_levels),
    x     = unname(.HG38_CEN[chrom_levels]),
    stringsAsFactors = FALSE
  )
  cen <- cen[!is.na(cen$x), , drop = FALSE]

  rb <- .RB_WINDOWS[.RB_WINDOWS$CHROM %in% chrom_levels, , drop = FALSE]
  if (nrow(rb)) rb$CHROM <- factor(rb$CHROM, levels = chrom_levels)

  p <- ggplot2::ggplot() +
    # RB windows first, underneath everything.
    { if (nrow(rb))
        ggplot2::geom_rect(
          data = rb,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0.05, ymax = 0.95),
          fill = "gold", alpha = 0.18, inherit.aes = FALSE) } +
    # Chromosome backbone.
    ggplot2::geom_rect(
      data = ext,
      ggplot2::aes(xmin = 0, xmax = len, ymin = 0.35, ymax = 0.65),
      fill = "gray93", colour = "gray70", linewidth = 0.2, inherit.aes = FALSE) +
    # Called segments.
    ggplot2::geom_rect(
      data = d[d$state != "neu", , drop = FALSE],
      ggplot2::aes(xmin = start_Mb, xmax = end_Mb, ymin = 0.35, ymax = 0.65,
                   fill = state),
      colour = NA, inherit.aes = FALSE) +
    ggplot2::geom_vline(data = cen, ggplot2::aes(xintercept = x),
                        linetype = "dashed", colour = "gray45", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = .CNV_PAL, name = "CNV state") +
    ggplot2::facet_wrap(~ CHROM, scales = "free_x", ncol = ncol) +
    ggplot2::scale_y_continuous(limits = c(0, 1.35), expand = c(0, 0)) +
    ggplot2::labs(
      title    = title,
      subtitle = paste("shaded = canonical RB window (1q+, 2p+, 6p+, 13q-, 16q-);",
                       "dashed = centromere; chr1/2/6/13/16 always shown"),
      x = "position (Mb)", y = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_blank(),
      axis.ticks.y    = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray92", colour = NA),
      strip.text      = ggplot2::element_text(face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = 7, colour = "gray35"),
      legend.position = "bottom"
    )

  if (nrow(nn)) {
    p <- p + ggrepel::geom_text_repel(
      data = nn,
      ggplot2::aes(x = (start_Mb + end_Mb) / 2, y = 0.68, label = lab),
      size = 2.1, direction = "y", nudge_y = 0.18, segment.size = 0.2,
      segment.colour = "gray55", min.segment.length = 0, box.padding = 0.12,
      max.overlaps = Inf, inherit.aes = FALSE)
  }
  p
}

#' Render the SCNA table as a page-sized graphic
#'
#' @param scna_tbl Output of [numbat_scna_table()].
#' @param title Title above the table.
#' @return A gtable, or NULL when there is nothing to show.
#' @export
plot_numbat_scna_table <- function(scna_tbl, title = "") {
  if (is.null(scna_tbl) || nrow(scna_tbl) == 0) return(NULL)
  th <- gridExtra::ttheme_minimal(
    base_size = 7,
    core    = list(fg_params = list(hjust = 1, x = 0.95)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
  gridExtra::arrangeGrob(
    gridExtra::tableGrob(scna_tbl, rows = NULL, theme = th),
    top = grid::textGrob(title, gp = grid::gpar(fontsize = 11, fontface = "bold"))
  )
}
