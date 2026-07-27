# Base-graphics PDF renderer for mc_report. No new dependencies: every page is a
# plot.new() on the pdf device, text drawn in a monospaced font, tables printed
# through print.data.frame and flowed line-by-line with automatic pagination.

# a paginating writer over the current device. emit() draws lines top-to-bottom
# and starts a new page when the current one fills.
rp_writer <- function(title) {
  st <- new.env(parent = emptyenv())
  st$lines <- 62L
  st$i <- st$lines # force a new page on the first emit
  st$top <- 0.94
  st$bot <- 0.05
  st$left <- 0.045

  page <- function() {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::text(
      st$left,
      0.975,
      title,
      adj = c(0, 0.5),
      font = 2,
      cex = 0.8,
      family = "sans"
    )
    graphics::segments(st$left, 0.958, 1 - st$left, 0.958, col = "grey75")
    st$i <- 0L
  }
  yfor <- function(i) st$top - i * (st$top - st$bot) / st$lines

  emit <- function(txt = "", font = 1L, cex = 0.6, col = "black", family = "mono") {
    for (t in txt) {
      if (st$i >= st$lines) {
        page()
      }
      graphics::text(
        st$left,
        yfor(st$i),
        t,
        adj = c(0, 1),
        font = font,
        cex = cex,
        col = col,
        family = family
      )
      st$i <- st$i + 1L
    }
    invisible()
  }

  heading <- function(txt) {
    emit("")
    emit(txt, font = 2L, cex = 0.85, family = "sans")
    emit(strrep("-", 92L), cex = 0.6, col = "grey60")
  }

  list(
    emit = emit,
    heading = heading,
    newpage = function() st$i <- st$lines
  )
}

# wrap a long string to a width in characters (for cli-style multi-line notes)
rp_wrap <- function(txt, width = 92L) {
  txt <- gsub("\\s+", " ", trimws(txt))
  strwrap(txt, width = width)
}

# a data.frame as aligned monospaced lines, numeric columns rounded, rows capped
rp_table_lines <- function(df, max_rows = 55L, digits = 3L) {
  if (is.null(df) || !nrow(df)) {
    return("  (none)")
  }
  disp <- df
  for (nm in names(disp)) {
    col <- disp[[nm]]
    if (is.numeric(col) && !is.integer(col)) {
      disp[[nm]] <- ifelse(
        is.finite(col),
        formatC(col, format = "f", digits = digits),
        "NA"
      )
    }
  }
  truncated <- nrow(disp) > max_rows
  if (truncated) {
    disp <- disp[seq_len(max_rows), , drop = FALSE]
  }
  lines <- utils::capture.output(print(disp, row.names = FALSE, right = FALSE))
  if (truncated) {
    lines <- c(lines, sprintf("... and %d more row(s)", nrow(df) - max_rows))
  }
  lines
}

# capped list of ids as a comma-wrapped block
rp_id_block <- function(ids, cap = 30L, width = 90L) {
  if (!length(ids)) {
    return("  (none)")
  }
  shown <- utils::head(ids, cap)
  txt <- paste(shown, collapse = ", ")
  if (length(ids) > cap) {
    txt <- sprintf("%s, ... (%d total)", txt, length(ids))
  }
  paste0("  ", strwrap(txt, width = width))
}

# one plot on its own page, confined to a sub-region (NDC box `fig`) so it is
# not stretched to the whole page; caller advances the writer afterwards
rp_plot_page <- function(fun, fig = c(0.10, 0.90, 0.54, 0.94), mar = c(4, 4, 2.5, 1)) {
  # order matters: setting mfrow resets fig, so clear any leftover grid FIRST,
  # then set fig + mar BEFORE plot.new (the writer leaves mar = 0, which
  # otherwise defeats fig); then reuse that frame (new = TRUE) for the plot
  graphics::par(mfrow = c(1, 1))
  graphics::par(fig = fig, mar = mar)
  graphics::plot.new()
  graphics::par(new = TRUE)
  fun()
  graphics::par(fig = c(0, 1, 0, 1), new = FALSE)
  invisible()
}

# ---- section drawers -------------------------------------------------------

rp_draw_overview <- function(w, rep) {
  m <- rep[["meta"]]
  kinds <- paste(c(
    if (isTRUE(m[["has_dnam"]])) "DNAm input QC",
    if (isTRUE(m[["has_score"]])) "score QC"
  ), collapse = " + ")
  w$emit("")
  w$emit("methylCIPHER QC report", font = 2L, cex = 1.4, family = "sans")
  w$emit("")
  w$emit(sprintf("Report type : %s", kinds))
  w$emit(sprintf("Samples     : %s", m[["n_samples"]]))
  if (isTRUE(m[["has_dnam"]])) {
    w$emit(sprintf("Probes      : %s", m[["n_probes"]]))
  }
  vd <- m[["verdict"]]
  if (!is.null(vd)) {
    w$emit("")
    w$emit(sprintf("Overall verdict : %s", vd[["overall"]]), font = 2L)
    for (nm in names(vd[["sections"]])) {
      w$emit(sprintf("  %-9s : %s", nm, vd[["sections"]][[nm]]))
    }
  }
}

rp_draw_dnam <- function(w, d) {
  # format
  w$heading("1. Data format")
  fmt <- d[["format"]]
  w$emit(sprintf("Dimensions: %d samples x %d probes.", fmt[["n_samples"]], fmt[["n_probes"]]))
  if (fmt[["ok"]]) {
    w$emit("No formatting problems detected.", col = "grey30")
  } else {
    w$emit("Issues:", font = 2L)
    for (iss in fmt[["issues"]]) {
      w$emit(paste0("  - ", rp_wrap(iss, 88L)))
    }
  }

  # array
  w$heading("2. Array")
  a <- d[["array"]]
  w$emit(sprintf("Best guess : %s", a[["guess"]]))
  w$emit(sprintf("cg probes  : %d of %d columns", a[["n_cg"]], a[["n_probes"]]))
  for (ln in rp_wrap(a[["note"]], 88L)) {
    w$emit(ln, col = "grey30")
  }

  # beta range
  w$heading("3. Beta values")
  b <- d[["beta"]]
  w$emit(sprintf("Range          : [%.4g, %.4g]", b[["min"]], b[["max"]]))
  w$emit(sprintf("Below 0        : %d value(s)", b[["n_below0"]]))
  w$emit(sprintf("Above 1        : %d value(s)", b[["n_above1"]]))
  w$emit(sprintf("NaN / Inf      : %d / %d", b[["n_nan"]], b[["n_inf"]]))
  w$emit(sprintf("Out-of-range   : %d probe(s)", b[["n_out_of_range"]]))
  for (ln in rp_wrap(b[["scale_note"]], 88L)) {
    w$emit(ln, col = "grey30")
  }
  if (b[["n_out_of_range"]] > 0L) {
    w$emit("Out-of-range probes:")
    for (ln in rp_id_block(b[["out_of_range_probes"]])) {
      w$emit(ln)
    }
  }

  # missingness
  w$heading("4. Missing probes")
  ms <- d[["missing"]]
  w$emit(sprintf("Probes needed across assessed clocks : %d", ms[["total_needed_union"]]))
  w$emit(sprintf("Absent from the matrix entirely      : %d", ms[["n_absent_from_matrix"]]))
  w$emit(sprintf("Present but all-NA (fully missing)    : %d", ms[["n_full_na_probes"]]))
  w$emit(sprintf("Present with some NA (partial)        : %d", ms[["n_partial_probes"]]))
  if (ms[["n_full_na_probes"]] > 0L) {
    w$emit("Fully-missing probes:")
    for (ln in rp_id_block(ms[["full_na_probes"]])) {
      w$emit(ln)
    }
  }
  if (ms[["n_partial_probes"]] > 0L) {
    w$emit("Partially-missing probes:")
    for (ln in rp_id_block(ms[["partial_probes"]])) {
      w$emit(ln)
    }
  }

  # per-sample QC
  w$heading("5. Per-sample QC")
  sm <- d[["samples"]]
  w$emit(sprintf(
    "%d sample(s); %d flagged as cohort outliers (missingness, mean-beta, distribution).",
    nrow(sm[["table"]]),
    sm[["n_flagged"]]
  ))
  if (!isTRUE(sm[["cohort_bimodal_ok"]])) {
    for (ln in rp_wrap(sprintf(
      "Cohort-wide: betas do not look bimodal -- the middle of the distribution
       sits at %.0f%% of the flanking modes (a bimodal dataset is well under
       60%%). The data may be mis-scaled or over-smoothed; check that these are
       raw beta values.",
      100 * sm[["cohort_bimodal_ratio"]]
    ), 88L)) {
      w$emit(ln, col = "firebrick")
    }
  }
  flagged <- sm[["table"]][nzchar(sm[["table"]][["flags"]]), , drop = FALSE]
  if (nrow(flagged)) {
    w$emit("")
    w$emit("Flagged samples:", font = 2L)
    for (ln in rp_table_lines(flagged, max_rows = 60L)) {
      w$emit(ln, cex = 0.55)
    }
  } else {
    w$emit("No samples flagged.", col = "grey30")
  }
  dups <- sm[["duplicates"]]
  w$emit("")
  if (nrow(dups)) {
    w$emit("Near-identical sample pairs (cor >= 0.999):", font = 2L)
    for (ln in rp_table_lines(dups, max_rows = 40L)) {
      w$emit(ln, cex = 0.55)
    }
  } else {
    w$emit("No near-identical sample pairs.", col = "grey30")
  }

  # per-clock coverage
  w$heading("6. Per-clock probe coverage")
  if (length(d[["skipped_external"]])) {
    w$emit(sprintf(
      "Not assessed (external pack not loaded): %s",
      paste(d[["skipped_external"]], collapse = ", ")
    ))
    w$emit("Pass assets= or pre-download the pack to include these.", col = "grey30")
    w$emit("")
  }
  cov <- d[["coverage"]]
  disp <- data.frame(
    clock_id = cov[["clock_id"]],
    needed = cov[["needed"]],
    present = cov[["present"]],
    absent = cov[["absent"]],
    all_na = cov[["all_na"]],
    cov_pct = as.integer(round(100 * cov[["coverage"]])),
    ext = ifelse(cov[["external"]], "y", ""),
    stringsAsFactors = FALSE
  )
  for (ln in rp_table_lines(disp, max_rows = 200L)) {
    w$emit(ln, cex = 0.55)
  }
}

rp_draw_score <- function(w, s) {
  w$heading("Score distribution")
  for (ln in rp_wrap(s[["note"]], 88L)) {
    w$emit(ln, col = "grey30")
  }
  w$emit("")
  for (ln in rp_table_lines(s[["summary"]], max_rows = 120L)) {
    w$emit(ln, cex = 0.55)
  }

  if (!is.null(s[["flags"]]) && nrow(s[["flags"]])) {
    w$heading("Flagged vs reference (|z| > 2)")
    for (ln in rp_table_lines(s[["flags"]], max_rows = 120L)) {
      w$emit(ln, cex = 0.55)
    }
  }

  # out-of-training-range flags (needs a reference with expect_lo/expect_hi)
  if (!is.null(s[["range_flags"]]) && nrow(s[["range_flags"]])) {
    w$heading("Scores outside the expected/training range")
    for (ln in rp_table_lines(s[["range_flags"]], max_rows = 120L)) {
      w$emit(ln, cex = 0.55)
    }
  }

  # samples that failed to score
  w$heading("Samples with NA scores")
  na <- s[["na_samples"]]
  if (!nrow(na)) {
    w$emit("Every sample scored for every clock.", col = "grey30")
  } else {
    w$emit(sprintf("%d sample(s) have >= 1 NA score:", nrow(na)))
    for (ln in rp_table_lines(na, max_rows = 60L)) {
      w$emit(ln, cex = 0.55)
    }
  }

  # age-association check vs the shipped reference (advisory)
  a <- s[["associations"]]
  if (!is.null(a)) {
    w$heading("Age association vs reference (advisory)")
    for (ln in rp_wrap(a[["note"]], 88L)) {
      w$emit(ln, col = "grey30")
    }
    w$emit("")
    off <- a[["table"]][a[["table"]][["outside"]] | a[["table"]][["wrong_sign"]], , drop = FALSE]
    if (nrow(off)) {
      w$emit(sprintf(
        "%d of %d checked clock(s) fall outside the expected age-correlation range:",
        nrow(off), nrow(a[["table"]])
      ))
      disp <- data.frame(
        clock = off[["clock"]],
        obs_r = off[["obs_age_r"]],
        exp_r = off[["exp_age_r"]],
        exp_lo = off[["exp_lo"]],
        exp_hi = off[["exp_hi"]],
        sign = ifelse(off[["wrong_sign"]], "wrong", ""),
        stringsAsFactors = FALSE
      )
      for (ln in rp_table_lines(disp, max_rows = 120L)) {
        w$emit(ln, cex = 0.55)
      }
    } else {
      w$emit(sprintf(
        "All %d checked clock(s) track age within the expected range.",
        nrow(a[["table"]])
      ), col = "grey30")
    }
  }

  w$heading("Clock coverage summary")
  cc <- s[["coverage"]]
  disp <- data.frame(
    clock_id = cc[["clock_id"]],
    role = cc[["role"]],
    score_needed = cc[["score_needed"]],
    score_used = cc[["score_used"]],
    imputed_full = cc[["score_imputed_full"]],
    dropped = cc[["score_dropped"]],
    stringsAsFactors = FALSE
  )
  for (ln in rp_table_lines(disp, max_rows = 120L)) {
    w$emit(ln, cex = 0.55)
  }
}

# coverage histogram page (only when there is something to show)
rp_draw_coverage_plot <- function(w, d) {
  cov <- d[["coverage"]][["coverage"]]
  cov <- cov[is.finite(cov)]
  if (length(cov) < 2L) {
    return(invisible())
  }
  rp_plot_page(function() {
    graphics::hist(
      100 * cov,
      breaks = seq(0, 100, by = 5),
      col = "grey80",
      border = "white",
      main = "Per-clock probe coverage",
      xlab = "Coverage (% of scoring panel present)",
      ylab = "Clocks"
    )
    graphics::abline(v = 75, col = "firebrick", lty = 2)
  })
  w$newpage()
}

# overlaid per-sample beta densities; flagged (non-bimodal) samples in red
rp_draw_sample_density <- function(w, d) {
  dens <- d[["samples"]][["densities"]]
  if (is.null(dens) || !ncol(dens[["y"]])) {
    return(invisible())
  }
  y <- dens[["y"]]
  keep <- which(vapply(seq_len(ncol(y)), function(j) any(is.finite(y[, j])), logical(1L)))
  if (!length(keep)) {
    return(invisible())
  }
  rp_plot_page(function() {
    graphics::matplot(
      dens[["x"]],
      y[, keep, drop = FALSE],
      type = "l",
      lty = 1,
      lwd = 0.7,
      col = ifelse(dens[["flagged"]][keep], "firebrick", "grey60"),
      main = "Per-sample beta distribution",
      xlab = "Beta value",
      ylab = "Density"
    )
    graphics::legend(
      "topright",
      legend = c("typical", "flagged (non-bimodal)"),
      col = c("grey60", "firebrick"),
      lwd = 1.2,
      bty = "n",
      cex = 0.8
    )
  })
  w$newpage()
}

# small-multiple histograms of each clock's scores
rp_draw_score_hists <- function(w, s) {
  scores <- s[["scores"]]
  if (is.null(scores) || !ncol(scores)) {
    return(invisible())
  }
  ids <- colnames(scores)
  show <- utils::head(ids, 12L)
  nc <- min(3L, length(show))
  nr <- ceiling(length(show) / nc)
  # confine the small-multiple grid to the top half of the page via outer margins
  graphics::par(
    fig = c(0, 1, 0, 1),
    mfrow = c(nr, nc),
    omi = c(4.4, 1.0, 0.9, 1.0),
    mar = c(3, 3, 2, 0.6)
  )
  for (id in show) {
    v <- scores[, id]
    v <- v[is.finite(v)]
    graphics::hist(
      v,
      breaks = "FD",
      col = "grey80",
      border = "white",
      main = id,
      xlab = "score",
      cex.main = 0.9
    )
  }
  graphics::mtext("Score distributions", outer = TRUE, cex = 1, line = 1)
  graphics::par(mfrow = c(1, 1), omi = c(0, 0, 0, 0), fig = c(0, 1, 0, 1))
  w$newpage()
}

# clock-clock correlation heatmap (blue -1 / white 0 / red +1)
rp_draw_cor_heatmap <- function(w, s) {
  cm <- s[["correlations"]]
  if (is.null(cm) || ncol(cm) < 2L) {
    return(invisible())
  }
  k <- ncol(cm)
  pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(51)
  cex_ax <- max(0.35, min(0.7, 26 / k))
  rp_plot_page(fig = c(0.14, 0.86, 0.44, 0.93), mar = c(6, 6, 2.5, 1), function() {
    graphics::image(
      1:k,
      1:k,
      t(cm[k:1, , drop = FALSE]),
      zlim = c(-1, 1),
      col = pal,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = "Clock-clock score correlation"
    )
    graphics::mtext("blue = -1   white = 0   red = +1", side = 3, cex = 0.7, col = "grey40")
    graphics::axis(1, at = 1:k, labels = colnames(cm), las = 2, cex.axis = cex_ax)
    graphics::axis(2, at = 1:k, labels = rev(colnames(cm)), las = 2, cex.axis = cex_ax)
  })
  w$newpage()
}

# ---- entry point -----------------------------------------------------------

# render an mc_report to a PDF file; returns the normalized path
render_report_pdf <- function(rep, file = NULL) {
  # default to a findable file in the working directory, not a temp dir that
  # gets wiped at session end; pass file= for anywhere else
  if (is.null(file)) {
    file <- file.path(getwd(), "methylCIPHER-report.pdf")
  }
  checkmate::assert_path_for_output(file, overwrite = TRUE)

  grDevices::pdf(file, width = 8.5, height = 11, title = "methylCIPHER QC report")
  on.exit(grDevices::dev.off(), add = TRUE)

  w <- rp_writer("methylCIPHER QC report")
  rp_draw_overview(w, rep)
  if (!is.null(rep[["dnam"]])) {
    rp_draw_dnam(w, rep[["dnam"]])
    rp_draw_coverage_plot(w, rep[["dnam"]])
    rp_draw_sample_density(w, rep[["dnam"]])
  }
  if (!is.null(rep[["score"]])) {
    rp_draw_score(w, rep[["score"]])
    rp_draw_score_hists(w, rep[["score"]])
    rp_draw_cor_heatmap(w, rep[["score"]])
  }

  normalizePath(file, winslash = "/", mustWork = FALSE)
}
