# report(): one entry point for both QC reports. It inspects what it is given --
# a DNAm matrix (input QC), an mc_result (score QC), or both -- builds a
# structured mc_report, flags headline problems to the console, and writes a PDF.

# expected per-clock score mean+SD reference. Returns NULL until the reference
# table is supplied (see dev/detail-plan.md); the comparison in report_score()
# is already wired to consume a data.frame(clock_id, ref_mean, ref_sd).
mc_score_reference <- function() {
  NULL
}

# shipped per-clock expectation table (meta-analyzed age/sex associations across
# public datasets); NULL if the package was installed without it
mc_clock_reference <- function() {
  path <- system.file("extdata", "clock_reference.csv", package = "methylCIPHERv2")
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# observed age correlation per clock vs the meta-analytic expectation. Advisory
# only: age correlation shrinks with a narrow cohort age range, so out-of-range
# values are a prompt to look, not a QC failure.
score_associations <- function(scores, pheno, pheno_id, ref = mc_clock_reference()) {
  if (is.null(ref) || is.null(pheno) || !("Age" %in% names(pheno))) {
    return(NULL)
  }
  idx <- match(rownames(scores), pheno[[pheno_id]])
  age <- suppressWarnings(as.numeric(pheno[["Age"]][idx]))
  if (sum(is.finite(age)) < 5L) {
    return(NULL)
  }
  ids <- intersect(colnames(scores), ref[["clock"]])
  if (!length(ids)) {
    return(NULL)
  }

  rows <- lapply(ids, function(cl) {
    v <- scores[, cl]
    ok <- is.finite(v) & is.finite(age)
    if (sum(ok) < 5L || stats::sd(v[ok]) == 0) {
      return(NULL)
    }
    obs <- suppressWarnings(stats::cor(v[ok], age[ok]))
    r <- ref[ref[["clock"]] == cl, ]
    outside <- is.finite(r[["age_r_lo"]]) &&
      (obs < r[["age_r_lo"]] | obs > r[["age_r_hi"]])
    wrong_sign <- is.finite(r[["age_r"]]) &&
      abs(r[["age_r"]]) > 0.3 && sign(obs) != sign(r[["age_r"]])
    data.frame(
      clock = cl,
      obs_age_r = round(obs, 3),
      exp_age_r = r[["age_r"]],
      exp_lo = r[["age_r_lo"]],
      exp_hi = r[["age_r_hi"]],
      outside = isTRUE(outside),
      wrong_sign = isTRUE(wrong_sign),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  tbl <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(tbl)) {
    return(NULL)
  }
  rng <- range(age, na.rm = TRUE)
  list(
    table = tbl[order(tbl[["obs_age_r"]] - tbl[["exp_age_r"]]), , drop = FALSE],
    age_range = rng,
    age_sd = stats::sd(age, na.rm = TRUE),
    n_outside = sum(tbl[["outside"]]),
    n_wrong_sign = sum(tbl[["wrong_sign"]]),
    note = sprintf(
      "Advisory. Age correlation depends on the cohort age spread (here %.0f-%.0f, SD %.1f); a narrow range lowers correlation for every clock, so read out-of-range values as a prompt to check, not a failure.",
      rng[1L], rng[2L], stats::sd(age, na.rm = TRUE)
    )
  )
}

# per-clock score distribution, plus deviation flags when a reference exists
report_score <- function(result, score_reference = NULL, pheno = NULL, pheno_id = "ID") {
  check_mc_result(result)
  scores <- result[["scores"]]
  ids <- colnames(scores)

  # fall back to the pheno carried on the record (only has Age if a clock needed it)
  if (is.null(pheno)) {
    pheno <- result[["pheno"]]
    pheno_id <- result[["provenance"]][["pheno_id"]]
  }

  summary <- do.call(rbind, lapply(ids, function(id) {
    v <- scores[, id]
    ok <- v[is.finite(v)]
    q <- if (length(ok)) {
      stats::quantile(ok, c(0.25, 0.5, 0.75), names = FALSE)
    } else {
      rep(NA_real_, 3L)
    }
    data.frame(
      clock_id = id,
      n = length(v),
      n_na = sum(!is.finite(v)),
      mean = if (length(ok)) mean(ok) else NA_real_,
      sd = if (length(ok) > 1L) stats::sd(ok) else NA_real_,
      min = if (length(ok)) min(ok) else NA_real_,
      median = q[2L],
      max = if (length(ok)) max(ok) else NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  # samples that failed to score for one or more clocks (NA / non-finite)
  na_mat <- !is.finite(scores)
  per_sample <- rowSums(na_mat)
  na_samples <- data.frame(
    sample_id = rownames(scores),
    n_na_clocks = as.integer(per_sample),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  na_samples <- na_samples[na_samples[["n_na_clocks"]] > 0L, , drop = FALSE]

  # clock-clock correlation (internal consistency); NULL when only one clock
  correlations <- if (length(ids) > 1L) {
    suppressWarnings(stats::cor(scores, use = "pairwise.complete.obs"))
  } else {
    NULL
  }

  # reference / expectation table (mean+SD and/or plausible score range). Not in
  # yet -- mc_score_reference() is NULL; the comparisons below light up when a
  # data.frame(clock_id, ref_mean, ref_sd, expect_lo, expect_hi) is supplied.
  ref <- if (is.null(score_reference)) mc_score_reference() else score_reference
  flags <- NULL
  range_flags <- NULL
  note <- "No score reference is available yet -- distributions are reported without comparison."
  if (!is.null(ref)) {
    note <- "Observed scores compared against the supplied reference."
    if (all(c("ref_mean", "ref_sd") %in% names(ref))) {
      m <- merge(
        summary[, c("clock_id", "mean", "sd")],
        ref[, c("clock_id", "ref_mean", "ref_sd")],
        by = "clock_id",
        all.x = FALSE
      )
      m[["z_mean"]] <- (m[["mean"]] - m[["ref_mean"]]) / m[["ref_sd"]]
      flags <- m[is.finite(m[["z_mean"]]) & abs(m[["z_mean"]]) > 2, , drop = FALSE]
      flags <- flags[order(-abs(flags[["z_mean"]])), , drop = FALSE]
    }
    if (all(c("expect_lo", "expect_hi") %in% names(ref))) {
      range_flags <- score_range_flags(scores, ref)
    }
  }

  list(
    summary = summary,
    reference = ref,
    flags = flags,
    range_flags = range_flags,
    na_samples = na_samples,
    correlations = correlations,
    associations = score_associations(scores, pheno, pheno_id),
    scores = scores,
    coverage = clocks_coverage(result),
    note = note
  )
}

# per-clock count of scores outside a declared plausible/training range
score_range_flags <- function(scores, ref) {
  ref <- ref[!is.na(ref[["expect_lo"]]) & !is.na(ref[["expect_hi"]]), , drop = FALSE]
  ref <- ref[ref[["clock_id"]] %in% colnames(scores), , drop = FALSE]
  if (!nrow(ref)) {
    return(NULL)
  }
  rows <- lapply(seq_len(nrow(ref)), function(i) {
    id <- ref[["clock_id"]][i]
    lo <- ref[["expect_lo"]][i]
    hi <- ref[["expect_hi"]][i]
    v <- scores[, id]
    n_out <- sum(is.finite(v) & (v < lo | v > hi))
    data.frame(
      clock_id = id,
      expect_lo = lo,
      expect_hi = hi,
      n_out_of_range = as.integer(n_out),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  out <- out[out[["n_out_of_range"]] > 0L, , drop = FALSE]
  if (!nrow(out)) NULL else out[order(-out[["n_out_of_range"]]), , drop = FALSE]
}

# per-section PASS/WARN/FAIL, plus the worst as an overall verdict
report_verdict <- function(rep) {
  rank <- c(PASS = 0L, WARN = 1L, FAIL = 2L)
  v <- character(0)

  d <- rep[["dnam"]]
  if (!is.null(d)) {
    fmt <- d[["format"]]
    v["format"] <- if (!fmt[["ok"]] &&
      (fmt[["transposed"]] || any(grepl("Duplicate", fmt[["issues"]])))) {
      "FAIL"
    } else if (!fmt[["ok"]]) {
      "WARN"
    } else {
      "PASS"
    }

    b <- d[["beta"]]
    v["beta"] <- if (b[["n_nan"]] + b[["n_inf"]] > 0L ||
      grepl("M-values|percentages", b[["scale_note"]])) {
      "FAIL"
    } else if (b[["n_out_of_range"]] > 0L) {
      "WARN"
    } else {
      "PASS"
    }

    # a few unscoreable clocks (e.g. EPIC-only panels on 450K data) is "partial",
    # not failure; FAIL only when coverage is broadly broken (wrong array / bad
    # data), i.e. most clocks cannot be scored.
    cov <- d[["coverage"]][["coverage"]]
    cov <- cov[is.finite(cov)]
    frac_under_75 <- if (length(cov)) mean(cov < 0.75) else 0
    frac_under_50 <- if (length(cov)) mean(cov < 0.5) else 0
    v["coverage"] <- if (frac_under_50 > 0.5) {
      "FAIL"
    } else if (frac_under_75 > 0) {
      "WARN"
    } else {
      "PASS"
    }

    s <- d[["samples"]]
    v["samples"] <- if (s[["n_flagged"]] > 0L ||
      nrow(s[["duplicates"]]) > 0L ||
      !s[["cohort_bimodal_ok"]]) {
      "WARN"
    } else {
      "PASS"
    }
  }

  sc <- rep[["score"]]
  if (!is.null(sc)) {
    bad <- nrow(sc[["na_samples"]]) > 0L ||
      (!is.null(sc[["flags"]]) && nrow(sc[["flags"]]) > 0L) ||
      (!is.null(sc[["range_flags"]]) && nrow(sc[["range_flags"]]) > 0L)
    v["scores"] <- if (bad) "WARN" else "PASS"
  }

  overall <- names(rank)[max(rank[v]) + 1L]
  list(sections = v, overall = overall)
}

# the important flags, surfaced to the console the moment report() runs. Plain
# sprintf strings (no cli {?} markers) so the bullet list needs no qty() juggling.
report_flag_console <- function(rep) {
  issues <- character(0)
  d <- rep[["dnam"]]
  if (!is.null(d)) {
    if (!d[["format"]][["ok"]]) {
      issues <- c(issues, sprintf(
        "Formatting: %d issue(s) -- %s",
        length(d[["format"]][["issues"]]),
        gsub("\\s+", " ", d[["format"]][["issues"]][[1L]])
      ))
    }
    b <- d[["beta"]]
    if (b[["n_out_of_range"]] > 0L || b[["n_nan"]] > 0L || b[["n_inf"]] > 0L) {
      issues <- c(issues, sprintf(
        "Betas: %d below 0, %d above 1 across %d probe(s); %s",
        b[["n_below0"]], b[["n_above1"]], b[["n_out_of_range"]], b[["scale_note"]]
      ))
    }
    low <- sum(d[["coverage"]][["coverage"]] < 0.75, na.rm = TRUE)
    if (low > 0L) {
      issues <- c(issues, sprintf("Coverage: %d clock(s) under 75%% probe coverage.", low))
    }
    sm <- d[["samples"]]
    if (!isTRUE(sm[["cohort_bimodal_ok"]])) {
      issues <- c(issues, "Distribution: betas are not bimodal cohort-wide (possible mis-scaling).")
    }
    if (sm[["n_flagged"]] > 0L) {
      issues <- c(issues, sprintf("Samples: %d flagged as outliers.", sm[["n_flagged"]]))
    }
    if (nrow(sm[["duplicates"]]) > 0L) {
      issues <- c(issues, sprintf(
        "Samples: %d near-identical pair(s) (possible swap/replicate).",
        nrow(sm[["duplicates"]])
      ))
    }
  }
  s <- rep[["score"]]
  if (!is.null(s)) {
    if (nrow(s[["na_samples"]]) > 0L) {
      issues <- c(issues, sprintf("Scores: %d sample(s) with NA scores.", nrow(s[["na_samples"]])))
    }
    if (!is.null(s[["flags"]]) && nrow(s[["flags"]]) > 0L) {
      issues <- c(issues, sprintf("Scores: %d clock(s) deviate from the reference.", nrow(s[["flags"]])))
    }
    if (!is.null(s[["range_flags"]]) && nrow(s[["range_flags"]]) > 0L) {
      issues <- c(issues, sprintf(
        "Scores: %d clock(s) have scores outside the expected range.",
        nrow(s[["range_flags"]])
      ))
    }
  }

  overall <- rep[["meta"]][["verdict"]][["overall"]]
  if (length(issues)) {
    cli::cli_warn(
      c("QC flags (verdict: {overall}):", bullets(issues)),
      call = NULL
    )
  } else {
    cli::cli_inform(c("v" = "QC checks passed (verdict: {overall})."))
  }

  # age-association check is advisory -- reported separately, never in the verdict
  a <- rep[["score"]][["associations"]]
  if (!is.null(a) && (a[["n_outside"]] > 0L || a[["n_wrong_sign"]] > 0L)) {
    msg <- sprintf(
      "Advisory: %d clock(s) track age outside the expected range%s. A narrow cohort age range can cause this.",
      a[["n_outside"]],
      if (a[["n_wrong_sign"]] > 0L) sprintf(" (%d with the wrong sign)", a[["n_wrong_sign"]]) else ""
    )
    cli::cli_inform(c("i" = msg))
  }
  invisible(rep)
}

# report(DNAm, result): either or both. report(a_result) also works positionally.
#' @export
report <- function(
  DNAm = NULL,
  result = NULL,
  file = NULL,
  clocks = "all",
  pheno = NULL,
  pheno_id = "ID",
  assets = NULL,
  ask = TRUE,
  score_reference = NULL
) {
  # let report(my_result) land on `result` without a named argument
  if (inherits(DNAm, "mc_result") && is.null(result)) {
    result <- DNAm
    DNAm <- NULL
  }
  if (is.null(DNAm) && is.null(result)) {
    cli::cli_abort(
      c(
        "report() needs {.arg DNAm}, a {.cls mc_result}, or both.",
        "i" = "Pass a DNAm matrix for an input-QC report, a {.fn calc_clocks}
               result for a score report, or both for both."
      ),
      call = NULL
    )
  }

  dnam_rep <- if (!is.null(DNAm)) report_dnam(DNAm, clocks, assets, ask) else NULL
  score_rep <- if (!is.null(result)) {
    report_score(result, score_reference, pheno, pheno_id)
  } else {
    NULL
  }

  rep <- structure(
    list(
      meta = list(
        has_dnam = !is.null(dnam_rep),
        has_score = !is.null(score_rep),
        n_samples = if (!is.null(DNAm)) nrow(DNAm) else nrow(result[["scores"]]),
        n_probes = if (!is.null(DNAm)) ncol(DNAm) else NA_integer_
      ),
      dnam = dnam_rep,
      score = score_rep
    ),
    class = "mc_report"
  )
  rep[["meta"]][["verdict"]] <- report_verdict(rep)

  report_flag_console(rep)

  path <- render_report_pdf(rep, file)
  rep[["meta"]][["file"]] <- path
  cli::cli_inform(c("v" = "Report written to {.path {path}}."))
  invisible(rep)
}

# console summary of a built report
#' @export
print.mc_report <- function(x, ...) {
  m <- x[["meta"]]
  cat(sprintf(
    "<mc_report> %s\n",
    paste(c(
      if (isTRUE(m[["has_dnam"]])) "DNAm QC",
      if (isTRUE(m[["has_score"]])) "score QC"
    ), collapse = " + ")
  ))
  cat(sprintf("  samples: %s\n", m[["n_samples"]]))
  if (!is.null(m[["verdict"]])) {
    cat(sprintf(
      "  verdict: %s (%s)\n",
      m[["verdict"]][["overall"]],
      paste(sprintf(
        "%s=%s",
        names(m[["verdict"]][["sections"]]),
        m[["verdict"]][["sections"]]
      ), collapse = " ")
    ))
  }

  d <- x[["dnam"]]
  if (!is.null(d)) {
    cat(sprintf(
      "  DNAm:    %d probes, array guess %s\n",
      d[["array"]][["n_probes"]],
      d[["array"]][["guess"]]
    ))
    cat(sprintf(
      "  betas:   range [%.3g, %.3g], %d out-of-range probe(s)\n",
      d[["beta"]][["min"]],
      d[["beta"]][["max"]],
      d[["beta"]][["n_out_of_range"]]
    ))
    cov <- d[["coverage"]][["coverage"]]
    cat(sprintf(
      "  clocks:  %d assessed, %d under 75%% coverage\n",
      nrow(d[["coverage"]]),
      sum(cov < 0.75, na.rm = TRUE)
    ))
  }
  s <- x[["score"]]
  if (!is.null(s)) {
    cat(sprintf("  scores:  %d clock(s)", nrow(s[["summary"]])))
    if (!is.null(s[["flags"]])) {
      cat(sprintf(", %d flagged vs reference", nrow(s[["flags"]])))
    }
    cat("\n")
  }
  if (!is.null(m[["file"]])) {
    cat(sprintf("  file:    %s\n", m[["file"]]))
  }
  invisible(x)
}
