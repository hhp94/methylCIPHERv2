# pre-score coverage gates over the resolved panels

WARN_COVERAGE_MARGIN <- 1.1

# requestable token for a compute-sequence id (alias, not routed member)
gate_label <- function(id, routed = sex_routed_members()) {
  if (!id %in% names(routed[["alias"]])) {
    return(id)
  }
  sprintf("%s (%s model)", routed[["alias"]][[id]], routed[["sex"]][[id]])
}

check_coverage <- function(cpg_list, threshold = 0.75) {
  checkmate::assert_number(threshold, lower = 0, upper = 1)
  warn_below <- min(1, threshold * WARN_COVERAGE_MARGIN)
  routed <- sex_routed_members()

  panel_line <- function(id, present, needed, label) {
    sprintf(
      "%s: %d/%d %s CpGs (%.1f%%)",
      gate_label(id, routed),
      length(present),
      length(needed),
      label,
      100 * length(present) / length(needed)
    )
  }

  classify <- function(x) {
    if (!length(x[["score_needed"]])) {
      return(list(level = "", line = NA_character_))
    }
    ratio <- length(x[["score_present"]]) / length(x[["score_needed"]])
    # empty panel is scoreable only under vendor_mean fill
    undefined <- ratio == 0 &&
      !identical(clock_impute(x[["clock_id"]])[["policy"]], "vendor_mean")
    level <- if (undefined || ratio < threshold) {
      "stop"
    } else if (ratio < warn_below) {
      "warn"
    } else {
      ""
    }
    list(
      level = level,
      line = panel_line(
        x[["clock_id"]],
        x[["score_present"]],
        x[["score_needed"]],
        "scoring"
      )
    )
  }

  graded <- lapply(cpg_list[["per_clock"]], classify)
  levels <- vapply(graded, function(g) g[["level"]], character(1L))
  lines_for <- function(lvl) {
    vapply(graded[levels == lvl], function(g) g[["line"]], character(1L))
  }

  fail <- lines_for("stop")
  if (length(fail)) {
    cli::cli_abort(
      c(
        "{length(fail)} clock{?s} {?doesn't/don't} have enough CpGs to score
         ({.arg min_clocks_coverage} = {format(threshold)}):",
        capped_bullets(fail),
        "i" = "Try dropping {cli::qty(fail)}{?it/them} from {.arg clocks}, or
               lower {.arg min_clocks_coverage} if you meant to allow thinner
               panels."
      ),
      call = NULL
    )
  }

  marginal <- lines_for("warn")
  if (length(marginal)) {
    cli::cli_warn(
      c(
        "{length(marginal)} clock{?s} only just clear{?s/}
         {.arg min_clocks_coverage} = {format(threshold)}:",
        capped_bullets(marginal),
        "i" = "Scoring continues, but more of the panel will be filled by
               imputation."
      ),
      call = NULL
    )
  }

  # thin QN backgrounds warn only
  thin <- vapply(
    cpg_list[["per_clock"]],
    function(x) {
      if (
        !length(x[["norm_needed"]]) ||
          length(x[["norm_present"]]) / length(x[["norm_needed"]]) >= threshold
      ) {
        return(NA_character_)
      }
      panel_line(
        x[["clock_id"]],
        x[["norm_present"]],
        x[["norm_needed"]],
        "normalization"
      )
    },
    character(1L)
  )
  thin <- thin[!is.na(thin)]
  if (length(thin)) {
    # QN fills absent background CpGs from the target, BMIQ does not
    thin_schemes <- unique(vapply(names(thin), clock_norm_scheme, character(1)))
    fate <- if (all(thin_schemes == "bmiq")) {
      "Absent background CpGs are dropped from the calibration fit."
    } else if (any(thin_schemes == "bmiq")) {
      "Absent background CpGs are dropped from a BMIQ fit, and filled from
       the reference mean for quantile normalization."
    } else {
      "Missing background CpGs are filled from the reference mean."
    }
    cli::cli_warn(
      c(
        "{length(thin)} clock{?s} {?has/have} a thin normalization background
         (under {.arg min_clocks_coverage} = {format(threshold)}):",
        capped_bullets(thin),
        "i" = fate,
        "i" = "See {.fn clocks_coverage} for the panel counts per clock."
      ),
      call = NULL
    )
  }

  invisible(NULL)
}

# per-sample observed fraction of the row-gate panel (norm if normalizes, else score)
row_coverage <- function(cov, score_miss, norm_miss) {
  if (is.null(cov)) {
    return(NULL)
  }
  qn <- cov[["normalizes"]]
  # needed is a scalar count, so 0 is the only empty
  needed <- if (qn) cov[["norm_needed"]] else cov[["score_needed"]]
  present <- if (qn) cov[["norm_present"]] else cov[["score_present"]]
  miss <- if (qn) norm_miss else score_miss
  if (is.null(miss) || needed == 0L) {
    return(NULL)
  }
  panel_ratio(present, miss, needed)
}

# per-sample coverage gate (warn only) over the hoisted coverage structure
check_row_coverage <- function(coverage, threshold = 0.75) {
  # threshold is validated at the calc_clocks() front door, before scoring
  routed <- sex_routed_members()

  line_for <- function(id) {
    rc <- row_coverage(
      coverage[["per_clock"]][[id]],
      coverage[["sample_miss"]][["score"]][[id]],
      coverage[["sample_miss"]][["norm"]][[id]]
    )
    if (is.null(rc)) {
      return(NA_character_)
    }
    cov <- rc[["cov"]]
    low <- !is.na(cov) & cov < threshold
    if (!any(low)) {
      return(NA_character_)
    }
    sprintf(
      "%s: %d of %d sample(s), worst %.1f%% of %d CpGs",
      gate_label(id, routed),
      sum(low),
      sum(!is.na(cov)),
      100 * min(cov[low]),
      rc[["needed"]]
    )
  }

  lines <- vapply(names(coverage[["per_clock"]]), line_for, character(1L))
  lines <- lines[!is.na(lines)]
  if (length(lines)) {
    cli::cli_warn(
      c(
        "{length(lines)} clock{?s} scored some samples under
         {.arg min_samples_coverage} = {format(threshold)}:",
        capped_bullets(lines),
        "i" = "Those sample scores rely more on imputed CpGs, so interpret
               them with a bit of care."
      ),
      call = NULL
    )
  }

  invisible(NULL)
}
