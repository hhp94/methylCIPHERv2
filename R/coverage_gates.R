# pre-score coverage gates over the resolved panels

# requestable token for a compute-sequence id (alias, not routed member).
gate_label <- function(id, routed = sex_routed_members()) {
  if (!id %in% names(routed[["alias"]])) {
    return(cli::format_inline("{.val {id}}"))
  }
  cli::format_inline(
    "{.val {routed[['alias']][[id]]}} ({routed[['sex']][[id]]} model)"
  )
}

check_coverage <- function(cpg_list, threshold = 0.75) {
  # threshold is min_clocks_coverage, already validated at the front door
  # warn within 10% of the floor, before the gate itself trips
  warn_below <- min(1, threshold * 1.1)
  routed <- sex_routed_members()
  per_clock <- cpg_list[["per_clock"]]

  # interpolated labels. braces cannot become a cli template.
  panel_line <- function(id, present, needed, label) {
    cli::format_inline(
      "{gate_label(id, routed)}: {length(present)}/{length(needed)}
       {label} CpGs ({round(100 * length(present) / length(needed), 1)}%)"
    )
  }
  score_lines <- function(ids) {
    vapply(
      ids,
      function(id) {
        x <- per_clock[[id]]
        panel_line(id, x[["score_present"]], x[["score_needed"]], "scoring")
      },
      character(1L)
    )
  }
  norm_lines <- function(ids) {
    vapply(
      ids,
      function(id) {
        x <- per_clock[[id]]
        panel_line(id, x[["norm_present"]], x[["norm_needed"]], "normalization")
      },
      character(1L)
    )
  }

  classify <- function(x) {
    if (!length(x[["score_needed"]])) {
      return("")
    }
    ratio <- length(x[["score_present"]]) / length(x[["score_needed"]])
    # empty panel is scoreable only under vendor_mean fill
    undefined <- ratio == 0 &&
      !identical(clock_impute(x[["clock_id"]])[["policy"]], "vendor_mean")
    if (undefined || ratio < threshold) {
      "stop"
    } else if (ratio < warn_below) {
      "warn"
    } else {
      ""
    }
  }

  graded <- vapply(per_clock, classify, character(1L))
  ids_for <- function(lvl) names(graded)[graded == lvl]

  fail <- ids_for("stop")
  if (length(fail)) {
    cli::cli_abort(
      c(
        "{length(fail)} clock{?s} {?has/have} too few CpGs in {.arg DNAm} to
         score ({.arg min_clocks_coverage} = {format(threshold)}):",
        capped_bullets(fail, score_lines),
        "i" = "Remove {cli::qty(fail)}{?it/them} from {.arg clocks}, or lower
               {.arg min_clocks_coverage}.",
        "i" = "Call {.fn clock_cpgs} with a clock id to list every CpG that
               clock needs."
      ),
      call = NULL
    )
  }

  marginal <- ids_for("warn")
  if (length(marginal)) {
    cli::cli_warn(
      c(
        "{length(marginal)} clock{?s} {?is/are} just above
         {.arg min_clocks_coverage} = {format(threshold)}:",
        capped_bullets(marginal, score_lines),
        "i" = "Call {.fn clock_cpgs} with a clock id to list every CpG that
               clock needs.",
        "i" = "See {.fn clocks_coverage} for the panel counts per clock."
      ),
      call = NULL
    )
  }

  # thin QN backgrounds warn only
  thin <- names(per_clock)[vapply(
    per_clock,
    function(x) {
      length(x[["norm_needed"]]) > 0L &&
        length(x[["norm_present"]]) / length(x[["norm_needed"]]) < threshold
    },
    logical(1L)
  )]
  if (length(thin)) {
    # qn fills absent background CpGs from the target, BMIQ does not
    thin_schemes <- unique(vapply(thin, clock_norm_scheme, character(1)))
    fate <- if (all(thin_schemes == "bmiq")) {
      "The absent CpGs are dropped from the BMIQ fit."
    } else if (any(thin_schemes == "bmiq")) {
      c(
        "The absent CpGs are dropped from the BMIQ fit.",
        "For quantile normalization, the absent CpGs are filled from the
         reference mean."
      )
    } else {
      "The absent CpGs are filled from the reference mean."
    }
    cli::cli_warn(
      c(
        "{length(thin)} clock{?s} {?has/have} too few normalization CpGs
         (below {.arg min_clocks_coverage} = {format(threshold)}):",
        capped_bullets(thin, norm_lines),
        stats::setNames(fate, rep("i", length(fate))),
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

  # counts per clock first, strings only for the ids that survive the cap
  stat_for <- function(id) {
    rc <- row_coverage(
      coverage[["per_clock"]][[id]],
      coverage[["sample_miss"]][["score"]][[id]],
      coverage[["sample_miss"]][["norm"]][[id]]
    )
    if (is.null(rc)) {
      return(NULL)
    }
    cov <- rc[["cov"]]
    low <- !is.na(cov) & cov < threshold
    if (!any(low)) {
      return(NULL)
    }
    list(
      low = sum(low),
      scored = sum(!is.na(cov)),
      worst = min(cov[low]),
      needed = rc[["needed"]]
    )
  }

  ids <- names(coverage[["per_clock"]])
  hits <- stats::setNames(lapply(ids, stat_for), ids)
  hits <- hits[!vapply(hits, is.null, logical(1L))]

  row_lines <- function(these) {
    vapply(
      these,
      function(id) {
        s <- hits[[id]]
        cli::format_inline(
          "{gate_label(id, routed)}: {s[['low']]} of {s[['scored']]}
           sample{?s}, worst {round(100 * s[['worst']], 1)}% of
           {s[['needed']]} CpGs"
        )
      },
      character(1L)
    )
  }

  if (length(hits)) {
    cli::cli_warn(
      c(
        "{length(hits)} clock{?s} scored some samples below
         {.arg min_samples_coverage} = {format(threshold)}:",
        capped_bullets(names(hits), row_lines),
        "i" = "Call {.fn samples_coverage} to see the coverage of every
               sample.",
        "i" = "Drop those samples from {.arg DNAm}, or lower
               {.arg min_samples_coverage}."
      ),
      call = NULL
    )
  }

  invisible(NULL)
}
