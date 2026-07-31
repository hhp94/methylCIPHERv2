# partial NA -> cohort mean, fully absent -> vendor ref

# value gates over one col_stats() sweep: overflow stops, everything else warns
check_col_values <- function(scan, cols) {
  at <- scan[["overflow_col"]]
  if (!is.null(at)) {
    cli::cli_abort(
      c(
        "DNAm column {.val {cols[at[[1L]]]}} does not sum to a finite value.",
        "i" = "Its entries look finite but are extremely large -- far outside
               the usual beta range of {.val {0}} to {.val {1}}. Please check
               that column before scoring.",
        "i" = "The scan stops at the first such column, so there may be others."
      ),
      call = NULL
    )
  }

  # an Inf is missing, not fatal -- but it is a data bug worth naming
  if (isTRUE(scan[["any_inf"]])) {
    cli::cli_warn(
      c(
        "DNAm contains infinite values.",
        "i" = "They are treated as missing: filled from the cohort mean where
               the probe is partly observed, and counted absent where it is
               not.",
        "i" = "An infinite beta is often an upstream divide-by-zero.
               {.fn clocks_coverage} shows what was imputed."
      ),
      call = NULL
    )
  }

  # range flags are global (a matrix can carry both sides)
  if (isTRUE(scan[["any_lt0"]])) {
    cli::cli_warn(
      c(
        "DNAm contains values below {.val {0}}.",
        "i" = "{.fn calc_clocks} expects beta values in {.val {0}} to
               {.val {1}}. An M-value matrix is a common cause -- it will
               score without error, but the ages won't be meaningful.",
        "i" = "If that sounds right, convert with
               {.code beta <- 2^m / (2^m + 1)}."
      ),
      call = NULL
    )
  }
  if (isTRUE(scan[["any_gt1"]])) {
    cli::cli_warn(
      c(
        "DNAm contains values above {.val {1}}.",
        "i" = "{.fn calc_clocks} expects beta values in {.val {0}} to
               {.val {1}}. You may want to double-check the scale of
               {.arg DNAm}."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# post-score value gate for nan/inf. na is legitimate (branch declined a sample)
check_score_values <- function(scores) {
  n_bad <- vapply(
    scores,
    function(v) sum(is.nan(v) | is.infinite(v)),
    integer(1L)
  )
  bad <- n_bad[n_bad > 0L]
  if (!length(bad)) {
    return(invisible(NULL))
  }

  lines <- sprintf(
    "%s: %d of %d sample(s)",
    names(bad),
    bad,
    lengths(scores[names(bad)])
  )
  # a full-panel clock divides by a per-sample sd, which can be 0 or undefined
  full <- names(bad)[vapply(names(bad), clock_needs_full_panel, logical(1))]
  hint <- if (length(full)) {
    c(
      "i" = "{.val {full}} divide{cli::qty(full)}{?s/} by a per-sample sd taken
             over every column of {.arg DNAm}, so a sample observing one value,
             or the same value everywhere, has no spread to scale by."
    )
  }

  cli::cli_warn(
    c(
      "{length(bad)} clock{?s} produced {cli::qty(sum(bad))}non-finite
       score{?s}:",
      capped_bullets(lines),
      hint,
      "i" = "{.code NaN} or {.code Inf} usually means a non-finite value
             reached the arithmetic. Please check {.arg DNAm} rather than
             the score itself."
    ),
    call = NULL
  )
  invisible(NULL)
}

# one col_stats() sweep for columns, means, value gates, row_obs, and optional moments
scan_missing_cpgs <- function(
  DNAm,
  needed_cpgs,
  score_cpgs,
  row_moments = FALSE
) {
  present_needed <- intersect(needed_cpgs, colnames(DNAm))
  # unique by construction. row_moments needs that or a repeated index double-counts
  needed_idx <- match(present_needed, colnames(DNAm))
  nr <- nrow(DNAm)

  # index into dnam, not a slice. with row_moments the complement only feeds row accumulators
  scan <- col_stats(DNAm, needed_idx, row_moments = row_moments)
  check_col_values(scan, present_needed)

  # dead samples checked on scoring panels only. identical inputs skip the rescan
  present_score <- if (identical(score_cpgs, needed_cpgs)) {
    present_needed
  } else {
    intersect(score_cpgs, colnames(DNAm))
  }
  if (length(present_score)) {
    obs <- if (identical(present_score, present_needed)) {
      scan[["row_obs"]]
    } else {
      row_observed(DNAm, match(present_score, colnames(DNAm)))
    }
    dead <- rownames(DNAm)[obs == 0L]
    if (length(dead)) {
      cli::cli_abort(
        c(
          "{length(dead)} sample{?s} {?has/have} no observed CpGs on any
           scoring panel: {.val {utils::head(dead, 10L)}}.",
          "i" = "Please remove or repair {cli::qty(dead)}{?it/them} before
                 scoring."
        ),
        call = NULL
      )
    }
  }

  # past the overflow gate, stats is populated -- and it is the panel's alone
  st <- scan[["stats"]]
  n_obs <- st["n_obs", ]
  all_na <- present_needed[n_obs == 0]
  partial <- present_needed[n_obs > 0 & n_obs < nr]
  i <- match(partial, present_needed)

  # moments span every column. sd's divisor is panel row_obs plus the complement
  moments <- if (row_moments) {
    n_all <- scan[["row_obs"]] + scan[["row_obs_complement"]]
    list(
      mean = scan[["row_mean"]],
      sd = sqrt(scan[["row_m2"]] / (n_all - 1))
    )
  }

  # only partial columns get a mean (all-NA columns are classified, not divided)
  list(
    usable_cols = setdiff(present_needed, all_na),
    partial_na_cols = partial,
    all_na_cols = all_na,
    col_mean = stats::setNames(st["sum", i] / st["n_obs", i], partial),
    sample_moments = moments
  )
}

# cohort-mean fill on a fresh slice (fill_imp_col mutates in place)
build_partial_cache <- function(DNAm, cols, partial_fill) {
  if (!length(partial_fill)) {
    return(NULL)
  }
  sub <- DNAm[, cols, drop = FALSE]
  fill_imp_col(sub, partial_fill)
  sub
}
