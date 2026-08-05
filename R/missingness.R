# partial NA -> cohort mean, fully absent -> vendor ref

# max this far past 1 is treated as percent methylation.
PERCENT_SCALE_AT <- 50

# value gates over one col_stats() sweep: overflow stops, everything else warns
check_col_values <- function(scan, cols) {
  at <- scan[["overflow_col"]]
  if (!is.null(at)) {
    cli::cli_abort(
      c(
        "{.arg DNAm} column {.val {cols[at[[1L]]]}} does not sum to a finite
         value.",
        "i" = "The entries are finite but very large. Beta values normally run
               from {.val {0}} to {.val {1}}.",
        "i" = "Use {.fn range} on that column before you score.",
        "i" = "The value check stops at the first such column, so the matrix
               may hold others."
      ),
      call = NULL
    )
  }

  # an Inf is missing, not fatal -- but it is a data bug worth naming
  if (isTRUE(scan[["any_inf"]])) {
    cli::cli_warn(
      c(
        "{.arg DNAm} contains infinite values.",
        "i" = "An infinite value counts as missing.",
        "i" = "The cohort mean fills it where the probe is partly observed.",
        "i" = "The probe counts as absent where no sample observes it.",
        "i" = "An infinite beta is often a divide by zero earlier in the
               pipeline.",
        "i" = "{.fn clocks_coverage} reports what the fill replaced."
      ),
      call = NULL
    )
  }

  # range is a running min/max seeded at beta bounds. only panel columns are scanned.
  lo <- scan[["min_val"]]
  if (lo < 0) {
    cli::cli_warn(
      c(
        "{.arg DNAm} contains values below {.val {0}}.",
        "x" = "The smallest is {.val {signif(lo, 4)}}, in column
               {.val {cols[scan[['min_col']]]}}.",
        "i" = "{.fn calc_clocks} expects beta values from {.val {0}} to
               {.val {1}}.",
        "i" = "An M-value matrix is a common cause. The resulting ages are not
               meaningful.",
        "i" = "Convert an M-value matrix with {.code beta <- 2^m / (2^m + 1)}."
      ),
      call = NULL
    )
  }
  hi <- scan[["max_val"]]
  if (hi > 1) {
    cli::cli_warn(
      c(
        "{.arg DNAm} contains values above {.val {1}}.",
        "x" = "The largest is {.val {signif(hi, 4)}}, in column
               {.val {cols[scan[['max_col']]]}}.",
        "i" = "{.fn calc_clocks} expects beta values from {.val {0}} to
               {.val {1}}.",
        "i" = "The resulting ages are not meaningful.",
        if (hi > PERCENT_SCALE_AT) {
          c(
            "i" = "Percent methylation is the usual cause at this size.",
            "i" = "Convert percent methylation with {.code DNAm / 100}."
          )
        } else {
          c("i" = "Check the scale of {.arg DNAm}.")
        }
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# post-score value gate for nan/inf. na is legitimate.
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

  bad_lines <- function(ids) {
    vapply(
      ids,
      function(id) {
        cli::format_inline(
          "{.val {id}}: {bad[[id]]} of {length(scores[[id]])} sample{?s}"
        )
      },
      character(1L)
    )
  }
  # sample_scale clocks divide by per-sample sd (may be 0 or undefined).
  scaled <- names(bad)[vapply(
    names(bad),
    function(id) !is.null(clock_moment_key(id)),
    logical(1)
  )]
  full <- scaled[vapply(scaled, clock_needs_full_panel, logical(1))]
  ref <- setdiff(scaled, full)
  hint <- c(
    if (length(full)) {
      c(
        "i" = "{.val {full}} divide{cli::qty(full)}{?s/} by a per-sample
               standard deviation over every column of {.arg DNAm}.",
        "i" = "A sample with one distinct value has no spread to divide by."
      )
    },
    if (length(ref)) {
      c(
        "i" = "{.val {ref}} divide{cli::qty(ref)}{?s/} by a per-sample standard
               deviation over {cli::qty(ref)}{?its/their} declared reference
               set.",
        "i" = "A sample with one distinct value on that set has no spread to
               divide by."
      )
    }
  )

  cli::cli_warn(
    c(
      "{length(bad)} clock{?s} produced {cli::qty(sum(bad))}non-finite
       score{?s}:",
      capped_bullets(names(bad), bad_lines),
      hint,
      "i" = "A {.code NaN} or an {.code Inf} usually means a non-finite value
             reached the score calculation.",
      "i" = "Check {.arg DNAm} rather than the score."
    ),
    call = NULL
  )
  invisible(NULL)
}

# max moment domains per col_stats() sweep (uint8_t mask width).
MAX_MOMENT_SETS <- 8L

# moment_sets for col_stats(): NULL, or a list of column-index vectors.
check_moment_sets <- function(sets, nc) {
  if (is.null(sets)) {
    return(NULL)
  }
  bug <- function(...) {
    stop(sprintf(...), call. = FALSE)
  }
  if (!is.list(sets) || !length(sets) || length(sets) > MAX_MOMENT_SETS) {
    bug(
      "moment_sets must be a list of 1 to %d index vectors, got %s of length %d.",
      MAX_MOMENT_SETS,
      class(sets)[[1L]],
      length(sets)
    )
  }
  labels <- names(sets) %||% rep("", length(sets))
  out <- lapply(seq_along(sets), function(k) {
    who <- if (nzchar(labels[[k]])) labels[[k]] else as.character(k)
    v <- sets[[k]]
    if (!is.numeric(v)) {
      bug("moment_sets[[%s]] must be numeric, got %s.", who, class(v)[[1L]])
    }
    if (anyNA(v)) {
      bug("moment_sets[[%s]] has missing values.", who)
    }
    if (any(v != trunc(v))) {
      bug("moment_sets[[%s]] has non-integer values.", who)
    }
    # 1-based column indices into DNAm -- the kernel does not range-check them
    if (length(v) && (min(v) < 1 || max(v) > nc)) {
      bug(
        "moment_sets[[%s]] indexes outside 1:%d (range %g to %g).",
        who,
        nc,
        min(v),
        max(v)
      )
    }
    as.integer(v)
  })
  names(out) <- names(sets)
  out
}

# domain cpgs -> column indices. NULL element is the whole matrix. declared refs keep only measured cols.
resolve_moment_sets <- function(domains, cpgs) {
  if (!length(domains)) {
    return(NULL)
  }
  lapply(domains, function(d) {
    if (is.null(d)) {
      return(seq_along(cpgs))
    }
    # one match() only. kernel ORs repeated indices, so no dedup needed.
    m <- match(d, cpgs)
    m[!is.na(m)]
  })
}

# per-output mean/sd. mean needs n >= 1, sd needs n >= 2.
split_moments <- function(scan, sets) {
  if (is.null(sets)) {
    return(NULL)
  }
  n_mom <- scan[["row_moment_obs"]]
  row_mean <- scan[["row_mean"]]
  row_m2 <- scan[["row_m2"]]
  # explicit [, k]: the kernel returns a matrix per output, one column per set
  out <- lapply(seq_along(sets), function(k) {
    nk <- n_mom[, k, drop = TRUE]
    mk <- row_mean[, k, drop = TRUE]
    sk <- sqrt(row_m2[, k, drop = TRUE] / (nk - 1))
    mk[nk < 1L] <- NA_real_
    sk[nk < 2L] <- NA_real_
    list(mean = mk, sd = sk)
  })
  names(out) <- names(sets)
  out
}

# one col_stats() sweep: columns, means, value gates, row_obs, moment domains.
scan_missing_cpgs <- function(
  DNAm,
  needed_cpgs,
  score_cpgs,
  moment_domains = NULL
) {
  present_needed <- intersect(needed_cpgs, colnames(DNAm))
  # unique by construction, or a repeated index double-counts the column stats
  needed_idx <- match(present_needed, colnames(DNAm))
  nr <- nrow(DNAm)

  # domains index DNAm directly and are validated before the kernel sees them
  sets <- check_moment_sets(
    resolve_moment_sets(moment_domains, colnames(DNAm)),
    ncol(DNAm)
  )

  # index into dnam, not a slice.
  scan <- col_stats(DNAm, needed_idx, sets)
  check_col_values(scan, present_needed)

  # dead samples checked on scoring panels only.
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
           scoring panel: {.val {capped_vals(dead)}}.",
          "i" = "Remove or repair {cli::qty(dead)}{?it/them} before you
                 score.",
          "i" = "{.code rowSums(!is.na(DNAm))} counts the observed values per
                 sample."
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

  # the kernel counts each domain itself, so the divisor is that domain's count
  moments <- split_moments(scan, sets)

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
