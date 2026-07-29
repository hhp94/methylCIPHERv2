# partial NA -> cohort mean, fully absent -> vendor ref

# value gates over one col_stats() sweep: overflow stops, everything else warns
check_col_values <- function(scan, cols) {
  at <- scan[["overflow_col"]]
  if (!is.null(at)) {
    cli::cli_abort(
      c(
        "DNAm column {.val {cols[at[[1L]]]}} does not sum to a finite value.",
        "i" = "Its entries are finite but astronomically large -- far outside
               the {.val {0}} to {.val {1}} beta range. Fix the matrix before
               scoring.",
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
        "i" = "An infinite beta is usually an upstream divide-by-zero.
               {.fn clocks_coverage} reports what was imputed."
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
               {.val {1}}. An M-value matrix is the usual cause -- it scores
               without error but is wrong.",
        "i" = "Convert with {.code beta <- 2^m / (2^m + 1)} if that is it."
      ),
      call = NULL
    )
  }
  if (isTRUE(scan[["any_gt1"]])) {
    cli::cli_warn(
      c(
        "DNAm contains values above {.val {1}}.",
        "i" = "{.fn calc_clocks} expects beta values in {.val {0}} to
               {.val {1}}."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# one col_stats() sweep: classify columns, means, value gates, row_obs
scan_missing_cpgs <- function(DNAm, needed_cpgs, score_cpgs) {
  present_needed <- intersect(needed_cpgs, colnames(DNAm))
  nr <- nrow(DNAm)

  # index, not a slice: the kernel strides over DNAm's own columns
  scan <- col_stats(DNAm, match(present_needed, colnames(DNAm)))
  check_col_values(scan, present_needed)

  # a dead row is dead on the scoring panels: a normalization CpG can never
  # score a sample. empty panel is the coverage gate's problem, not this one.
  present_score <- intersect(score_cpgs, colnames(DNAm))
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
          "i" = "Remove or fix {cli::qty(dead)}{?it/them} before scoring."
        ),
        call = NULL
      )
    }
  }

  # past the Inf gate, stats is populated
  st <- scan[["stats"]]
  n_obs <- st["n_obs", ]
  all_na <- present_needed[n_obs == 0]
  partial <- present_needed[n_obs > 0 & n_obs < nr]
  i <- match(partial, present_needed)

  # only partial columns get a mean (all-NA columns are classified, not divided)
  list(
    usable_cols = setdiff(present_needed, all_na),
    partial_na_cols = partial,
    all_na_cols = all_na,
    col_mean = stats::setNames(st["sum", i] / st["n_obs", i], partial)
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
