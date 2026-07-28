# partial NA -> cohort mean, fully absent -> vendor ref

# value gates over one col_stats() sweep: Inf stops, out-of-[0,1] warns
check_col_values <- function(scan, DNAm, cols) {
  at <- scan[["inf_at"]]
  if (!is.null(at)) {
    cpg <- cols[at[[2L]]]
    # row NA means overflow (finite entries, non-finite sum), not Inf
    if (is.na(at[[1L]])) {
      cli::cli_abort(
        c(
          "DNAm column {.val {cpg}} does not sum to a finite value.",
          "i" = "Its entries are finite but astronomically large -- far outside
                 the {.val {0}} to {.val {1}} beta range. Fix the matrix before
                 scoring.",
          "i" = "The scan stops at the first such column, so there may be
                 others."
        ),
        call = NULL
      )
    }
    cli::cli_abort(
      c(
        "DNAm contains an infinite value: sample {.val {rownames(DNAm)[at[[1L]]]}},
         CpG {.val {cpg}}.",
        "i" = "Infinite values are not missing values -- no imputation can
               repair one. Fix the matrix before scoring.",
        "i" = "The scan stops at the first one, so there may be others."
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
scan_missing_cpgs <- function(DNAm, needed_cpgs) {
  present_needed <- intersect(needed_cpgs, colnames(DNAm))
  nr <- nrow(DNAm)

  scan <- col_stats(DNAm[, present_needed, drop = FALSE])
  check_col_values(scan, DNAm, present_needed)

  # empty panel is the coverage gate's problem, not the all-NA row check
  if (length(present_needed)) {
    dead <- rownames(DNAm)[scan[["row_obs"]] == 0L]
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
build_partial_cache <- function(DNAm, partial_fill) {
  if (!length(partial_fill)) {
    return(NULL)
  }
  sub <- DNAm[, names(partial_fill), drop = FALSE]
  fill_imp_col(sub, partial_fill)
  sub
}
