# partial NA -> cohort mean, fully absent -> vendor ref

# Value gates over one `col_stats()` sweep. The kernel counts, this decides:
# +/-Inf is bad data no fill can repair, so it stops; a value outside [0, 1] is
# usually an M-value matrix handed over as betas, which scores plausible
# garbage, so it warns rather than stops -- the range is a strong tell, not a
# definition.
#
# `inf_at` is checked first and nothing else in `scan` is read in that branch:
# the kernel bails at the first Inf, so `stats` is NULL and the range flags
# only cover the prefix it managed to scan.
check_col_values <- function(scan, DNAm, cols) {
  at <- scan[["inf_at"]]
  if (!is.null(at)) {
    cpg <- cols[at[[2L]]]
    # the kernel finds a bad column by its sum going non-finite, so it reports
    # a column always and a row only when an actual Inf sits in it. A missing
    # row means the column is finite but overflowed, which is a different
    # defect and needs a different sentence.
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

  # two flags, two warnings: a matrix can carry both, and neither implies the
  # other. Neither names a column -- the kernel keeps one global flag per
  # side rather than a per-column range.
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

# NA scan over needed columns. One `col_stats()` traversal yields the column
# classification, the means the fill uses, the value gates and the per-row
# observed count, so none of the four can disagree with the others. It runs
# unconditionally: gating it on `anyNA()` would skip the value gates on a
# matrix whose only defect is an Inf, since `anyNA()` does not see one.
#
# The row half (all-NA sample abort) is per row and stays correct against any
# subset of the cohort. It reads `row_obs`, which is counted over the panel
# slice, not the full array: a sample with nothing observed on any scoring
# panel is unscoreable however many off-panel CpGs it carries, and the old
# full-width test let exactly that case through.
scan_missing_cpgs <- function(DNAm, needed_cpgs) {
  present_needed <- intersect(needed_cpgs, colnames(DNAm))
  nr <- nrow(DNAm)

  scan <- col_stats(DNAm[, present_needed, drop = FALSE])
  check_col_values(scan, DNAm, present_needed)

  # no panel columns at all is the coverage gate's diagnostic, not this one's
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

  # past the Inf gate, so `stats` is populated
  st <- scan[["stats"]]
  n_obs <- st["n_obs", ]
  all_na <- present_needed[n_obs == 0]
  partial <- present_needed[n_obs > 0 & n_obs < nr]
  i <- match(partial, present_needed)

  # `n_obs == 0` is handled by classification, not by division: an all-missing
  # column is an ordinary, expected case -- it lands in `all_na_cols`, leaves
  # `usable_cols`, and takes the clock's vendored ref or is dropped by policy.
  # Only `partial` columns get a mean, and `n_obs > 0` holds for every one of
  # them by construction, so `sum / n_obs` is never 0/0.
  list(
    usable_cols = setdiff(present_needed, all_na),
    partial_na_cols = partial,
    all_na_cols = all_na,
    col_mean = stats::setNames(st["sum", i] / st["n_obs", i], partial)
  )
}

# apply cohort means to this matrix's NAs, one column per `partial_fill` entry.
# `sub` is a fresh slice and that is load-bearing: fill_imp_col() mutates its
# argument in place, so it must never be handed the caller's own DNAm.
build_partial_cache <- function(DNAm, partial_fill) {
  if (!length(partial_fill)) {
    return(NULL)
  }
  sub <- DNAm[, names(partial_fill), drop = FALSE]
  fill_imp_col(sub, partial_fill)
  sub
}
