# Advisory age-association check: a cohort's observed score-age correlation
# against the shipped per-clock reference (DECISIONS 2026-08-02).

MIN_ASSOC_N <- 5L

# below this reference |r| the expected sign carries no information
SIGN_FLAG_MIN_R <- 0.3

ASSOC_COLS <- c(
  "clock_id",
  "n",
  "obs_age_r",
  "exp_age_r",
  "exp_lo",
  "exp_hi",
  "outside",
  "wrong_sign"
)

mc_clock_reference <- function() {
  path <- system.file(
    "extdata",
    "clock_reference.csv",
    package = "methylCIPHERv2"
  )
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

assoc_age <- function(x, age) {
  n <- nrow(x[["scores"]])
  if (is.null(age)) {
    pheno <- x[["pheno"]]
    if (!("Age" %in% names(pheno))) {
      cli::cli_abort(
        c(
          "{.arg x} has no {.field Age} column, and {.arg age} is
           {.code NULL}.",
          "i" = "Pass a numeric vector to {.arg age}.",
          "i" = "Or score with a {.arg pheno} that has an {.field Age}
                 column."
        ),
        call = NULL
      )
    }
    age <- pheno[["Age"]]
  }
  age <- suppressWarnings(as.numeric(age))
  if (length(age) != n) {
    cli::cli_abort(
      "{.arg age} has {length(age)} value{?s}. {.arg x} has {n} sample{?s}.",
      call = NULL
    )
  }
  age
}

assoc_row <- function(id, v, age, r) {
  ok <- is.finite(v) & is.finite(age)
  if (sum(ok) < MIN_ASSOC_N) {
    return(NULL)
  }
  if (stats::sd(v[ok]) == 0 || stats::sd(age[ok]) == 0) {
    return(NULL)
  }
  obs <- stats::cor(v[ok], age[ok])
  lo <- r[["age_r_lo"]]
  hi <- r[["age_r_hi"]]
  exp_r <- r[["age_r"]]
  data.frame(
    clock_id = id,
    n = sum(ok),
    obs_age_r = round(obs, 3),
    exp_age_r = exp_r,
    exp_lo = lo,
    exp_hi = hi,
    outside = isTRUE(is.finite(lo) && is.finite(hi) && (obs < lo || obs > hi)),
    wrong_sign = isTRUE(
      is.finite(exp_r) &&
        abs(exp_r) > SIGN_FLAG_MIN_R &&
        sign(obs) != sign(exp_r)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

assoc_empty <- function() {
  out <- data.frame(
    clock_id = character(0),
    n = integer(0),
    obs_age_r = numeric(0),
    exp_age_r = numeric(0),
    exp_lo = numeric(0),
    exp_hi = numeric(0),
    outside = logical(0),
    wrong_sign = logical(0),
    stringsAsFactors = FALSE
  )
  out[ASSOC_COLS]
}

#' Clock Age Associations
#'
#' Compares each clock's observed correlation with age against a reference
#' range shipped with the package.
#'
#' @inheritParams mc-params
#' @param age A numeric vector. The age for each sample, in the same order as
#'   the samples in `x`. Default is `NULL`, which uses the `Age` column of the
#'   `pheno` in `x`.
#'
#' @details
#' Each row compares a clock's observed score-age correlation with its
#' reference correlation and expected range. The `outside` column marks a
#' clock whose observed correlation falls outside that range. The
#' `wrong_sign` column marks a clock whose observed correlation has the
#' opposite sign from a reference correlation stronger than 0.3.
#'
#' A clock needs at least 5 samples with a finite score and age, and
#' variation in both, to appear in the result. A clock with no entry in the
#' shipped reference table is left out.
#'
#' This function recalculates any clock that depends on sample-wise
#' information, such as a z-score, from all the available samples when `x`
#' holds more than one batch. This is the same calculation as
#' [refinalize_clocks()].
#'
#' Rows are ordered by the gap between the observed and the reference
#' correlation, most negative first.
#'
#' @returns A data.frame. One row for each clock with a reference entry and
#'   enough samples. Columns are the clock id, the sample count, the
#'   observed and reference age correlations, the reference range, and the
#'   two flags described above.
#'
#' @seealso
#' [calc_accel()] for the age acceleration of each sample.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 20)
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#' score_associations(res, age = runif(20, 20, 80))
#'
#' @export
score_associations <- function(x, age = NULL) {
  # the reference table is package data, not a swappable argument
  assoc_report(x, age, mc_clock_reference())
}

# the body, with the reference injectable so a test can control it
assoc_report <- function(x, age, ref) {
  check_mc_result(x)
  if (is.null(ref)) {
    stop(
      "No clock reference table installed with methylCIPHERv2.",
      call. = FALSE
    )
  }
  x <- finalized(x)
  scores <- x[["scores"]]
  age <- assoc_age(x, age)

  ids <- intersect(colnames(scores), ref[["clock"]])
  rows <- lapply(ids, function(id) {
    assoc_row(id, scores[, id], age, ref[match(id, ref[["clock"]]), ])
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) {
    return(assoc_empty())
  }
  # widest gap between observed and expected first
  out <- out[order(out[["obs_age_r"]] - out[["exp_age_r"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}
