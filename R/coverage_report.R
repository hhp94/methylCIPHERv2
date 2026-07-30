# formatters over a finished record's $coverage (no re-touch of beta)

check_mc_result <- function(x, arg = "x") {
  if (!inherits(x, "mc_result")) {
    stop(
      sprintf("%s must be an mc_result from calc_clocks().", arg),
      call. = FALSE
    )
  }
  invisible(x)
}

# per-sample miss from the finished panel matrix. column is always there by
# construction, so absence is a bug (panel_ratio would silently return numeric(0))
miss_vec <- function(x, id, panel = c("score", "norm")) {
  panel <- match.arg(panel)
  m <- x[["coverage"]][["sample_miss"]][[panel]]
  if (is.null(m) || !(id %in% colnames(m))) {
    stop(
      sprintf(
        paste0(
          "%s has no %s-panel miss column. ",
          "This is a package bug -- please report it."
        ),
        id,
        panel
      ),
      call. = FALSE
    )
  }
  m[, id]
}

# one row per clock computed (aliases have NA panels)
#' @export
clocks_coverage <- function(x) {
  check_mc_result(x)
  per_clock <- x[["coverage"]][["per_clock"]]
  returned <- x[["provenance"]][["clocks"]]
  ids <- names(per_clock)

  int_field <- function(nm) {
    unname(vapply(
      per_clock,
      function(r) if (is.null(r)) NA_integer_ else as.integer(r[[nm]]),
      integer(1L)
    ))
  }

  out <- data.frame(
    clock_id = ids,
    # from the catalog, not the record -- a NULL record still has a group
    group_id = unname(vapply(ids, clock_group_id, character(1L))),
    role = ifelse(ids %in% returned, "returned", "routing_target"),
    policy = unname(vapply(
      per_clock,
      function(r) if (is.null(r)) NA_character_ else r[["policy"]],
      character(1L)
    )),
    normalizes = unname(vapply(
      per_clock,
      function(r) if (is.null(r)) NA else r[["normalizes"]],
      logical(1L)
    )),
    score_needed = int_field("score_needed"),
    score_present = int_field("score_present"),
    score_used = int_field("score_used"),
    score_imputed_partial = int_field("score_imputed_partial"),
    score_imputed_full = int_field("score_imputed_full"),
    score_dropped = int_field("score_dropped"),
    norm_needed = int_field("norm_needed"),
    norm_present = int_field("norm_present"),
    norm_imputed_partial = int_field("norm_imputed_partial"),
    norm_imputed_full = int_field("norm_imputed_full"),
    norm_dropped = int_field("norm_dropped"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # absent-probe list stays a list-column
  out[["missing_cpgs"]] <- unname(lapply(
    per_clock,
    function(r) if (is.null(r)) character(0) else r[["missing_cpgs"]]
  ))
  out
}

# one panel's per-sample rows for a non-alias returned clock
panel_rows <- function(id, panel, ratio, sample_id) {
  data.frame(
    id = sample_id,
    clock_id = id,
    panel = panel,
    n_observed = as.integer(ratio[["n_observed"]]),
    n_needed = as.integer(ratio[["needed"]]),
    coverage = ratio[["cov"]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# score-panel rows, plus a norm-panel row when the clock normalizes
clock_sample_rows <- function(x, id, sample_id) {
  rec <- x[["coverage"]][["per_clock"]][[id]]

  rows <- list()
  rows[["score"]] <- panel_rows(
    id,
    "score",
    panel_ratio(
      rec[["score_present"]],
      miss_vec(x, id, "score"),
      rec[["score_needed"]]
    ),
    sample_id
  )
  # norm panel only when the clock normalizes (same condition as having a
  # norm column), so the read sits inside the guard
  if (rec[["normalizes"]]) {
    rows[["norm"]] <- panel_rows(
      id,
      "norm",
      panel_ratio(
        rec[["norm_present"]],
        miss_vec(x, id, "norm"),
        rec[["norm_needed"]]
      ),
      sample_id
    )
  }
  do.call(rbind, rows)
}

# zero-row frame in samples_coverage() shape -- nothing to report is not an error
empty_sample_rows <- function() {
  data.frame(
    id = character(0),
    clock_id = character(0),
    panel = character(0),
    n_observed = integer(0),
    n_needed = integer(0),
    coverage = numeric(0),
    stringsAsFactors = FALSE
  )
}

# one row per (sample, returned clock, panel)
#' @export
samples_coverage <- function(x) {
  check_mc_result(x)
  sample_id <- x[["provenance"]][["sample_id"]]
  per_clock <- x[["coverage"]][["per_clock"]]

  # no record means no CpGs of its own. honest per-sample rows are its
  # descendants' -- routing targets appear here, pure composites do not
  ids <- names(per_clock)[!vapply(per_clock, is.null, logical(1L))]
  parts <- lapply(ids, clock_sample_rows, x = x, sample_id = sample_id)

  # seed with the empty frame so a run of pure composites keeps the shape
  out <- do.call(rbind, c(list(empty_sample_rows()), parts))

  # NA coverage is a routed member masked on a row its sex did not score.
  # drop those rows -- sample was never scored against that panel
  out <- out[!is.na(out[["coverage"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}
