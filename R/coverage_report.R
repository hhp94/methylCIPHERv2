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

# per-sample miss from the finished panel matrix. a missing column is a bug
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

# one batch's rows (aliases have NA panels)
batch_coverage <- function(per_clock, batch, returned) {
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
  # batch last: a hash next to clock_id reads as noise, but it is still the join key
  out[["batch"]] <- batch
  out
}

# one row per (clock, batch). counts are only true of the batch that produced them
#' @export
clocks_coverage <- function(x) {
  check_mc_result(x)
  batches <- x[["coverage"]][["per_clock"]]
  returned <- x[["provenance"]][["clocks"]]
  out <- do.call(
    rbind,
    lapply(names(batches), function(b) {
      batch_coverage(batches[[b]], b, returned)
    })
  )
  rownames(out) <- NULL
  out
}

# one panel's per-sample rows for a non-alias returned clock
panel_rows <- function(id, panel, batch, ratio, sample_id) {
  data.frame(
    id = sample_id,
    clock_id = id,
    panel = panel,
    n_observed = as.integer(ratio[["n_observed"]]),
    n_needed = as.integer(ratio[["needed"]]),
    coverage = ratio[["cov"]],
    # last, like clocks_coverage() -- the column the two frames join on
    batch = batch,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# score-panel rows, plus a norm row when the clock normalizes
clock_sample_rows <- function(x, id, rec, batch, rows) {
  sample_id <- x[["provenance"]][["sample_id"]][rows]

  out <- list()
  out[["score"]] <- panel_rows(
    id,
    "score",
    batch,
    panel_ratio(
      rec[["score_present"]],
      miss_vec(x, id, "score")[rows],
      rec[["score_needed"]]
    ),
    sample_id
  )
  # norm panel only when the clock normalizes
  if (rec[["normalizes"]]) {
    out[["norm"]] <- panel_rows(
      id,
      "norm",
      batch,
      panel_ratio(
        rec[["norm_present"]],
        miss_vec(x, id, "norm")[rows],
        rec[["norm_needed"]]
      ),
      sample_id
    )
  }
  do.call(rbind, out)
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
    batch = character(0),
    stringsAsFactors = FALSE
  )
}

# one row per (sample, returned clock, panel)
#' @export
samples_coverage <- function(x) {
  check_mc_result(x)
  batch <- x[["provenance"]][["batch"]]

  parts <- list()
  for (b in names(x[["coverage"]][["per_clock"]])) {
    per_clock <- x[["coverage"]][["per_clock"]][[b]]
    rows <- batch == b
    # no record means no cpgs of its own. read descendants for pure composites
    ids <- names(per_clock)[!vapply(per_clock, is.null, logical(1L))]
    parts <- c(
      parts,
      lapply(ids, function(id) {
        clock_sample_rows(x, id, per_clock[[id]], b, rows)
      })
    )
  }

  # seed with the empty frame so a run of pure composites keeps the shape
  out <- do.call(rbind, c(list(empty_sample_rows()), parts))

  # drop na coverage rows (routed member on a sex it did not score)
  out <- out[!is.na(out[["coverage"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}
