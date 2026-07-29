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

# per-sample miss vector from the finished panel matrix, or NULL
miss_vec <- function(x, id, panel = c("score", "norm")) {
  panel <- match.arg(panel)
  m <- x[["coverage"]][["sample_miss"]][[panel]]
  if (!is.null(m) && id %in% colnames(m)) m[, id] else NULL
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
  sm <- miss_vec(x, id, "score")
  nm <- miss_vec(x, id, "norm")

  rows <- list()
  rows[["score"]] <- panel_rows(
    id,
    "score",
    panel_ratio(rec[["score_present"]], sm, rec[["score_needed"]]),
    sample_id
  )
  # norm panel only when the clock normalizes
  if (rec[["normalizes"]]) {
    rows[["norm"]] <- panel_rows(
      id,
      "norm",
      panel_ratio(rec[["norm_present"]], nm, rec[["norm_needed"]]),
      sample_id
    )
  }
  do.call(rbind, rows)
}

# alias sample rows: denominators from the member that scored each row
alias_sample_rows <- function(x, alias, sample_id) {
  route <- clock_routing(alias)
  pheno_id <- x[["provenance"]][["pheno_id"]]
  female <- rep(NA_integer_, length(sample_id))
  pheno <- x[["pheno"]]
  if (!is.null(pheno) && "Female" %in% names(pheno)) {
    idx <- match(sample_id, pheno[[pheno_id]])
    female <- as.integer(pheno[["Female"]])[idx]
  }
  miss <- miss_vec(x, alias, "score")

  n_observed <- rep(NA_integer_, length(sample_id))
  n_needed <- rep(NA_integer_, length(sample_id))
  coverage <- rep(NA_real_, length(sample_id))

  rows <- sex_rows(female, length(sample_id))
  for (sx in names(rows)) {
    member <- as.character(route[[sx]])
    rec <- x[["coverage"]][["per_clock"]][[member]]
    if (is.null(rec)) {
      next
    }
    applies <- rows[[sx]]
    if (!any(applies)) {
      next
    }
    member_miss <- rep(NA_integer_, length(sample_id))
    member_miss[applies] <- miss[applies]
    rc <- row_coverage(rec, member_miss, NULL)
    coverage[applies] <- rc[["cov"]][applies]
    n_needed[applies] <- rc[["needed"]]
    n_observed[applies] <- as.integer(rc[["n_observed"]][applies])
  }

  data.frame(
    id = sample_id,
    clock_id = alias,
    panel = "score",
    n_observed = n_observed,
    n_needed = n_needed,
    coverage = coverage,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# one row per (sample, returned clock, panel)
#' @export
samples_coverage <- function(x) {
  check_mc_result(x)
  sample_id <- x[["provenance"]][["sample_id"]]
  returned <- x[["provenance"]][["clocks"]]
  per_clock <- x[["coverage"]][["per_clock"]]

  parts <- lapply(returned, function(id) {
    if (is.null(per_clock[[id]])) {
      alias_sample_rows(x, id, sample_id)
    } else {
      clock_sample_rows(x, id, sample_id)
    }
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
