# shared scoring helpers

# name a character vector as cli "*" bullets
bullets <- function(x) {
  stats::setNames(x, rep("*", length(x)))
}

# betanorm is a soft dep; every normalizing branch needs it
require_betanorm <- function(id) {
  if (!requireNamespace("betanorm", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "{.val {id}} needs the {.pkg betanorm} package for normalization.",
        "i" = "Install it from GitHub: {.code pak::pak(\"hhp94/betanorm\")}."
      ),
      call = NULL
    )
  }
  invisible(TRUE)
}

# present CpGs covered by the cohort-mean cache
cached_cols <- function(present, partial_cache) {
  if (is.null(partial_cache)) {
    character(0)
  } else {
    intersect(present, colnames(partial_cache))
  }
}

# observed betas for `present`: cohort-mean-filled columns first, then raw
observed_panel <- function(present, DNAm, partial_cache = NULL) {
  cached <- cached_cols(present, partial_cache)
  raw <- setdiff(present, cached)
  list(
    cols = c(cached, raw),
    values = cbind(
      partial_cache[, cached, drop = FALSE],
      DNAm[, raw, drop = FALSE]
    )
  )
}

# per-sample count of cohort-mean fills (always integer)
count_sample_miss <- function(DNAm, cached) {
  out <- if (length(cached)) {
    as.integer(slideimp::mat_miss(DNAm[, cached, drop = FALSE], col = FALSE))
  } else {
    integer(nrow(DNAm))
  }
  names(out) <- rownames(DNAm)
  out
}

# one panel's per-sample coverage from present/needed scalars and miss counts
panel_ratio <- function(present, miss, needed) {
  n_observed <- present - miss
  list(n_observed = n_observed, cov = n_observed / needed, needed = needed)
}

# vendor-mean fill for fully absent CpGs
vendor_offset <- function(coef, absent, ref, id) {
  miss_ref <- setdiff(absent, names(ref))
  if (length(miss_ref)) {
    cli::cli_abort(
      c(
        "{.val {id}}: no vendor mean for {length(miss_ref)} absent CpG{?s}.",
        "x" = "{.val {utils::head(miss_ref, 5L)}}"
      ),
      call = NULL
    )
  }
  sum(coef[absent] * ref[absent])
}

# one clock's coverage record (the only writer of record fields)
coverage_record <- function(cpgs, score_miss, norm_miss = NULL) {
  policy <- clock_impute(cpgs$clock_id)[["policy"]]
  fill <- identical(policy, "vendor_mean")
  n_absent <- length(cpgs$score_absent)
  list(
    clock_id = cpgs$clock_id,
    policy = policy,
    normalizes = isTRUE(cpgs$normalizes),
    score_needed = length(cpgs$score_needed),
    score_present = length(cpgs$score_present),
    score_used = length(cpgs$score_present) + if (fill) n_absent else 0L,
    score_imputed_partial = sum(score_miss, na.rm = TRUE),
    score_imputed_full = if (fill) n_absent else 0L,
    score_dropped = if (fill) 0L else n_absent,
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    norm_imputed_partial = if (is.null(norm_miss)) {
      0L
    } else {
      sum(norm_miss, na.rm = TRUE)
    },
    missing_cpgs = cpgs$score_absent
  )
}

# n x 1 score matrix every branch returns
score_matrix <- function(values, sample_id, id) {
  matrix(
    as.numeric(values),
    nrow = length(sample_id),
    ncol = 1L,
    dimnames = list(sample_id, id)
  )
}
