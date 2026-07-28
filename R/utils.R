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
  invisible(NULL)
}

# Scoring-time failures, per clock. Coverage is computed before any score, so
# it cannot see a sample the scorer itself could not fit -- that sample gets an
# NA score with full coverage reported, indistinguishable on a saved record
# from any other NA. A branch still returns only its score: it writes the
# sample ids here instead, and score_cohort() harvests the block's collector
# once. Nothing routes on this and no score reads it.
new_notes <- function() {
  new.env(parent = emptyenv())
}

# record that `sample_id` could not be scored for clock `id`
note_scoring_failure <- function(block, id, sample_id) {
  notes <- block[["notes"]]
  if (!is.environment(notes) || !length(sample_id)) {
    return(invisible(NULL))
  }
  notes[[id]] <- union(notes[[id]], sample_id)
  invisible(NULL)
}

# collector -> plain named list, clock order stable; empty when nothing failed
collect_notes <- function(notes) {
  ids <- sort(names(notes))
  if (!length(ids)) {
    return(list())
  }
  stats::setNames(lapply(ids, function(id) notes[[id]]), ids)
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
observed_panel <- function(present, block) {
  cache <- block[["partial_cache"]]
  cached <- cached_cols(present, cache)
  raw <- setdiff(present, cached)
  list(
    cols = c(cached, raw),
    values = cbind(
      cache[, cached, drop = FALSE],
      block[["DNAm"]][, raw, drop = FALSE]
    )
  )
}

# Which rows each sex-routed member scored, as two masks over the n samples in
# hand. The score and every per-sample quantity derived from it read this same
# split, which is what makes the alias's stitched `sample_miss` true row for row.
#
# Unknown sex is "neither", and so is no sex column at all: the alias leaves
# those rows NA rather than routing them to a member, so a member must not be
# credited with having covered them either.
sex_rows <- function(female, n) {
  none <- rep(FALSE, n)
  if (!length(female)) {
    return(list(female = none, male = none))
  }
  f <- as.numeric(female)
  known <- !is.na(f)
  list(female = known & f == 1, male = known & f == 0)
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

# What a clock's fully-absent CpGs contribute, under the policy the catalog
# declares: `vendor_mean` folds each one in at its vendored value and counts it
# as a term, anything else drops it so it contributes nothing and counts for
# nothing. Every scoring site reads the policy through here, because
# coverage_record() reports score_used / score_imputed_full / score_dropped off
# the same declaration -- a sync that flips a policy has to move the score and
# the counts together, not one without the other.
#
# `ref` is resolved only on the fill path: a clock that drops has no vendored
# ref, and clock_impute_ref() stops rather than inventing one.
absent_fill <- function(id, coef, absent, ref = NULL, label = id) {
  no_fill <- list(offset = 0, filled = character(0))
  if (!length(absent)) {
    return(no_fill)
  }
  if (!identical(clock_impute(id)[["policy"]], "vendor_mean")) {
    return(no_fill)
  }
  if (is.null(ref)) {
    ref <- clock_impute_ref(id)
  }
  list(offset = vendor_offset(coef, absent, ref, label), filled = absent)
}

# one clock's coverage record (the only writer of record fields)
coverage_record <- function(cpgs, score_miss, norm_miss = NULL) {
  policy <- clock_impute(cpgs[["clock_id"]])[["policy"]]
  fill <- identical(policy, "vendor_mean")
  n_absent <- length(cpgs[["score_absent"]])
  list(
    clock_id = cpgs[["clock_id"]],
    policy = policy,
    normalizes = cpgs[["normalizes"]],
    score_needed = length(cpgs[["score_needed"]]),
    score_present = length(cpgs[["score_present"]]),
    score_used = length(cpgs[["score_present"]]) + if (fill) n_absent else 0L,
    score_imputed_partial = sum(score_miss, na.rm = TRUE),
    score_imputed_full = if (fill) n_absent else 0L,
    score_dropped = if (fill) 0L else n_absent,
    norm_needed = length(cpgs[["norm_needed"]]),
    norm_present = length(cpgs[["norm_present"]]),
    norm_imputed_partial = if (is.null(norm_miss)) {
      0L
    } else {
      sum(norm_miss, na.rm = TRUE)
    },
    missing_cpgs = cpgs[["score_absent"]]
  )
}

# y = sum_k coef[k+1] * x^k, lowest degree first. Accumulated rather than built
# as a power matrix: the matrix form collapses to a vector at length(x) == 1 and
# takes the score with it.
poly_eval <- function(x, coef) {
  out <- rep(0, length(x))
  for (k in seq_along(coef)) {
    out <- out + coef[[k]] * x^(k - 1L)
  }
  out
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
