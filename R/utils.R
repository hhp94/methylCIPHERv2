# shared scoring helpers

# name a character vector as cli "*" bullets
bullets <- function(x) {
  stats::setNames(x, rep("*", length(x)))
}

# betanorm is a soft dep (every normalizing branch needs it)
require_betanorm <- function(id) {
  if (!requireNamespace("betanorm", quietly = TRUE)) {
    stop(
      sprintf(
        paste0(
          "%s needs the betanorm package for normalization. ",
          "Install it from GitHub: pak::pak(\"hhp94/betanorm\")."
        ),
        id
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

# scoring-time failures per clock (coverage cannot see these)
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

# collector -> plain named list (empty when nothing failed)
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

# male/female masks for sex-routed members (unknown sex is neither)
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
    stop(
      sprintf(
        "%s: no vendor mean for %d absent CpG(s): %s.",
        id,
        length(miss_ref),
        paste(utils::head(miss_ref, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sum(coef[absent] * ref[absent])
}

# absent-CpG contribution under the declared policy (vendor_mean or drop)
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

# polynomial eval, lowest degree first (horner-style, 1-row safe)
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
