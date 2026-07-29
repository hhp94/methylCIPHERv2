# coverage/QC, once upstream of scoring, keyed by clock id

# per-sample observed count over column positions (callers are past the value gate)
row_observed <- function(DNAm, idx) {
  obs <- col_stats(DNAm, idx)[["row_obs"]]
  # NULL means the kernel bailed on a column whose sum overflowed, which the
  # value gate in scan_missing_cpgs() already ruled out for these columns
  if (is.null(obs)) {
    stop(
      "row_observed: col_stats bailed on a non-finite value past the value
       gate. This is a package bug -- please report it.",
      call. = FALSE
    )
  }
  obs
}

# per-sample count of cohort-mean fills over cached column positions (integer)
count_sample_miss <- function(DNAm, cached_idx) {
  out <- if (length(cached_idx)) {
    length(cached_idx) - row_observed(DNAm, cached_idx)
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

# per-sample partial-fill count over one panel's present CpGs
panel_sample_miss <- function(present, block) {
  # counted on the raw matrix -- the cache has the NAs already filled
  count_sample_miss(
    block[["DNAm"]],
    block_cols(cached_cols(present, block[["partial_cache"]]), block)
  )
}

# alias score-panel miss: stitch each row from the member that scored it
stitch_routed_sample_miss <- function(alias, score_miss, female, sample_id) {
  miss <- rep(NA_integer_, length(sample_id))
  names(miss) <- sample_id
  route <- clock_routing(alias)
  rows <- sex_rows(female, length(sample_id))
  for (key in names(rows)) {
    i <- rows[[key]]
    if (!any(i)) {
      next
    }
    sm <- score_miss[[as.character(route[[key]])]]
    if (!is.null(sm)) {
      miss[i] <- as.integer(sm)[i]
    }
  }
  miss
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

# full coverage for the compute sequence (per_clock records + per-panel sample_miss)
compute_coverage <- function(clock_sequence, cpg_list, block) {
  sample_id <- block[["sample_id"]]
  routed <- sex_routed_members()
  is_alias <- vapply(
    clock_sequence,
    function(p) identical(clock_kind(p), "sex_routed_alias"),
    logical(1L)
  )
  seqi <- seq_along(clock_sequence)
  pidx <- cpg_list[["panel_index"]]

  # count each distinct panel once, then fan out via the index
  score_part_miss <- lapply(
    pidx[["score"]][["parts"]],
    function(p) panel_sample_miss(p[["present"]], block)
  )
  norm_part_miss <- lapply(pidx[["norm"]][["parts"]], function(p) {
    if (!length(p[["needed"]])) {
      NULL
    } else {
      panel_sample_miss(p[["present"]], block)
    }
  })

  per_clock <- stats::setNames(
    vector("list", length(clock_sequence)),
    clock_sequence
  )
  score_miss <- per_clock
  norm_miss <- per_clock

  # raw per-panel miss for every non-alias clock
  for (i in seqi[!is_alias]) {
    id <- clock_sequence[[i]]
    score_miss[[id]] <- score_part_miss[[pidx[["score"]][["idx"]][[i]]]]
    norm_miss[[id]] <- norm_part_miss[[pidx[["norm"]][["idx"]][[i]]]]
  }

  pheno <- block[["pheno"]]
  female <- if (is.null(pheno)) NULL else as.numeric(pheno[["Female"]])

  # aliases: stitch score-panel miss from member counts
  for (i in seqi[is_alias]) {
    id <- clock_sequence[[i]]
    score_miss[[id]] <- stitch_routed_sample_miss(
      id,
      score_miss,
      female,
      sample_id
    )
  }

  # blank member rows this sex did not score
  rows <- sex_rows(female, length(sample_id))
  for (id in intersect(clock_sequence, names(routed[["sex"]]))) {
    sm <- score_miss[[id]]
    sm[!rows[[routed[["sex"]][[id]]]]] <- NA_integer_
    score_miss[[id]] <- sm
  }

  # records last (aliases keep NULL -- panels differ by member)
  for (i in seqi[!is_alias]) {
    id <- clock_sequence[[i]]
    per_clock[[id]] <- coverage_record(
      cpg_list[["per_clock"]][[id]],
      score_miss[[id]],
      norm_miss[[id]]
    )
  }

  list(
    per_clock = per_clock,
    sample_miss = list(score = score_miss, norm = norm_miss)
  )
}
