# coverage/QC, once upstream of scoring, keyed by clock id

# per-sample observed count over column positions (callers are past the value gate)
row_observed <- function(DNAm, idx) {
  obs <- col_stats(DNAm, idx)[["row_obs"]]
  # NULL means the kernel bailed on overflow, which the value gate in
  # scan_missing_cpgs() already ruled out for these columns
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

# per-sample partial-fill count over a panel's cached (partly filled) CpGs
panel_sample_miss <- function(cached, block) {
  # counted on the raw matrix -- the cache has the NAs already filled
  count_sample_miss(block[["DNAm"]], block_cols(cached, block))
}

# a routed member's count on a row its sex did not score is not its count
mask_routed_rows <- function(miss, sex_key, rows) {
  if (is.null(miss) || is.na(sex_key)) {
    return(miss)
  }
  miss[!rows[[sex_key]]] <- NA_integer_
  miss
}

# one clock's coverage record (only writer of record fields). every count is a
# CpG count, partials included. per-sample axis lives in sample_miss, not here
coverage_record <- function(cpgs, score_partial, norm_partial = 0L) {
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
    score_imputed_partial = score_partial,
    score_imputed_full = if (fill) n_absent else 0L,
    score_dropped = if (fill) 0L else n_absent,
    norm_needed = length(cpgs[["norm_needed"]]),
    norm_present = length(cpgs[["norm_present"]]),
    norm_imputed_partial = norm_partial,
    missing_cpgs = cpgs[["score_absent"]]
  )
}

# full coverage for the compute sequence (per_clock records + per-panel sample_miss)
compute_coverage <- function(clock_sequence, cpg_list, block) {
  sample_id <- block[["sample_id"]]
  routed <- sex_routed_members()
  # a clock that reads no betas gets neither a record nor a per-sample count
  reads_cpgs <- vapply(clock_sequence, clock_reads_cpgs, logical(1L))
  seqi <- seq_along(clock_sequence)
  pidx <- cpg_list[["panel_index"]]

  # one cached-column set per distinct panel. "partly filled" for both axes,
  # so both derive from it and cannot disagree
  score_cached <- lapply(
    pidx[["score"]][["parts"]],
    function(p) cached_cols(p[["present"]], block[["partial_cache"]])
  )
  norm_cached <- lapply(pidx[["norm"]][["parts"]], function(p) {
    if (!length(p[["needed"]])) {
      NULL
    } else {
      cached_cols(p[["present"]], block[["partial_cache"]])
    }
  })

  # count each distinct panel once, then fan out via the index
  score_part_miss <- lapply(score_cached, panel_sample_miss, block = block)
  norm_part_miss <- lapply(norm_cached, function(cc) {
    if (is.null(cc)) NULL else panel_sample_miss(cc, block)
  })
  # probe axis: a panel CpG that was filled for any sample (NULL panel -> 0)
  score_part_cpgs <- lengths(score_cached)
  norm_part_cpgs <- lengths(norm_cached)

  per_clock <- stats::setNames(
    vector("list", length(clock_sequence)),
    clock_sequence
  )
  score_miss <- per_clock
  norm_miss <- per_clock

  pheno <- block[["pheno"]]
  female <- if (is.null(pheno)) NULL else as.numeric(pheno[["Female"]])
  rows <- sex_rows(female, length(sample_id))

  # every beta-reading clock: per-sample miss (masked for routed members on
  # rows their sex did not score) plus the record on the probe axis
  for (i in seqi[reads_cpgs]) {
    id <- clock_sequence[[i]]
    # NA when the clock routes no sex -- masks nothing
    sex_key <- routed[["sex"]][id]
    score_miss[[id]] <- mask_routed_rows(
      score_part_miss[[pidx[["score"]][["idx"]][[i]]]],
      sex_key,
      rows
    )
    # single bracket + list(): no-norm clock keeps its NULL element.
    # `[[<-` would delete the name
    norm_miss[id] <- list(mask_routed_rows(
      norm_part_miss[[pidx[["norm"]][["idx"]][[i]]]],
      sex_key,
      rows
    ))
    per_clock[[id]] <- coverage_record(
      cpg_list[["per_clock"]][[id]],
      score_part_cpgs[[pidx[["score"]][["idx"]][[i]]]],
      norm_part_cpgs[[pidx[["norm"]][["idx"]][[i]]]]
    )
  }

  list(
    per_clock = per_clock,
    sample_miss = list(score = score_miss, norm = norm_miss)
  )
}
