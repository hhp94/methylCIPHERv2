# coverage/QC, once upstream of scoring, keyed by clock id

# per-sample partial-fill count over one panel's present CpGs
panel_sample_miss <- function(present, block) {
  count_sample_miss(
    block[["DNAm"]],
    cached_cols(present, block[["partial_cache"]])
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
