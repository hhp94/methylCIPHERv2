# Dunedin pace-of-aging family (PoAm linear; PACE QN then linear).
score_dunedin <- function(
  id,
  cpgs,
  DNAm,
  partial_cache = NULL,
  min_coverage = 0.75
) {
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)

  coef <- clock_coefs(id)
  intercept <- clock_intercept(id)
  model_needed <- cpgs$score_needed
  model_present <- cpgs$score_present
  model_absent <- cpgs$score_absent

  # Panel for impute/(PACE) QN: gold panel for PACE, model CpGs for PoAm.
  qn <- identical(clock_norm_scheme(id), "quantile")
  if (qn) {
    fill_ref <- dunedin_gold_means(id)
    panel_needed <- cpgs$norm_needed
    panel_present <- cpgs$norm_present
    panel_absent <- cpgs$norm_absent
  } else {
    fill_ref <- clock_impute_ref(id)
    panel_needed <- model_needed
    panel_present <- model_present
    panel_absent <- model_absent
  }

  score <- matrix(NA_real_, nrow = n, ncol = 1L, dimnames = list(sample_id, id))

  cached <- if (is.null(partial_cache)) {
    character(0)
  } else {
    intersect(panel_present, colnames(partial_cache))
  }
  raw_cols <- setdiff(panel_present, cached)
  sample_miss <- if (length(cached)) {
    slideimp::mat_miss(DNAm[, cached, drop = FALSE], col = FALSE)
  } else {
    integer(n)
  }
  names(sample_miss) <- sample_id

  # Per-sample NA-gate. Whole-clock coverage already stopped in check_coverage().
  raw_miss <- if (length(panel_present)) {
    slideimp::mat_miss(DNAm[, panel_present, drop = FALSE], col = FALSE)
  } else {
    integer(n)
  }
  not_obs <- raw_miss + if (qn) 0L else length(panel_absent)
  sample_cov <- 1 - not_obs / length(panel_needed)
  low_sample <- if (min_coverage > 0) sample_cov < min_coverage else logical(n)

  # Assemble panel: present from cache/raw, absent from fill ref.
  panel <- matrix(
    0,
    nrow = n,
    ncol = length(panel_needed),
    dimnames = list(sample_id, panel_needed)
  )
  if (length(cached)) {
    panel[, cached] <- partial_cache[, cached, drop = FALSE]
  }
  if (length(raw_cols)) {
    panel[, raw_cols] <- DNAm[, raw_cols, drop = FALSE]
  }
  if (length(panel_absent)) {
    panel[, panel_absent] <- rep(fill_ref[panel_absent], each = n)
  }

  scored <- if (qn) {
    if (!requireNamespace("betanorm", quietly = TRUE)) {
      stop(
        "score_dunedin(): '",
        id,
        "' quantile-normalizes and needs the 'betanorm' package. Install it ",
        "(Remotes: hhp94/betanorm) to score this clock.",
        call. = FALSE
      )
    }
    norm <- betanorm::quantile_norm(
      panel,
      target = as.numeric(fill_ref[panel_needed])
    )
    dimnames(norm) <- dimnames(panel)
    norm
  } else {
    panel
  }
  score[, 1] <- as.numeric(
    intercept + scored[, model_needed, drop = FALSE] %*% coef[model_needed]
  )

  # Mask samples that observe too little of the panel.
  if (any(low_sample)) {
    score[low_sample, 1] <- NA_real_
    warning(
      id,
      ": ",
      sum(low_sample),
      " of ",
      n,
      " sample(s) returned NA -- each observes under ",
      round(100 * min_coverage),
      "% of the scoring panel.",
      call. = FALSE
    )
  }

  coverage <- list(
    clock_id = id,
    policy = clock_impute(id)[["policy"]],
    score_needed = length(model_needed),
    score_present = length(model_present),
    score_used = length(model_needed),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = length(model_absent),
    score_dropped = 0L,
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = model_absent
  )

  list(score = score, coverage = coverage, sample_miss = sample_miss)
}
