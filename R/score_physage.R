# DNAmPhysAge: surrogate means, reverse-code, cohort z-score, row sum.
score_physage <- function(id, cpgs, DNAm, partial_cache = NULL) {
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)
  if (n < 2L) {
    stop(
      "score_physage(): '",
      id,
      "' is batch-dependent (cohort z-score) and needs >= 2 samples; got ",
      n,
      ".",
      call. = FALSE
    )
  }

  surrogates <- physage_surrogates(id)
  ref <- clock_impute_ref(id)

  # Each surrogate: mean of coef*beta over present + vendor-filled absent.
  raws <- vapply(
    surrogates,
    function(s) {
      coef <- s$coef
      present <- intersect(names(coef), cpgs$score_present)
      absent <- setdiff(names(coef), present)
      miss_ref <- setdiff(absent, names(ref))
      if (length(miss_ref)) {
        stop(
          "score_physage(): '",
          id,
          "' surrogate '",
          s$name,
          "' absent CpG(s) lack a vendor mean: ",
          paste(utils::head(miss_ref, 5L), collapse = ", "),
          call. = FALSE
        )
      }
      lp <- linear_predictor(
        coef = coef,
        intercept = 0,
        cov_coefs = numeric(0),
        score_present = present,
        DNAm = DNAm,
        partial_cache = partial_cache,
        id = s$name
      )
      absent_offset <- if (length(absent)) {
        sum(coef[absent] * ref[absent])
      } else {
        0
      }
      raw <- (as.numeric(lp$cpg_contrib) + absent_offset) / length(coef)
      if (s$negate) -raw else raw
    },
    numeric(n)
  )

  z <- scale(raws)
  phys <- rowSums(z)

  poly <- physage_poly_coef(id)
  score_vec <- if (is.null(poly)) {
    phys
  } else {
    poly_eval(as.numeric(scale(phys)), poly)
  }

  score <- matrix(
    as.numeric(score_vec),
    nrow = n,
    ncol = 1L,
    dimnames = list(sample_id, id)
  )

  cached <- if (is.null(partial_cache)) {
    character(0)
  } else {
    intersect(cpgs$score_present, colnames(partial_cache))
  }
  sample_miss <- if (length(cached)) {
    slideimp::mat_miss(DNAm[, cached, drop = FALSE], col = FALSE)
  } else {
    integer(n)
  }
  names(sample_miss) <- sample_id

  coverage <- list(
    clock_id = id,
    policy = clock_impute(id)[["policy"]],
    score_needed = length(cpgs$score_needed),
    score_present = length(cpgs$score_present),
    score_used = length(cpgs$score_needed),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = length(cpgs$score_absent),
    score_dropped = 0L,
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = cpgs$score_absent
  )

  list(score = score, coverage = coverage, sample_miss = sample_miss)
}

# y = sum_k coef[k+1] * x^k (lowest degree first).
poly_eval <- function(x, coef) {
  powers <- vapply(seq_along(coef) - 1L, function(k) x^k, numeric(length(x)))
  as.numeric(powers %*% coef)
}
