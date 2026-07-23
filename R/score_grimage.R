# GrimAgeV1/V2: Cox stack of surrogates + Age/Female, then rescale to years.
score_grimage <- function(
  id,
  cpgs,
  results,
  usable,
  DNAm,
  partial_cache = NULL,
  pheno = NULL
) {
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)

  group_id <- clock_group_bundle(id)[["group_id"]]
  cox <- grimage_cox_coef(id)
  comps <- clock_components(id)
  stack_names <- names(cox)

  X <- matrix(
    0,
    nrow = n,
    ncol = length(stack_names),
    dimnames = list(sample_id, stack_names)
  )
  for (nm in stack_names) {
    if (nm %in% c("Age", "Female")) {
      if (is.null(pheno) || !nm %in% names(pheno)) {
        stop(
          "score_grimage(): '",
          id,
          "' needs covariate '",
          nm,
          "' but it is absent from `pheno`.",
          call. = FALSE
        )
      }
      X[, nm] <- as.numeric(pheno[[nm]])
    } else if (startsWith(nm, "_internal_")) {
      # V2 surrogate scored inline; absent CpGs dropped.
      comp <- Filter(function(c) identical(c[["name"]], nm), comps)
      if (length(comp) != 1L) {
        stop(
          "score_grimage(): ",
          id,
          " is missing internal surrogate '",
          nm,
          "'.",
          call. = FALSE
        )
      }
      comp <- comp[[1]]
      coef <- bundle_tensor(group_id, comp[["file"]])
      intercept <- if (is.null(comp[["intercept"]])) 0 else comp[["intercept"]]
      lp <- linear_predictor(
        coef = coef,
        intercept = intercept,
        cov_coefs = covariate_coefs_from(comp[["covariates"]]),
        score_present = intersect(names(coef), usable),
        DNAm = DNAm,
        partial_cache = partial_cache,
        pheno = pheno,
        id = nm
      )
      X[, nm] <- lp$linpred
    } else {
      # Dependency surrogate already scored upstream.
      r <- results[[nm]]
      if (is.null(r)) {
        stop(
          "score_grimage(): ",
          id,
          " depends on surrogate '",
          nm,
          "' but it was not computed upstream.",
          call. = FALSE
        )
      }
      X[, nm] <- as.numeric(r$score)
    }
  }

  cox_score <- X[, stack_names, drop = FALSE] %*% cox
  score <- matrix(
    as.numeric(grimage_rescale(cox_score, grimage_rescale_params(id))),
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
    score_used = length(cpgs$score_present),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = 0L,
    score_dropped = length(cpgs$score_absent),
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = cpgs$score_absent
  )

  list(score = score, coverage = coverage, sample_miss = sample_miss)
}

# Cox scale -> years.
grimage_rescale <- function(cox_score, params) {
  (cox_score - params[["m_cox"]]) /
    params[["sd_cox"]] *
    params[["sd_age"]] +
    params[["m_age"]]
}
