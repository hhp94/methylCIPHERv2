# Shared linear scorers. linear and linear_transformed differ only by output_transform.

# Horvath age back-transform (adult.age = 20).
anti_trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}

resolve_output_transform <- function(name) {
  switch(
    name,
    identity = function(x) x,
    anti.trafo = anti_trafo,
    stop(
      "Unknown output_transform '",
      name,
      "' -- add it to the registry in score.R.",
      call. = FALSE
    )
  )
}

# linpred = intercept + sum(coef * beta) + covariates.
linear_predictor <- function(
  coef,
  intercept,
  cov_coefs,
  score_present,
  DNAm,
  partial_cache = NULL,
  pheno = NULL,
  id = "<component>"
) {
  n <- nrow(DNAm)

  cached <- if (is.null(partial_cache)) {
    character(0)
  } else {
    intersect(score_present, colnames(partial_cache))
  }
  raw <- setdiff(score_present, cached)
  used_cols <- c(cached, raw)

  if (length(used_cols)) {
    sub <- cbind(
      partial_cache[, cached, drop = FALSE],
      DNAm[, raw, drop = FALSE]
    )
    cpg_contrib <- sub %*% coef[used_cols]
  } else {
    cpg_contrib <- matrix(0, nrow = n, ncol = 1L)
  }

  cov_contrib <- 0
  if (length(cov_coefs)) {
    need <- names(cov_coefs)
    if (is.null(pheno) || !all(need %in% names(pheno))) {
      stop(
        "linear_predictor(): '",
        id,
        "' needs covariate(s) ",
        paste(need, collapse = ", "),
        " but they are absent from `pheno`.",
        call. = FALSE
      )
    }
    cov_mat <- as.matrix(pheno[, need, drop = FALSE])
    cov_contrib <- cov_mat %*% cov_coefs[need]
  }

  linpred <- cpg_contrib + cov_contrib + intercept
  list(
    linpred = linpred,
    cpg_contrib = cpg_contrib,
    cov_contrib = cov_contrib,
    used_cols = used_cols,
    cached = cached
  )
}

# Linear engine for one cpg_coefficient clock.
linear_score <- function(
  cpgs,
  DNAm,
  partial_cache = NULL,
  pheno = NULL,
  packs = NULL
) {
  id <- cpgs$clock_id
  policy <- clock_impute(id)[["policy"]]
  reduction <- clock_reduction(id)
  coef <- clock_coefs(id, packs)
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)

  absent <- cpgs$score_absent
  vendor_mean <- length(absent) && identical(policy, "vendor_mean")
  if (length(absent) && !policy %in% c("omit", "drop", "vendor_mean")) {
    stop(
      "linear_score(): clock '",
      id,
      "' has unsupported imputation policy '",
      policy,
      "' for ",
      length(absent),
      " absent CpG(s).",
      call. = FALSE
    )
  }

  lp <- linear_predictor(
    coef = coef,
    intercept = clock_intercept(id),
    cov_coefs = clock_covariate_coefs(id),
    score_present = cpgs$score_present,
    DNAm = DNAm,
    partial_cache = partial_cache,
    pheno = pheno,
    id = id
  )

  if (vendor_mean) {
    ref <- clock_impute_ref(id, packs)
    miss_ref <- setdiff(absent, names(ref))
    if (length(miss_ref)) {
      stop(
        "linear_score(): '",
        id,
        "' absent CpG(s) lack a vendor mean (cannot fill): ",
        paste(utils::head(miss_ref, 5L), collapse = ", "),
        call. = FALSE
      )
    }
    absent_offset <- sum(coef[absent] * ref[absent])
    vendor_filled <- absent
    dropped <- character(0)
  } else {
    absent_offset <- 0
    vendor_filled <- character(0)
    dropped <- absent
  }

  # Mean reduces only the CpG term; sum reuses linpred + vendor offset.
  if (identical(reduction, "mean")) {
    n_terms <- length(cpgs$score_present) + length(vendor_filled)
    cpg_num <- lp$cpg_contrib + absent_offset
    linpred <- cpg_num / n_terms + lp$cov_contrib + clock_intercept(id)
  } else {
    linpred <- lp$linpred + absent_offset
  }

  transform <- resolve_output_transform(clock_output_transform(id))
  score <- matrix(
    as.numeric(transform(linpred)),
    nrow = n,
    ncol = 1L,
    dimnames = list(sample_id, id)
  )

  sample_miss <- if (length(lp$cached)) {
    slideimp::mat_miss(DNAm[, lp$cached, drop = FALSE], col = FALSE)
  } else {
    integer(n)
  }
  names(sample_miss) <- sample_id

  coverage <- list(
    clock_id = id,
    policy = policy,
    score_needed = length(cpgs$score_needed),
    score_present = length(cpgs$score_present),
    score_used = length(lp$used_cols) + length(vendor_filled),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = length(vendor_filled),
    score_dropped = length(dropped),
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = absent
  )

  list(score = score, coverage = coverage, sample_miss = sample_miss)
}
