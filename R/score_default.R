# shared linear scorers

# Horvath age back-transform (adult.age = 20)
anti_trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}

resolve_output_transform <- function(name) {
  switch(
    name,
    identity = function(x) x,
    anti.trafo = anti_trafo,
    cli::cli_abort(
      "Unknown output_transform {.val {name}}.",
      call = NULL
    )
  )
}

# linpred = intercept + sum(coef * beta) + covariates
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
  obs <- observed_panel(score_present, DNAm, partial_cache)
  cpg_contrib <- obs$values %*% coef[obs$cols]

  cov_contrib <- 0
  if (length(cov_coefs)) {
    need <- names(cov_coefs)
    if (is.null(pheno) || !all(need %in% names(pheno))) {
      cli::cli_abort(
        c(
          "{.val {id}} needs {cli::qty(need)} pheno column{?s}
           {.field {need}}.",
          "i" = "Add {cli::qty(need)}{?it/them} to {.arg pheno}."
        ),
        call = NULL
      )
    }
    cov_mat <- as.matrix(pheno[, need, drop = FALSE])
    cov_contrib <- cov_mat %*% cov_coefs[need]
  }

  linpred <- cpg_contrib + cov_contrib + intercept
  list(
    linpred = linpred,
    cpg_contrib = cpg_contrib,
    cov_contrib = cov_contrib
  )
}

# linear engine for one cpg_coefficient clock
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

  absent <- cpgs$score_absent
  vendor_mean <- length(absent) && identical(policy, "vendor_mean")
  if (length(absent) && !policy %in% c("omit", "drop", "vendor_mean")) {
    cli::cli_abort(
      "Clock {.val {id}} has unsupported imputation policy {.val {policy}}
       for {length(absent)} absent CpG{?s}.",
      call = NULL
    )
  }

  lp <- linear_predictor(
    coef = coef,
    intercept = clock_intercept(id),
    cov_coefs = clock_covariates_coefs(id),
    score_present = cpgs$score_present,
    DNAm = DNAm,
    partial_cache = partial_cache,
    pheno = pheno,
    id = id
  )

  if (vendor_mean) {
    absent_offset <- vendor_offset(
      coef,
      absent,
      clock_impute_ref(id, packs),
      id
    )
    vendor_filled <- absent
  } else {
    absent_offset <- 0
    vendor_filled <- character(0)
  }

  if (identical(reduction, "mean")) {
    n_terms <- length(cpgs$score_present) + length(vendor_filled)
    cpg_num <- lp$cpg_contrib + absent_offset
    linpred <- cpg_num / n_terms + lp$cov_contrib + clock_intercept(id)
  } else {
    linpred <- lp$linpred + absent_offset
  }

  transform <- resolve_output_transform(clock_output_transform(id))
  score_matrix(transform(linpred), sample_id, id)
}
