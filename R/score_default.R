# shared linear scorers

# Horvath age back-transform, at the published adult age of 20
ADULT_AGE <- 20

anti_trafo <- function(x) {
  ifelse(x < 0, (1 + ADULT_AGE) * exp(x) - 1, (1 + ADULT_AGE) * x + ADULT_AGE)
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
  block,
  id,
  observed = NULL
) {
  # `observed` lets a pre-transform branch supply already-normalized betas
  obs <- if (is.null(observed)) {
    observed_panel(score_present, block)
  } else {
    observed
  }
  cpg_contrib <- obs[["values"]] %*% coef[obs[["cols"]]]

  cov_contrib <- 0
  if (length(cov_coefs)) {
    need <- names(cov_coefs)
    pheno <- block[["pheno"]]
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
linear_score <- function(cpgs, block, observed = NULL) {
  id <- cpgs[["clock_id"]]
  reduction <- clock_reduction(id)
  coef <- clock_coefs(id)
  fill <- absent_fill(id, coef, cpgs[["score_absent"]])

  lp <- linear_predictor(
    coef = coef,
    intercept = clock_intercept(id),
    cov_coefs = clock_covariates_coefs(id),
    score_present = cpgs[["score_present"]],
    block = block,
    id = id,
    observed = observed
  )

  if (identical(reduction, "mean")) {
    n_terms <- length(cpgs[["score_present"]]) + length(fill[["filled"]])
    cpg_num <- lp[["cpg_contrib"]] + fill[["offset"]]
    linpred <- cpg_num / n_terms + lp[["cov_contrib"]] + clock_intercept(id)
  } else {
    linpred <- lp[["linpred"]] + fill[["offset"]]
  }

  transform <- resolve_output_transform(clock_output_transform(id))
  score_matrix(transform(linpred), block[["sample_id"]], id)
}
