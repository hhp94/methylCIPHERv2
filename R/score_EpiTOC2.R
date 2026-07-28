# epiTOC2 (tnsc): mean of 2*(beta-beta0)/(delta*(1-beta0)) over present CpGs
score_EpiTOC2 <- function(id, cpgs, block, results) {
  params <- epitoc2_params(id)
  present <- cpgs[["score_present"]]
  delta <- params[["delta"]][present]
  beta0 <- params[["beta0"]][present]
  coef <- 1 / (delta * (1 - beta0))

  lp <- linear_predictor(
    coef = coef,
    intercept = 0,
    cov_coefs = numeric(0),
    score_present = present,
    block = block,
    id = id
  )

  ground <- sum(coef * beta0)
  score_matrix(
    2 * (as.numeric(lp[["cpg_contrib"]]) - ground) / length(present),
    block[["sample_id"]],
    id
  )
}

# epiTOC2 per-CpG ground state: named delta (de-novo rate) and beta0 vectors
epitoc2_params <- function(id) {
  tab <- require_fields(
    component_tensor(id, "cpg"),
    c("cpg", "delta", "beta0"),
    "params tensor",
    id
  )
  list(
    delta = stats::setNames(as.numeric(tab[["delta"]]), tab[["cpg"]]),
    beta0 = stats::setNames(as.numeric(tab[["beta0"]]), tab[["cpg"]])
  )
}
