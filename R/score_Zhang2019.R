# zhang2019: full-matrix sample z-score, then linear sum. both arms
score_Zhang2019 <- function(id, cpgs, block, results) {
  coef <- clock_coefs(id, block[["packs"]])
  present <- cpgs[["score_present"]]

  # sample_scale over the full input matrix, banked by mc_cohort()
  mom <- block_domain_moments(block, id)
  m <- mom[["mean"]]
  s <- mom[["sd"]]

  lp <- linear_predictor(
    coef = coef,
    intercept = 0,
    cov_coefs = numeric(0),
    score_present = present,
    score_idx = cpgs[["score_present_idx"]],
    block = block
  )
  csum <- sum(coef[present])
  z_sum <- (as.numeric(lp[["cpg_contrib"]]) - m * csum) / s
  score_matrix(clock_intercept(id) + z_sum, block[["sample_id"]], id)
}
