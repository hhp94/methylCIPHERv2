# zhang2019: full-matrix sample z-score, then linear sum over the arm's present
# CpGs. routes both arms -- EN bundled, BLUP coefs from its pack
score_Zhang2019 <- function(id, cpgs, block, results) {
  coef <- clock_coefs(id, block[["packs"]])
  present <- cpgs[["score_present"]]

  # sample_scale is over the full input matrix, not the panel union -- banked by
  # mc_cohort(), the one place that still sees every column
  mom <- block[["sample_moments"]]
  m <- mom[["mean"]]
  s <- mom[["sd"]]

  lp <- linear_predictor(
    coef = coef,
    intercept = 0,
    cov_coefs = numeric(0),
    score_present = present,
    block = block
  )
  csum <- sum(coef[present])
  z_sum <- (as.numeric(lp[["cpg_contrib"]]) - m * csum) / s
  score_matrix(clock_intercept(id) + z_sum, block[["sample_id"]], id)
}
