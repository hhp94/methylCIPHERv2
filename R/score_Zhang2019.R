# zhang2019: full-matrix sample z-score, then linear sum over present EN CpGs
score_Zhang2019 <- function(id, cpgs, block, results) {
  coef <- clock_coefs(id)
  present <- cpgs[["score_present"]]

  # sample_scale is over the full input matrix, not the panel union
  full <- block[["DNAm_full"]]
  m <- matrixStats::rowMeans2(full, na.rm = TRUE)
  s <- matrixStats::rowSds(full, na.rm = TRUE)

  lp <- linear_predictor(
    coef = coef,
    intercept = 0,
    cov_coefs = numeric(0),
    score_present = present,
    block = block,
    id = id
  )
  csum <- sum(coef[present])
  z_sum <- (as.numeric(lp[["cpg_contrib"]]) - m * csum) / s
  score_matrix(clock_intercept(id) + z_sum, block[["sample_id"]], id)
}
