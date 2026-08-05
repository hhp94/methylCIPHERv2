# epiTOC2 branch. parity owns the clean-panel golden. this owns degraded-panel and coverage.

test_that("absent EpiTOC2 CpGs drop out rather than fill", {
  skip_on_cran()
  panel <- clock_scoring_cpgs("EpiTOC2")
  DNAm <- random_betas(clock_cpgs(c("EpiTOC2", "HypoClock")), n = 6L)
  drop <- panel[1:10]
  DNAm2 <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]

  res <- calc_clocks(DNAm2, "EpiTOC2")
  expect_true(all(is.finite(res$scores[, "EpiTOC2"])))

  cov <- res$coverage$per_clock[[1]][["EpiTOC2"]]
  expect_equal(cov$score_dropped, 10L)
  expect_equal(cov$score_imputed_full, 0L)
})
