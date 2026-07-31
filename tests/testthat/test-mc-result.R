# mc_result S3 verb behavior

test_that("calc_clocks returns an mc_result, and as.matrix gives naked scores", {
  DNAm <- random_betas(clock_scoring_cpgs("Hannum"), n = 3L)
  res <- calc_clocks(DNAm, "Hannum")

  expect_s3_class(res, "mc_result")
  expect_equal(as.matrix(res), res$scores)
})
