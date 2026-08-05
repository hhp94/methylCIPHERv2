# mc_result S3 verb behavior

test_that("calc_clocks returns a well-formed mc_result", {
  DNAm <- random_betas(clock_scoring_cpgs("Hannum"), n = 3L)
  res <- calc_clocks(DNAm, "Hannum")

  # a record over list, never a matrix subclass
  expect_s3_class(res, "mc_result")
  expect_type(res, "list")
  expect_true(all(
    c("scores", "pheno", "coverage", "provenance") %in% names(res)
  ))
  expect_equal(dim(res$scores), c(3L, 1L))
  expect_equal(nrow(res$pheno), 3L)

  # as.matrix gives the naked scores
  expect_equal(as.matrix(res), res$scores)
})
