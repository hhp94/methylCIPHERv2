# dunedin: degraded-coverage and PACE QN paths not covered by parity

# vendor fill for fully-absent model CpGs. parity owns the value. this owns count and NA absence.
test_that("DunedinPoAm38 vendor-fills fully-absent CpGs (score_imputed_full)", {
  skip_on_cran()
  cpgs <- clock_scoring_cpgs("DunedinPoAm38")
  keep <- cpgs[seq_len(length(cpgs) - 2L)] # drop 2 of 46 -> still >= 80% covered
  DNAm <- random_betas(keep, n = 5L)
  res <- calc_clocks(DNAm, "DunedinPoAm38")

  cov <- res$coverage$per_clock[[1]]$DunedinPoAm38
  expect_equal(cov$score_imputed_full, 2L)
  expect_false(anyNA(res$scores[, "DunedinPoAm38"]))
})

# always-on proof that normalization is applied. reference golden is in test-fixtures-parity.R.
test_that("DunedinPACE quantile-normalizes the gold panel before the linear score", {
  skip_on_cran()
  skip_if_not_installed("betanorm")
  gold <- clock_norm_target("DunedinPACE")
  panel <- names(gold)
  DNAm <- random_betas(panel, n = 5L) # full gold-panel coverage, no gates/fill
  res <- calc_clocks(DNAm, "DunedinPACE")

  norm <- betanorm::quantile_norm(
    DNAm[, panel, drop = FALSE],
    target = as.numeric(gold[panel])
  )
  colnames(norm) <- panel
  coef <- clock_coefs("DunedinPACE")
  golden <- as.numeric(
    clock_intercept("DunedinPACE") + norm[, names(coef)] %*% coef
  )
  expect_equal(as.numeric(res$scores[, "DunedinPACE"]), golden)

  linear <- as.numeric(
    clock_intercept("DunedinPACE") + DNAm[, names(coef)] %*% coef
  )
  expect_false(isTRUE(all.equal(golden, linear)))
})

# normalizing clock keeps score- and norm-panel partial fills apart
test_that("DunedinPACE reports score and norm panel miss separately", {
  skip_on_cran()
  skip_if_not_installed("betanorm")
  norm_panel <- names(clock_norm_target("DunedinPACE"))
  score_panel <- clock_scoring_cpgs("DunedinPACE")
  norm_only <- setdiff(norm_panel, score_panel)

  DNAm <- random_betas(norm_panel, n = 4L)
  DNAm[1, norm_only[1]] <- NA_real_ # norm panel only
  DNAm[2, score_panel[1]] <- NA_real_ # in both panels (score subset of norm)

  res <- calc_clocks(DNAm, "DunedinPACE", min_samples_coverage = 0)
  sm <- res$coverage$sample_miss
  expect_equal(unname(sm$norm[, "DunedinPACE"]), c(1L, 1L, 0L, 0L))
  expect_equal(unname(sm$score[, "DunedinPACE"]), c(0L, 1L, 0L, 0L))

  cov <- res$coverage$per_clock[[1]]$DunedinPACE
  expect_true(cov$normalizes)
  # score panel holds one of the two holed cpgs, norm panel both.
  expect_equal(cov$score_imputed_partial, 1L)
  expect_equal(cov$norm_imputed_partial, 2L)
})

# qn needs the whole background panel, so an absent one is filled, not dropped
test_that("DunedinPACE counts fully-absent norm CpGs as norm_imputed_full", {
  skip_on_cran()
  skip_if_not_installed("betanorm")
  norm_panel <- names(clock_norm_target("DunedinPACE"))
  score_panel <- clock_scoring_cpgs("DunedinPACE")
  norm_only <- setdiff(norm_panel, score_panel)

  DNAm <- random_betas(setdiff(norm_panel, norm_only[1:3]), n = 4L)
  res <- calc_clocks(DNAm, "DunedinPACE")

  cov <- res$coverage$per_clock[[1]]$DunedinPACE
  expect_equal(cov$norm_imputed_full, 3L)
  expect_equal(cov$norm_dropped, 0L)
  # score panel untouched, and no absent norm CpG leaks into its counts
  expect_equal(cov$score_present, cov$score_needed)
  expect_false(anyNA(res$scores[, "DunedinPACE"]))
})

# coverage floors live in test-coverage-gate.R for all clocks
