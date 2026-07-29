# dunedin: degraded-coverage and PACE QN paths not covered by parity

# vendor fill for fully-absent model CpGs
test_that("DunedinPoAm38 vendor-fills fully-absent CpGs (score_imputed_full)", {
  cpgs <- clock_scoring_cpgs("DunedinPoAm38")
  keep <- cpgs[seq_len(length(cpgs) - 2L)] # drop 2 of 46 -> still >= 80% covered
  DNAm <- random_betas(keep, n = 5L)
  res <- calc_clocks(DNAm, "DunedinPoAm38")

  cov <- res$coverage$per_clock$DunedinPoAm38
  expect_equal(cov$score_imputed_full, 2L)
  expect_false(anyNA(res$scores[, "DunedinPoAm38"]))

  coef <- clock_coefs("DunedinPoAm38")
  ref <- clock_impute_ref("DunedinPoAm38")
  absent <- setdiff(names(coef), keep)
  golden <- as.numeric(
    clock_intercept("DunedinPoAm38") +
      DNAm[, keep] %*% coef[keep] +
      sum(coef[absent] * ref[absent])
  )
  expect_equal(as.numeric(res$scores[, "DunedinPoAm38"]), golden)
})

# PACE QN golden (always-on proof that normalization runs)
test_that("DunedinPACE quantile-normalizes the gold panel before the linear score", {
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
  skip_if_not_installed("betanorm")
  norm_panel <- names(clock_norm_target("DunedinPACE"))
  score_panel <- clock_scoring_cpgs("DunedinPACE")
  norm_only <- setdiff(norm_panel, score_panel)

  DNAm <- random_betas(norm_panel, n = 4L)
  DNAm[1, norm_only[1]] <- NA_real_ # norm panel only
  DNAm[2, score_panel[1]] <- NA_real_ # in both panels (score subset of norm)

  res <- calc_clocks(DNAm, "DunedinPACE", min_samples_coverage = 0)
  sm <- res$coverage$sample_miss
  expect_equal(colnames(sm$score), "DunedinPACE")
  expect_equal(colnames(sm$norm), "DunedinPACE")

  expect_equal(unname(sm$norm[, "DunedinPACE"]), c(1L, 1L, 0L, 0L))
  expect_equal(unname(sm$score[, "DunedinPACE"]), c(0L, 1L, 0L, 0L))

  cov <- res$coverage$per_clock$DunedinPACE
  expect_true(cov$normalizes)
  # score panel holds one of the two holed CpGs, norm panel both (score is a
  # subset of norm)
  expect_equal(cov$score_imputed_partial, 1L)
  expect_equal(cov$norm_imputed_partial, 2L)
})

# coverage floors live in test-coverage-gate.R for all clocks
