# epiTOC2 branch. group request also scores HypoClock -- cover both panels.

group_panel <- function() {
  unique(unlist(lapply(c("EpiTOC2", "HypoClock"), clock_scoring_cpgs)))
}

# 2 * mean((beta - beta0) / (delta * (1 - beta0))) over `cpgs`.
epitoc2_form <- function(DNAm, cpgs) {
  p <- epitoc2_params("EpiTOC2")
  delta <- p$delta[cpgs]
  beta0 <- p$beta0[cpgs]
  coef <- 1 / (delta * (1 - beta0))
  2 * (as.numeric(DNAm[, cpgs] %*% coef) - sum(coef * beta0)) / length(cpgs)
}

test_that("EpiTOC2 scores the ground-state-corrected mean", {
  panel <- clock_scoring_cpgs("EpiTOC2")
  DNAm <- random_betas(group_panel(), n = 6L)

  got <- calc_clocks(DNAm, "EpiTOC2")$scores[, "EpiTOC2"]
  expect_equal(unname(got), epitoc2_form(DNAm, panel), tolerance = 1e-10)

  p <- epitoc2_params("EpiTOC2")
  coef <- 1 / p$delta[panel]
  no_ground <- 2 * as.numeric(DNAm[, panel] %*% coef) / length(panel)
  expect_false(isTRUE(all.equal(unname(got), no_ground)))
})

test_that("absent EpiTOC2 CpGs drop out of the mean", {
  panel <- clock_scoring_cpgs("EpiTOC2")
  DNAm <- random_betas(group_panel(), n = 6L)
  drop <- panel[1:10]
  present <- setdiff(panel, drop)
  DNAm2 <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]

  res <- calc_clocks(DNAm2, "EpiTOC2")
  expect_equal(
    unname(res$scores[, "EpiTOC2"]),
    epitoc2_form(DNAm2, present),
    tolerance = 1e-10
  )

  cov <- res$coverage$per_clock[["EpiTOC2"]]
  expect_equal(cov$score_dropped, 10L)
  expect_equal(cov$score_imputed_full, 0L)
})
