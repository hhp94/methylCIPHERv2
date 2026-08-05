# the family's two arms, in the order the karyotype spec declares them
WANG <- unname(karyotype_inputs(karyotype_spec()))

# panels plus a slice of the shared z-score reference.
wang_DNAm <- function(n = 6L, n_ref = 300L) {
  ref <- clock_sample_scale_ref(WANG[[1]])
  cpgs <- union(
    unlist(lapply(WANG, clock_scoring_cpgs)),
    sample(ref, n_ref)
  )
  random_betas(cpgs, n = n)
}

test_that("both Wang members score off one shared moment domain", {
  skip_on_cran()
  DNAm <- wang_DNAm()
  res <- calc_clocks(DNAm, WANG, pheno = mc_pheno(rownames(DNAm)))

  expect_equal(colnames(res$scores), WANG)
  expect_true(all(is.finite(res$scores)))
  # the ref is not a panel, so it never widens coverage
  cov <- res$coverage$per_clock[[1]]
  expect_equal(cov[[WANG[[1]]]]$score_needed, cov[[WANG[[1]]]]$score_present)
  expect_equal(cov[[WANG[[2]]]]$score_dropped, 0)
})

test_that("clock_cpgs reports the declared moment ref, so a sim scores", {
  skip_on_cran()
  # the ref is input, not a panel. matrix must carry it or the score is NA.
  ref <- clock_sample_scale_ref(WANG[[1]])
  expect_true(all(ref %in% clock_cpgs(WANG)))

  sim <- sim_DNAm(WANG, n = 4L)
  res <- calc_clocks(sim$DNAm, WANG, pheno = sim$pheno)
  expect_true(all(is.finite(res$scores)))
})

test_that("a mixed request keeps each clock on its own moment domain", {
  skip_on_cran()
  # full-matrix domain and a declared ref in one sweep. cross-contamination is only visible here.
  ids <- c("Zhang2019EN", WANG[[1]])
  sim <- sim_DNAm(ids, n = 4L)
  # the full-matrix note rides along on every call. Not what is under test.
  both <- suppressMessages(calc_clocks(sim$DNAm, ids, pheno = sim$pheno))

  for (id in ids) {
    alone <- suppressMessages(calc_clocks(sim$DNAm, id, pheno = sim$pheno))
    expect_equal(both$scores[, id], alone$scores[, id])
  }
})

test_that("a sample with no z-score reference is scored NA, not silently", {
  skip_on_cran()
  # panels only: the declared ref meets nothing that was measured
  cpgs <- unlist(lapply(WANG, clock_scoring_cpgs))
  DNAm <- random_betas(cpgs, n = 4L)
  pheno <- mc_pheno(rownames(DNAm))

  expect_warning(
    res <- calc_clocks(DNAm, WANG[[1]], pheno = pheno)
  )
  expect_true(all(is.na(res$scores)))
  # coverage cannot see this -- the scoring panel itself is complete
  expect_equal(
    res$provenance$scoring_failures[[WANG[[1]]]],
    rownames(DNAm)
  )
})
