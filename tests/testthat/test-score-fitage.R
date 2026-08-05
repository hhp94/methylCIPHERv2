# dnamFitAge engine wiring (synthetic betas). alias routing is only gated here.

# hand-compute one member (intercept + betas %*% coef + Age).
member_expected <- function(id, DNAm, rows, age) {
  coef <- clock_coefs(id)
  cov <- clock_covariates_coefs(id)
  out <- clock_intercept(id) +
    as.numeric(DNAm[rows, names(coef), drop = FALSE] %*% coef)
  if (length(cov)) {
    out <- out + cov[["Age"]] * age[rows]
  }
  out
}

test_that("the alias routes each sample to its own sex's model", {
  skip_on_cran()
  fx <- grip_fixture(c(1, 1, 1, 0, 0, 0))
  DNAm <- fx$DNAm
  female <- fx$female
  age <- fx$age

  got <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = fx$pheno)$scores
  f <- which(female == 1)
  m <- which(female == 0)

  expected <- numeric(6)
  expected[f] <- member_expected("DNAmGrip_wAge_Female", DNAm, f, age)
  expected[m] <- member_expected("DNAmGrip_wAge_Male", DNAm, m, age)
  expect_equal(
    unname(got[, "DNAmGrip_wAge"]),
    unname(expected),
    tolerance = 1e-9
  )

  # one column for the whole cohort -- the members stay internal
  expect_true(all(is.finite(got[, "DNAmGrip_wAge"])))
  expect_false(any(grepl("_(Female|Male)$", colnames(got))))
})

test_that("panel coverage lands on the members, never on the alias", {
  skip_on_cran()
  fx <- grip_fixture()
  fem <- fx$fem
  mal <- fx$mal

  res <- calc_clocks(fx$DNAm, "DNAmGrip_wAge", pheno = fx$pheno)
  cov <- res$coverage$per_clock[[1]]

  # the two panels differ in size, so no one count is true of every sample.
  expect_null(cov[["DNAmGrip_wAge"]])
  expect_equal(cov[["DNAmGrip_wAge_Female"]]$score_needed, length(fem))
  expect_equal(cov[["DNAmGrip_wAge_Male"]]$score_needed, length(mal))

  # sample_miss spans clocks that read cpgs, not returned columns.
  expect_setequal(
    colnames(res$coverage$sample_miss$score),
    c("DNAmGrip_wAge_Female", "DNAmGrip_wAge_Male")
  )
  # nothing normalizes here, so the norm matrix has no columns
  expect_equal(ncol(res$coverage$sample_miss$norm), 0L)
})

test_that("per-sample QC routes with the score; panel counts do not", {
  skip_on_cran()
  fx <- grip_fixture(c(1, 1, 1, 0, 0, 0))
  DNAm <- fx$DNAm

  # a CpG only the female model uses, blanked for one female and one male.
  fem_only <- setdiff(names(fx$fem), names(fx$mal))[1]
  DNAm[c(1, 4), fem_only] <- NA_real_

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = fx$pheno)

  # the alias reads no CpGs, so it counts nothing and gets no column
  miss <- res$coverage$sample_miss$score
  expect_false("DNAmGrip_wAge" %in% colnames(miss))

  # cohort-mean count sits on the member that owns the panel.
  expect_equal(
    unname(miss[, "DNAmGrip_wAge_Female"]),
    c(1L, 0L, 0L, NA, NA, NA)
  )

  # the same fill, counted per panel, stays on that member too.
  cov <- res$coverage$per_clock[[1]]
  expect_equal(cov[["DNAmGrip_wAge_Female"]]$score_imputed_partial, 1)
  expect_equal(cov[["DNAmGrip_wAge_Male"]]$score_imputed_partial, 0)
})

# front-door refusal. regex tells a routed member apart from an unknown token.
test_that("a routed member is not callable and the refusal names its alias", {
  expect_error(sim_DNAm("DNAmGrip_wAge_Female", n = 2L), "DNAmGrip_wAge")
})

test_that("routed members stay out of the callable pool", {
  skip_on_cran()
  pool <- resolve_clocks("all")
  expect_true("DNAmGrip_wAge" %in% pool)
  expect_false(any(names(sex_routed_members()$alias) %in% pool))
})

test_that("absent member CpGs vendor-fill from that sex's medians", {
  skip_on_cran()
  fx <- gait_holed_fixture()
  coef <- fx$coef
  drop <- fx$drop
  res <- calc_clocks(fx$DNAm, "DNAmGait_noAge", pheno = fx$pheno)

  # golden: fill drew on this sex's medians, not the sibling model
  present <- setdiff(names(coef), drop)
  expected <- clock_intercept(fx$id) +
    as.numeric(fx$DNAm[, present, drop = FALSE] %*% coef[present]) +
    sum(coef[drop] * fx$medians[drop])
  expect_equal(
    unname(res$scores[, "DNAmGait_noAge"]),
    unname(expected),
    tolerance = 1e-9
  )

  cov <- res$coverage$per_clock[[1]][[fx$id]]
  expect_equal(cov$score_imputed_full, 5L)
  expect_equal(cov$score_dropped, 0L)
})

# KDM stack is parity-owned. this checks the alias assembles one column over its deps.
test_that("DNAmFitAge assembles from its dependencies on every sample", {
  skip_on_cran()
  DNAm <- random_betas(clock_cpgs("DNAmFitAge"), n = 6L)
  pheno <- mc_pheno(
    rownames(DNAm),
    Age = mc_ages(6),
    Female = c(1, 1, 1, 0, 0, 0)
  )

  sc <- calc_clocks(DNAm, "DNAmFitAge", pheno = pheno)$scores
  expect_true(all(is.finite(sc[, "DNAmFitAge"])))
  expect_true(all(is.finite(sc[, "GrimAgeV1"])))
  expect_false(any(names(sex_routed_members()$alias) %in% colnames(sc)))
})

# a clock assembled from other clocks' scores counts nothing of its own
test_that("composites report no coverage; the CpG readers under them do", {
  skip_on_cran()
  seq_ids <- resolve_clocks_sequence(resolve_clocks("DNAmFitAge"))
  DNAm <- random_betas(clock_cpgs("DNAmFitAge"), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- mc_pheno(rownames(DNAm), Age = mc_ages(6), Female = female)

  # drop 3 CpGs the female fitness components use
  drop <- intersect(
    clock_scoring_cpgs("DNAmGait_noAge_Female"),
    colnames(DNAm)
  )[
    1:3
  ]
  DNAm2 <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]

  res <- calc_clocks(DNAm2, "DNAmFitAge", pheno = pheno)
  cov <- res$coverage$per_clock[[1]]

  # composites read no betas. no record and no samples_coverage rows.
  composites <- Filter(Negate(clock_reads_cpgs), seq_ids)
  expect_gt(length(composites), 0L)
  for (id in composites) {
    expect_null(cov[[id]])
  }
  expect_false(any(composites %in% samples_coverage(res)$clock_id))

  # the fill lands on the component that declared the panel
  gait <- cov[["DNAmGait_noAge_Female"]]
  expect_equal(gait$score_imputed_full, 3L)
  expect_equal(gait$score_dropped, 0L)
  expect_equal(gait$score_used, gait$score_needed)

  # and the composite still scores, on every sample
  expect_true(all(is.finite(res$scores[, "DNAmFitAge"])))
})
