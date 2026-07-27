# DNAmFitAge engine wiring (synthetic betas).

fitage_pheno <- function(ids, female, age = NULL) {
  ph <- data.frame(
    ID = ids,
    Female = as.integer(female),
    stringsAsFactors = FALSE
  )
  if (!is.null(age)) {
    ph$Age <- age
  }
  ph
}

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
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  age <- seq(45, 70, length.out = 6)
  pheno <- fitage_pheno(rownames(DNAm), female, age)

  expect_equal(names(clock_covariates_coefs("DNAmGrip_wAge_Female")), "Age")

  got <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)$scores
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
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 4L)
  pheno <- fitage_pheno(rownames(DNAm), c(1, 1, 0, 0), rep(50, 4))

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)
  cov <- res$coverage$per_clock

  # the two panels differ in size, so no one count is true of every sample.
  expect_null(cov[["DNAmGrip_wAge"]])
  expect_equal(cov[["DNAmGrip_wAge_Female"]]$score_needed, length(fem))
  expect_equal(cov[["DNAmGrip_wAge_Male"]]$score_needed, length(mal))
  expect_false(length(fem) == length(mal))

  # members are scored but never returned as columns.
  expect_false(any(
    c("DNAmGrip_wAge_Female", "DNAmGrip_wAge_Male") %in%
      colnames(res$scores)
  ))
  # sample_miss is per panel -- the score matrix spans every returned column
  expect_equal(
    colnames(res$coverage$sample_miss$score),
    colnames(res$scores)
  )
  # no returned clock normalizes here, so the norm matrix has no columns
  expect_equal(ncol(res$coverage$sample_miss$norm), 0L)
})

test_that("per-sample QC routes with the score; panel counts do not", {
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- fitage_pheno(rownames(DNAm), female, rep(50, 6))

  # A CpG only the female model uses, blanked for one female and one male.
  fem_only <- setdiff(names(fem), names(mal))[1]
  DNAm[c(1, 4), fem_only] <- NA_real_

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)

  # only the female with a blanked value leans on a cohort mean
  expect_equal(
    unname(res$coverage$sample_miss$score[, "DNAmGrip_wAge"]),
    c(1L, 0L, 0L, 0L, 0L, 0L)
  )

  # the same fill, counted per panel, stays on the member that owns the panel.
  cov <- res$coverage$per_clock
  expect_equal(cov[["DNAmGrip_wAge_Female"]]$score_imputed_partial, 1)
  expect_equal(cov[["DNAmGrip_wAge_Male"]]$score_imputed_partial, 0)
})

test_that("routed members are not directly callable and name their alias", {
  expect_true("DNAmGrip_wAge" %in% resolve_clocks("all"))
  expect_false("DNAmGrip_wAge_Female" %in% resolve_clocks("all"))
  expect_false(any(grepl("_(Female|Male)$", resolve_clocks("DNAmFitAge"))))
  expect_error(resolve_clocks("DNAmGrip_wAge_Female"), "DNAmGrip_wAge")
})

test_that("absent member CpGs vendor-fill from that sex's medians", {
  id <- "DNAmGait_noAge_Female"
  coef <- clock_coefs(id)
  medians <- clock_impute_ref(id)
  full <- union(names(coef), names(clock_coefs("DNAmGait_noAge_Male")))
  DNAm <- random_betas(full, n = 4L)
  pheno <- fitage_pheno(rownames(DNAm), rep(1L, 4L))

  drop <- intersect(names(coef), names(medians))[1:5]
  DNAm2 <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]
  res <- calc_clocks(DNAm2, "DNAmGait_noAge", pheno = pheno)

  # golden: fill drew on this sex's medians, not the sibling model
  present <- setdiff(names(coef), drop)
  expected <- clock_intercept(id) +
    as.numeric(DNAm2[, present, drop = FALSE] %*% coef[present]) +
    sum(coef[drop] * medians[drop])
  expect_equal(
    unname(res$scores[, "DNAmGait_noAge"]),
    unname(expected),
    tolerance = 1e-9
  )

  cov <- res$coverage$per_clock[[id]]
  expect_equal(cov$score_imputed_full, 5L)
  expect_equal(cov$score_dropped, 0L)
})


test_that("DNAmFitAge mixes same-sex members by KDM and carries no batch stamp", {
  seq_ids <- resolve_clocks_sequence(resolve_clocks("DNAmFitAge"))
  DNAm <- random_betas(panels_union(clock_panels(seq_ids)), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- fitage_pheno(rownames(DNAm), female, seq(40, 65, length.out = 6))

  res <- calc_clocks(DNAm, "DNAmFitAge", pheno = pheno)
  sc <- res$scores
  expect_true(all(is.finite(sc[, "DNAmFitAge"])))

  expect_true(all(is.finite(sc[, "GrimAgeV1"])))

  # a sex-routed KDM component is read through its alias
  alias <- sex_routed_members()$alias
  column_for <- function(component) {
    if (component %in% names(alias)) alias[[component]] else component
  }

  expected <- rep(NA_real_, nrow(sc))
  for (sx in c("female", "male")) {
    member <- clock_routing("DNAmFitAge")[[sx]]
    rows <- if (sx == "female") which(female == 1) else which(female == 0)
    kdm <- fitage_kdm_params(member)
    acc <- numeric(length(rows))
    for (i in seq_len(nrow(kdm))) {
      acc <- acc +
        kdm$weight[i] *
          (sc[rows, column_for(kdm$component[i])] - kdm$center[i]) /
          kdm$scale[i]
    }
    expected[rows] <- acc
  }
  expect_equal(unname(sc[, "DNAmFitAge"]), unname(expected), tolerance = 1e-9)
})

test_that("the alias declares the routing covariate and its members do not", {
  expect_equal(clock_covariates_required("DNAmFitAge"), "Female")
  expect_equal(clock_covariates_required("DNAmFitAge_Female"), character(0))
})

# composite vendor-fills over its own panel, not the family prep panel
test_that("the composite vendor-fills over its own panel", {
  seq_ids <- resolve_clocks_sequence(resolve_clocks("DNAmFitAge"))
  full <- panels_union(clock_panels(seq_ids))
  DNAm <- random_betas(full, n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- fitage_pheno(rownames(DNAm), female, seq(40, 65, length.out = 6))

  # drop 3 CpGs from the female composite panel so imputed_full > 0
  drop <- intersect(clock_scoring_cpgs("DNAmFitAge_Female"), colnames(DNAm))[
    1:3
  ]
  DNAm2 <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]

  res <- calc_clocks(DNAm2, "DNAmFitAge", pheno = pheno)

  for (id in c("DNAmFitAge_Female", "DNAmFitAge_Male")) {
    cov <- res$coverage$per_clock[[id]]
    expect_gt(cov$score_imputed_full, 0L)
    expect_equal(cov$score_dropped, 0L)
    expect_equal(cov$score_used, cov$score_needed)
  }
})
