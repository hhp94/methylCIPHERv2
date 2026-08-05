# grimAge: missing Age/Female yields NA scores and an up-front warning.

grim_pheno <- function(ids, na_row = integer(0), na_cols = character(0)) {
  ph <- mc_pheno(
    ids,
    Age = mc_ages(length(ids)),
    Female = rep(c(1L, 0L), length.out = length(ids))
  )
  if ("Age" %in% na_cols) {
    ph$Age[na_row] <- NA_real_
  }
  if ("Female" %in% na_cols) {
    ph$Female[na_row] <- NA_integer_
  }
  ph
}

grim_betas <- function(id, n = 4L) random_betas(clock_cpgs(id), n = n)

# family is derived. one na_cols set only. Age-only arm is the next test.
for (id in Filter(
  function(i) score_type(i) == "GrimAge",
  resolve_clocks("all")
)) {
  test_that(paste(id, "NA-s out a sample with missing Age and Female"), {
    skip_on_cran()
    DNAm <- grim_betas(id)
    ids <- rownames(DNAm)
    bad <- 3L

    expect_warning(
      res <- calc_clocks(
        DNAm,
        id,
        pheno = grim_pheno(ids, bad, c("Age", "Female"))
      )
    )
    got <- res$scores[, id]

    expect_equal(rownames(res$scores), ids)
    expect_true(is.na(got[bad]))
    expect_true(all(is.finite(got[-bad])))
  })
}

# surrogate path resolves cpgs from the declared panel.
test_that("a component coefficient outside the declared panel is a hard stop", {
  skip_on_cran()
  cpgs <- list(
    clock_id = "fake",
    score_needed = c("cg1", "cg2", "cg3"),
    score_present = c("cg1", "cg3")
  )
  coef <- c(cg1 = 0.5, cg3 = -0.5)

  expect_equal(unname(component_present(coef, cpgs, "fake")), c("cg1", "cg3"))
  expect_error(component_present(c(coef, cg9 = 1), cpgs, "fake"))
})

test_that("GrimAge surrogates go NA only where they consume the missing covariate", {
  skip_on_cran()
  id <- "GrimAgeV1"
  deps <- clock_depends_on(id)
  DNAm <- grim_betas(id)
  bad <- 3L

  expect_warning(
    res <- calc_clocks(DNAm, id, pheno = grim_pheno(rownames(DNAm), bad, "Age"))
  )
  sc <- res$scores

  uses_age <- vapply(
    deps,
    function(d) "Age" %in% names(clock_covariates_coefs(d)),
    logical(1L)
  )
  expect_true(any(uses_age) && !all(uses_age)) # both arms are exercised

  expect_true(all(is.na(sc[bad, deps[uses_age]])))
  expect_true(all(is.finite(sc[bad, deps[!uses_age]])))
})

test_that("the missing-covariate warning only fires for scored rows", {
  skip_on_cran()
  id <- "GrimAgeV1"
  DNAm <- grim_betas(id)
  ids <- rownames(DNAm)

  expect_no_warning(calc_clocks(DNAm, id, pheno = grim_pheno(ids)))

  # an NA on a pheno row no sample matches is not this cohort's problem
  extra <- grim_pheno(c(ids, "unscored"))
  extra$Age[nrow(extra)] <- NA_real_
  expect_no_warning(calc_clocks(DNAm, id, pheno = extra))
})
