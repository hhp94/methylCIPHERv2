# grimAge: missing Age/Female yields NA scores and an up-front warning.

grim_pheno <- function(ids, na_row = integer(0), na_cols = character(0)) {
  ph <- data.frame(
    ID = ids,
    Age = seq(40, 70, length.out = length(ids)),
    Female = rep(c(1L, 0L), length.out = length(ids)),
    stringsAsFactors = FALSE
  )
  if ("Age" %in% na_cols) {
    ph$Age[na_row] <- NA_real_
  }
  if ("Female" %in% na_cols) {
    ph$Female[na_row] <- NA_integer_
  }
  ph
}

grim_betas <- function(id, n = 4L) {
  random_betas(
    panels_union(clock_panels(resolve_clocks_sequence(resolve_clocks(id)))),
    n = n
  )
}

for (id in c("GrimAgeV1", "GrimAgeV2")) {
  test_that(paste(id, "NA-s out a sample with missing Age and/or Female"), {
    DNAm <- grim_betas(id)
    ids <- rownames(DNAm)
    bad <- 3L

    for (na_cols in list("Age", "Female", c("Age", "Female"))) {
      res <- NULL # expect_warning() returns the condition, so capture by assignment
      expect_warning(
        res <- calc_clocks(DNAm, id, pheno = grim_pheno(ids, bad, na_cols))
      )
      got <- res$scores[, id]

      expect_equal(rownames(res$scores), ids)
      expect_true(is.na(got[bad]))
      expect_true(all(is.finite(got[-bad])))
    }
  })
}

# surrogate path resolves CpGs from the declared panel (as PhysAge does), so a
# component can never score a CpG coverage did not count
test_that("a component coefficient outside the declared panel is a hard stop", {
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
  id <- "GrimAgeV1"
  deps <- clock_depends_on(id)
  DNAm <- grim_betas(id)
  bad <- 3L

  res <- NULL
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

test_that("the missing-covariate warning names each affected column, and only fires for scored rows", {
  id <- "GrimAgeV1"
  DNAm <- grim_betas(id)
  ids <- rownames(DNAm)

  expect_warning(
    calc_clocks(
      DNAm,
      id,
      pheno = grim_pheno(ids, c(2L, 3L), c("Age", "Female"))
    ),
    "Age.*Female"
  )
  expect_no_warning(calc_clocks(DNAm, id, pheno = grim_pheno(ids)))

  extra <- grim_pheno(c(ids, "unscored"))
  extra$Age[nrow(extra)] <- NA_real_
  expect_no_warning(calc_clocks(DNAm, id, pheno = extra))
})
