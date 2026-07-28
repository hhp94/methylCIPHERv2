# PhysAge: vendor-mean fill + cohort_zscore composites (synthetic).

physage_union <- function() {
  members <- mc_groups[["PhysAge"]]$members
  unique(unlist(lapply(members, clock_scoring_cpgs)))
}

# fill offset lands inside the mean's divisor (plain mean golden is elsewhere)
test_that("absent surrogate CpGs vendor-fill by offset under mean reduction", {
  DNAm <- random_betas(physage_union(), n = 6L)
  co <- clock_coefs("DNAmCRP")
  ic <- clock_intercept("DNAmCRP")

  ref <- clock_impute_ref("DNAmCRP")
  drop <- intersect(names(co), names(ref))[1:5]
  DNAm2 <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]
  res2 <- calc_clocks(DNAm2, "DNAmCRP")
  present <- setdiff(names(co), drop)
  filled <- ic +
    (as.numeric(DNAm2[, present] %*% co[present]) + sum(co[drop] * ref[drop])) /
      length(co)
  expect_equal(
    unname(res2$scores[, "DNAmCRP"]),
    unname(filled),
    tolerance = 1e-10
  )

  cov <- res2$coverage$per_clock[["DNAmCRP"]]
  expect_equal(cov$score_imputed_full, 5L)
  expect_equal(cov$score_dropped, 0L)
})

test_that("PhysAge composites run and need >= 2 samples", {
  DNAm <- random_betas(physage_union(), n = 6L)
  res <- calc_clocks(DNAm, c("DNAmPhysAge", "DNAmPhysAge_years"))

  expect_true(all(
    c("DNAmPhysAge", "DNAmPhysAge_years") %in% colnames(res$scores)
  ))
  expect_true(all(is.finite(res$scores[, "DNAmPhysAge"])))
  expect_true(all(is.finite(res$scores[, "DNAmPhysAge_years"])))

  one <- random_betas(physage_union(), n = 1L)
  expect_error(calc_clocks(one, "DNAmPhysAge"))
})

test_that("a surrogate that goes constant stops instead of NaN-ing the cohort", {
  # DNAmPulsePr is 60 CpGs of a 1711-CpG panel, so losing all of it still
  # clears the default min_clocks_coverage. Every sample then vendor-fills to
  # the same value, scale() returns a NaN column, and rowSums() spreads it to
  # every sample of every other surrogate.
  surr <- physage_surrogates("DNAmPhysAge")
  flat <- names(surr[[which(vapply(
    surr,
    function(s) identical(s[["name"]], "raw_DNAmPulsePr"),
    logical(1)
  ))]][["coef"]])

  panel <- setdiff(physage_union(), flat)
  DNAm <- random_betas(panel, n = 6L)
  expect_error(suppressWarnings(calc_clocks(DNAm, "DNAmPhysAge")))
})
