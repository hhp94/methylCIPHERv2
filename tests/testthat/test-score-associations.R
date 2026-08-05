# score_associations(): advisory age-correlation frame. delete with R/score_associations.R.

assoc_fixture <- function(clocks = "Hannum", n = 20L) {
  DNAm <- random_betas(clock_cpgs(clocks), n = n)
  list(
    result = calc_clocks(DNAm, clocks),
    age = stats::runif(n, 20, 80)
  )
}

# a one-row stand-in reference with the expectation we want to test against
assoc_ref <- function(id, age_r, lo, hi) {
  data.frame(
    clock = id,
    age_r = age_r,
    age_r_lo = lo,
    age_r_hi = hi,
    stringsAsFactors = FALSE
  )
}

test_that("one row per clock the reference covers", {
  fx <- assoc_fixture()
  out <- score_associations(fx$result, age = fx$age)

  expect_equal(
    names(out),
    c(
      "clock_id",
      "n",
      "obs_age_r",
      "exp_age_r",
      "exp_lo",
      "exp_hi",
      "outside",
      "wrong_sign"
    )
  )
  expect_equal(out$clock_id, "Hannum")
  expect_equal(out$n, 20L)
  expect_true(is.finite(out$obs_age_r))

  # age falls back to the Age column the record carries
  n <- 20L
  DNAm <- random_betas(clock_cpgs("DNAmADM"), n = n)
  pheno <- mc_pheno(rownames(DNAm), Age = stats::runif(n, 20, 80))
  carried <- score_associations(calc_clocks(DNAm, "DNAmADM", pheno = pheno))
  expect_equal(carried$clock_id, "DNAmADM")
})

test_that("degenerate cohorts drop out rather than correlating", {
  fx <- assoc_fixture(n = 20L)
  # fewer than MIN_ASSOC_N usable pairs
  age <- rep(NA_real_, 20L)
  age[1:3] <- c(30, 50, 70)
  expect_equal(nrow(score_associations(fx$result, age = age)), 0L)

  # constant age has no correlation to take
  expect_equal(nrow(score_associations(fx$result, age = rep(50, 20L))), 0L)
})

test_that("bad input is refused", {
  fx <- assoc_fixture()

  expect_error(score_associations(fx$result$scores, age = fx$age))
  expect_error(score_associations(fx$result))
  expect_error(score_associations(fx$result, age = fx$age[-1]))
})

test_that("the flags read off the supplied reference", {
  skip_on_cran()
  fx <- assoc_fixture()
  obs <- score_associations(fx$result, age = fx$age)$obs_age_r

  # outside is the interval test
  covered <- assoc_ref("Hannum", obs, obs - 0.2, obs + 0.2)
  missed <- assoc_ref("Hannum", obs, obs + 0.1, obs + 0.3)
  expect_false(assoc_report(fx$result, fx$age, covered)$outside)
  expect_true(assoc_report(fx$result, fx$age, missed)$outside)

  # wrong_sign needs a reference strong enough to have a sign
  strong <- assoc_ref("Hannum", -sign(obs) * 0.9, NA_real_, NA_real_)
  weak <- assoc_ref("Hannum", -sign(obs) * 0.1, NA_real_, NA_real_)
  expect_true(assoc_report(fx$result, fx$age, strong)$wrong_sign)
  expect_false(assoc_report(fx$result, fx$age, weak)$wrong_sign)
})
