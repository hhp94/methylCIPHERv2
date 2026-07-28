# engine reduction (sum vs mean) via calc_clocks() output.

test_that("linear_mean clocks reduce by mean, not sum (EpiTOC regression)", {
  co <- clock_coefs("EpiTOC")
  ic <- clock_intercept("EpiTOC")
  DNAm <- random_betas(names(co), n = 6L)

  got <- calc_clocks(DNAm, "EpiTOC")$scores[, "EpiTOC"]
  mean_form <- ic + as.numeric(DNAm[, names(co)] %*% co) / length(co)
  expect_equal(unname(got), unname(mean_form), tolerance = 1e-10)

  sum_form <- ic + as.numeric(DNAm[, names(co)] %*% co)
  expect_false(isTRUE(all.equal(unname(got), unname(sum_form))))
})

test_that("plain linear clocks still reduce by sum (unchanged)", {
  co <- clock_coefs("Hannum")
  ic <- clock_intercept("Hannum")
  DNAm <- random_betas(names(co), n = 6L) # all present -> no absent, policy irrelevant

  got <- calc_clocks(DNAm, "Hannum")$scores[, "Hannum"]
  sum_form <- ic + as.numeric(DNAm[, names(co)] %*% co)
  expect_equal(unname(got), unname(sum_form), tolerance = 1e-10)
})

test_that("every catalog clock maps to a known score_type tag", {
  known <- c(
    "linear",
    "GrimAge",
    "DNAmFitAge",
    "PhysAge",
    "pack_linear",
    "pack_systemsage",
    "Dunedin",
    "normalized",
    "EpiTOC2",
    "MiAge",
    "Zhang2019",
    "sex_routed"
  )
  # score_type() routes every catalog clock (or stops)
  tags <- vapply(mc_index$clock_id, score_type, character(1))
  expect_true(all(tags %in% known))
})

test_that("a clock no branch claims is a hard stop, not a silent tag", {
  local_mocked_bindings(clock_type = function(id) "some_future_computation")
  expect_error(score_type("Hannum"))
})
