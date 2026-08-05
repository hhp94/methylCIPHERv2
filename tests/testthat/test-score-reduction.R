# engine routing. sum-vs-mean goldens live in parity.

test_that("every catalog clock maps to a known score_type tag", {
  skip_on_cran()
  # score_type() stops on a pair no branch claims. routing every clock without error is the invariant.
  expect_no_error(vapply(mc_index$clock_id, score_type, character(1)))
})

test_that("a clock no branch claims is a hard stop, not a silent tag", {
  skip_on_cran()
  local_mocked_bindings(clock_type = function(id) "some_future_computation")
  expect_error(score_type("Hannum"))
})
