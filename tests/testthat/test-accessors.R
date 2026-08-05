# declaration resolution via declared pointers. internal accessors. off CRAN.

CATALOG_IDS <- names(mc_catalog)

SCALE_IDS <- CATALOG_IDS[vapply(
  CATALOG_IDS,
  function(i) length(sample_scale_steps(i)) > 0L,
  logical(1)
)]

# the sample_scale clocks, split on whether they declare a moment ref
WITH_REF <- SCALE_IDS[vapply(
  SCALE_IDS,
  function(i) !is.null(clock_sample_scale_ref(i)),
  logical(1)
)]
WITHOUT_REF <- setdiff(SCALE_IDS, WITH_REF)

test_that("full-panel is exactly a sample_scale step with no declared ref", {
  skip_on_cran()
  expect_gt(length(SCALE_IDS), 0L) # the loop is vacuous otherwise
  for (id in SCALE_IDS) {
    expect_equal(
      clock_needs_full_panel(id),
      is.null(clock_sample_scale_ref(id))
    )
  }

  # a clock with no sample_scale step has none of the three facts
  plain <- setdiff(CATALOG_IDS, SCALE_IDS)[[1]]
  expect_false(clock_needs_full_panel(plain))
  expect_null(clock_sample_scale_ref(plain))
  expect_null(clock_moment_domain(plain))
})

test_that("a moment domain keys on the declaration, and a shared ref collapses", {
  skip_on_cran()
  expect_equal(resolve_moment_domains(character(0)), list())

  # a ref-less clock spells its domain "full", with no CpGs to declare
  for (id in WITHOUT_REF) {
    d <- clock_moment_domain(id)
    expect_equal(d[["key"]], "full")
    expect_null(d[["cpgs"]])
  }
  # a declared ref is a closed set, so it gets its own key
  for (id in WITH_REF) {
    d <- clock_moment_domain(id)
    expect_match(d[["key"]], ":", fixed = TRUE)
    expect_equal(d[["cpgs"]], clock_sample_scale_ref(id))
  }

  skip_if(!length(WITHOUT_REF) || !length(WITH_REF))
  keys <- vapply(WITH_REF, function(i) clock_moment_domain(i)[["key"]], "")
  got <- resolve_moment_domains(c(WITHOUT_REF[[1]], WITH_REF))
  expect_equal(names(got), c("full", unique(unname(keys))))
})

test_that("an unresolvable shared name errors rather than searching", {
  skip_on_cran()
  entry <- list(shared = list(a = list(name = "a", file = "weights/a.csv.gz")))
  expect_equal(shared_named(entry, "a", "X")[["file"]], "weights/a.csv.gz")
  expect_error(shared_named(entry, "missing", "X"))
  expect_error(shared_named(list(), "a", "X"))
})

test_that("a full-panel clock announces that it reads every column", {
  skip_on_cran()
  expect_message(say_full_panel_clocks("Zhang2019EN"))
  expect_equal(
    suppressMessages(say_full_panel_clocks("Zhang2019EN")),
    "Zhang2019EN"
  )
})
