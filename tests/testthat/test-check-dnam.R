# check_DNAm(): orientation, EPICv2/MSA replicate suffixes, and shape refusals.

cpg_ids <- function(n) sprintf("cg%08d", seq_len(n))

test_that("a well-formed matrix passes quietly", {
  expect_silent(check_DNAm(random_betas(cpg_ids(400), n = 20)))
})

test_that("probe ids in the rows are reported as transposed", {
  x <- t(random_betas(cpg_ids(400), n = 20))
  expect_warning(check_DNAm(x), "transposed")
})

test_that("more rows than columns is not on its own transposed", {
  # one clock over a large cohort is legitimately taller than it is wide
  expect_silent(check_DNAm(random_betas(cpg_ids(353), n = 1000)))
})

test_that("unrecognizable column names on a tall matrix warn", {
  x <- random_betas(cpg_ids(20), n = 400)
  colnames(x) <- paste0("v", seq_len(ncol(x)))
  expect_warning(check_DNAm(x))
})

test_that("EPICv2/MSA replicate suffixes are reported", {
  ids <- paste0(cpg_ids(400), rep(c("_TC11", "_BC21"), length.out = 400))
  expect_warning(check_DNAm(random_betas(ids, n = 20)), "suffix")
})

test_that("a suffix is found even when few columns carry one", {
  # the scan is a bounded sample, so a sparse case is the one worth pinning
  ids <- cpg_ids(400)
  ids[c(5, 200, 399)] <- paste0(ids[c(5, 200, 399)], "_TC11")
  expect_warning(check_DNAm(random_betas(ids, n = 20)), "suffix")
})

test_that("unsuffixed arrays and underscored non-CpG probes stay quiet", {
  # the Retroelement panels carry ch... ids with underscores that are not
  # replicate addresses
  ids <- c(
    cpg_ids(20),
    "ch.13.39564907R_II_R_O_37491",
    "ch.2.30415474F_II_F_O_37488"
  )
  expect_silent(check_DNAm(random_betas(ids, n = 20)))
})

test_that("a data.frame is refused", {
  expect_error(check_DNAm(as.data.frame(random_betas(cpg_ids(400), n = 20))))
})

test_that("a transposed data.frame reports orientation before refusing", {
  x <- as.data.frame(t(random_betas(cpg_ids(400), n = 20)))
  expect_warning(expect_error(check_DNAm(x)), "transposed")
})

test_that("objects without two dimensions are refused", {
  expect_error(check_DNAm(seq_len(10)))
})

test_that("sample ids remain mandatory", {
  x <- random_betas(cpg_ids(400), n = 20)
  rownames(x) <- NULL
  expect_error(check_DNAm(x))
})
