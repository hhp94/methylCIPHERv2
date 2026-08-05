# check_DNAm(): orientation, EPICv2/MSA suffixes, and shape refusals.

test_that("a malformed DNAm is refused", {
  x <- random_betas(mc_fake_cpgs(40), n = 6)
  expect_error(check_DNAm(as.data.frame(x)))
  expect_error(check_DNAm(seq_len(10)))

  # sample ids are mandatory
  no_ids <- x
  rownames(no_ids) <- NULL
  expect_error(check_DNAm(no_ids))
})

test_that("orientation is reported, and a tall cohort is not mistaken for it", {
  expect_silent(check_DNAm(random_betas(mc_fake_cpgs(40), n = 6)))

  # regex kept: transposed and unrecognizable-names are two distinct warnings.
  expect_warning(
    check_DNAm(t(random_betas(mc_fake_cpgs(40), n = 6))),
    "transposed"
  )

  # one clock over a large cohort is legitimately taller than it is wide
  expect_silent(check_DNAm(random_betas(mc_fake_cpgs(3), n = 20)))
  unnamed <- random_betas(mc_fake_cpgs(3), n = 20)
  colnames(unnamed) <- paste0("v", seq_len(ncol(unnamed)))
  expect_warning(check_DNAm(unnamed))
})

test_that("replicate suffixes are reported, sparse ones included", {
  skip_on_cran()
  # the scan is a bounded sample, so a sparse case is the one worth pinning
  ids <- mc_fake_cpgs(400)
  ids[c(5, 200, 399)] <- paste0(ids[c(5, 200, 399)], "_TC11")
  expect_warning(check_DNAm(random_betas(ids, n = 3)))

  # retroelement panels carry ch... ids with underscores that are not replicate addresses.
  quiet <- c(
    mc_fake_cpgs(20),
    "ch.13.39564907R_II_R_O_37491",
    "ch.2.30415474F_II_F_O_37488"
  )
  expect_silent(check_DNAm(random_betas(quiet, n = 6)))
})
