# normalize=: per-clock normalization decision, resolved before any DNAm read

# BMIQ needs a genuinely multi-modal panel: on U(0,1) noise the three fitted
# components overlap and the boundary between them has no density crossing, so
# samples fail. Jitter the gold standard instead, so the input is
# methylation-shaped.
methylation_betas <- function(gold, n = 4L) {
  panel <- names(gold)
  m <- matrix(
    rep(as.numeric(gold), each = n),
    nrow = n,
    dimnames = list(paste0("sample", seq_len(n)), panel)
  )
  pmin(pmax(m + stats::rnorm(length(m), sd = 0.05), 0.001), 0.999)
}

# Horvath1 and Knight declare BMIQ but default to declining it, so the default
# score is the plain linear one. Parity skips Horvath1 (the horvath block), so
# this re-derivation is its only numeric gate.
test_that("Horvath1 defaults to no normalization and scores raw", {
  panel <- clock_scoring_cpgs("Horvath1")
  DNAm <- random_betas(panel, n = 5L)
  res <- calc_clocks(DNAm, "Horvath1")

  coef <- clock_coefs("Horvath1")
  tf <- resolve_output_transform(clock_output_transform("Horvath1"))
  golden <- as.numeric(tf(
    clock_intercept("Horvath1") + DNAm[, names(coef)] %*% coef
  ))
  expect_equal(as.numeric(res$scores[, "Horvath1"]), golden)

  cov <- res$coverage$per_clock$Horvath1
  expect_false(cov$normalizes)
  expect_equal(cov$norm_imputed_partial, 0L)
})

# declining is the absence of a panel, not a special case downstream
test_that("a declined clock asks for no normalization panel", {
  expect_equal(length(clock_norm_cpgs("Horvath1", FALSE)), 0L)
  expect_equal(length(clock_norm_cpgs("Horvath1", TRUE)), 21368L)

  # the 21k gold panel never reaches the required CpG set
  DNAm <- random_betas(clock_scoring_cpgs("Knight"), n = 4L)
  res <- calc_clocks(DNAm, "Knight")
  expect_false(res$coverage$per_clock$Knight$normalizes)
  expect_false("Knight" %in% colnames(res$coverage$sample_miss$norm))
})

# sim_DNAm builds over the same resolved decision
test_that("sim_DNAm omits the gold panel when normalization is declined", {
  sim <- sim_DNAm("Horvath1", n = 3L)
  expect_equal(ncol(sim$DNAm), length(clock_scoring_cpgs("Horvath1")))
})

test_that("constitutive normalization cannot be declined", {
  expect_error(resolve_normalize(c(DunedinPACE = FALSE), "DunedinPACE"))
  expect_true(resolve_normalize(NULL, "DunedinPACE")[["DunedinPACE"]])
})

test_that("normalize= refuses requests the catalog cannot express", {
  # no scheme declared at all
  expect_error(resolve_normalize(c(Hannum = TRUE), "Hannum"))
  # noob is an IDAT-level correction, unreachable from a beta matrix
  expect_error(resolve_normalize(c(Horvath2 = TRUE), "Horvath2"))
  # declining a scheme the clock never declared is merely redundant
  expect_false(resolve_normalize(c(Hannum = FALSE), "Hannum")[["Hannum"]])
})

test_that("normalize= refuses clocks outside the run and bad shapes", {
  expect_error(resolve_normalize(c(Horvath1 = TRUE), "Hannum"))
  expect_error(resolve_normalize(c(TRUE, FALSE), c("Hannum", "Horvath1")))
  expect_error(resolve_normalize(c(Hannum = NA), "Hannum"))
})

# a bare policy is a wish, not a claim: it passes over clocks that cannot honor
test_that("an unnamed policy reaches only the clocks that can honor it", {
  got <- resolve_normalize(TRUE, c("Horvath1", "DunedinPACE", "Hannum"))
  expect_equal(unname(got), c(TRUE, TRUE, FALSE))
})

# BMIQ golden: proves calc_clocks applies the calibration over the declared
# panel at the fixed settings, then scores the calibrated betas
test_that("Horvath1 BMIQ-calibrates the gold panel before the linear score", {
  skip_if_not_installed("betanorm")
  gold <- clock_norm_target("Horvath1")
  panel <- names(gold)
  DNAm <- methylation_betas(gold, n = 4L)
  res <- calc_clocks(DNAm, "Horvath1", normalize = c(Horvath1 = TRUE))

  fit <- betanorm::bmiq_calibration(
    DNAm[, panel, drop = FALSE],
    goldstandard.beta = as.numeric(gold[panel]),
    nfit = length(panel),
    verbose = FALSE,
    on.sample.error = "continue",
    failed.sample = "NA"
  )
  coef <- clock_coefs("Horvath1")
  tf <- resolve_output_transform(clock_output_transform("Horvath1"))
  golden <- as.numeric(tf(
    clock_intercept("Horvath1") + fit$calibrated[, names(coef)] %*% coef
  ))
  expect_equal(as.numeric(res$scores[, "Horvath1"]), golden)

  # calibration is not a no-op
  raw <- as.numeric(tf(
    clock_intercept("Horvath1") + DNAm[, names(coef)] %*% coef
  ))
  expect_false(isTRUE(all.equal(golden, raw)))
})

# the record must be able to say which Horvath1 it holds
test_that("provenance records which clocks were normalized", {
  DNAm <- random_betas(clock_scoring_cpgs("Horvath1"), n = 3L)
  res <- calc_clocks(DNAm, "Horvath1")
  expect_equal(res$provenance$normalized, character(0))

  skip_if_not_installed("betanorm")
  gold <- clock_norm_target("Horvath1")
  norm <- calc_clocks(
    methylation_betas(gold, n = 3L),
    "Horvath1",
    normalize = c(Horvath1 = TRUE)
  )
  expect_equal(norm$provenance$normalized, "Horvath1")

  # constitutive normalization is recorded even though nobody asked for it
  pace <- calc_clocks(
    random_betas(names(clock_norm_target("DunedinPACE")), n = 3L),
    "DunedinPACE",
    min_samples_coverage = 0
  )
  expect_equal(pace$provenance$normalized, "DunedinPACE")
})

test_that("a normalizing clock gains a norm panel and its coverage column", {
  skip_if_not_installed("betanorm")
  gold <- clock_norm_target("Horvath1")
  DNAm <- methylation_betas(gold, n = 4L)
  res <- calc_clocks(DNAm, "Horvath1", normalize = c(Horvath1 = TRUE))

  cov <- res$coverage$per_clock$Horvath1
  expect_true(cov$normalizes)
  expect_true("Horvath1" %in% colnames(res$coverage$sample_miss$norm))
})

# Coverage is computed before any score, so a sample whose mixture will not fit
# gets an NA score with its panel still reported whole -- indistinguishable on a
# saved record from any other NA. A constant sample has no three-component
# mixture to fit, so the failure is provoked rather than mocked. The warning
# comes off the same branch as the record, so asserting the record covers both.
test_that("a sample BMIQ cannot fit is on the record, not a bare NA", {
  skip_if_not_installed("betanorm")
  gold <- clock_norm_target("Horvath1")
  DNAm <- methylation_betas(gold, n = 4L)
  DNAm[2, ] <- 0.5

  res <- suppressWarnings(
    calc_clocks(DNAm, "Horvath1", normalize = c(Horvath1 = TRUE))
  )
  got <- res$scores[, "Horvath1"]

  expect_equal(res$provenance$scoring_failures$Horvath1, rownames(DNAm)[2])
  expect_true(is.na(got[[2]]))
  expect_false(anyNA(got[-2]))

  # coverage says the panel was whole, which is exactly why provenance has to
  # carry this: nothing in the counts distinguishes that NA from absent input
  cov <- res$coverage$per_clock$Horvath1
  expect_equal(cov$score_used, cov$score_needed)
})

# the empty case is the normal one, and it has to be an empty list rather than
# an absent field -- a reader must not have to test for NULL
test_that("a clean run records no scoring failures", {
  DNAm <- random_betas(clock_scoring_cpgs("Hannum"), n = 3L)
  res <- calc_clocks(DNAm, "Hannum")
  expect_equal(res$provenance$scoring_failures, list())
})

# absent background CpGs are dropped from the fit, never filled from the target
test_that("BMIQ drops absent background CpGs rather than filling them", {
  skip_if_not_installed("betanorm")
  gold <- clock_norm_target("Horvath1")
  full <- methylation_betas(gold, n = 4L)
  # drop background-only probes: the scoring panel stays whole, so no gate fires
  norm_only <- setdiff(names(gold), clock_scoring_cpgs("Horvath1"))
  dropped <- norm_only[seq_len(2000L)]
  thin <- full[, setdiff(names(gold), dropped), drop = FALSE]

  res <- calc_clocks(thin, "Horvath1", normalize = c(Horvath1 = TRUE))
  cov <- res$coverage$per_clock$Horvath1
  expect_equal(cov$norm_needed, length(gold))
  expect_equal(cov$norm_present, length(gold) - 2000L)
  expect_false(anyNA(res$scores[, "Horvath1"]))

  coef <- clock_coefs("Horvath1")
  tf <- resolve_output_transform(clock_output_transform("Horvath1"))
  bmiq_score <- function(m) {
    fit <- betanorm::bmiq_calibration(
      m,
      goldstandard.beta = as.numeric(gold[colnames(m)]),
      nfit = ncol(m),
      verbose = FALSE,
      on.sample.error = "continue",
      failed.sample = "NA"
    )
    as.numeric(tf(
      clock_intercept("Horvath1") + fit$calibrated[, names(coef)] %*% coef
    ))
  }
  got <- as.numeric(res$scores[, "Horvath1"])
  expect_equal(got, bmiq_score(thin))

  # the rejected alternative: filling those 2000 from the target would bias the
  # sample's own mixture fit toward the gold standard
  filled <- full
  filled[, dropped] <- rep(as.numeric(gold[dropped]), each = nrow(full))
  expect_false(isTRUE(all.equal(got, bmiq_score(filled[, names(gold)]))))
})
