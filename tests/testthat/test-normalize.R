# normalize=: per-clock normalization decision, resolved before any DNAm read

GOLD <- clock_norm_target("Horvath1")

# bmiq needs multi-modal input -- jitter the gold standard, not U(0,1).
# `background` thins the gold panel to scoring panel plus that many background probes.
methylation_betas <- function(gold = GOLD, n = 4L, background = NULL) {
  panel <- names(gold)
  if (!is.null(background)) {
    score <- clock_scoring_cpgs("Horvath1")
    panel <- c(score, setdiff(panel, score)[seq_len(background)])
  }
  m <- matrix(
    rep(as.numeric(gold[panel]), each = n),
    nrow = n,
    dimnames = list(paste0("sample", seq_len(n)), panel)
  )
  pmin(pmax(m + stats::rnorm(length(m), sd = 0.05), 0.001), 0.999)
}

# the linear half of Horvath1, over whatever matrix it is handed
horvath1_score <- function(m) {
  coef <- clock_coefs("Horvath1")
  tf <- resolve_output_transform(clock_output_transform("Horvath1"))
  as.numeric(tf(clock_intercept("Horvath1") + m[, names(coef)] %*% coef))
}

bmiq_calibrated <- function(m) {
  betanorm::bmiq_calibration(
    m,
    goldstandard.beta = as.numeric(GOLD[colnames(m)]),
    nfit = ncol(m),
    verbose = FALSE,
    on.sample.error = "continue",
    failed.sample = "NA"
  )$calibrated
}

# front-door refusals: the user's own normalize= argument
test_that("normalize= refuses what the catalog cannot express", {
  # constitutive normalization cannot be declined
  expect_error(resolve_normalize(c(DunedinPACE = FALSE), "DunedinPACE"))
  # no scheme declared at all
  expect_error(resolve_normalize(c(Hannum = TRUE), "Hannum"))
  # noob is an IDAT-level correction, unreachable from a beta matrix
  expect_error(resolve_normalize(c(Horvath2 = TRUE), "Horvath2"))
  # a clock outside the run, an unnamed vector, a non-logical
  expect_error(resolve_normalize(c(Horvath1 = TRUE), "Hannum"))
  expect_error(resolve_normalize(c(TRUE, FALSE), c("Hannum", "Horvath1")))
  expect_error(resolve_normalize(c(Hannum = NA), "Hannum"))
})

test_that("normalize= resolves per clock, and a bare policy is only a wish", {
  expect_true(resolve_normalize(NULL, "DunedinPACE")[["DunedinPACE"]])
  # declining a scheme the clock never declared is merely redundant
  expect_false(resolve_normalize(c(Hannum = FALSE), "Hannum")[["Hannum"]])
  # an unnamed policy passes over the clocks that cannot honor it
  got <- resolve_normalize(TRUE, c("Horvath1", "DunedinPACE", "Hannum"))
  expect_equal(unname(got), c(TRUE, TRUE, FALSE))
})

test_that("Horvath1 defaults to declining normalization", {
  DNAm <- random_betas(clock_scoring_cpgs("Horvath1"), n = 5L)
  res <- calc_clocks(DNAm, "Horvath1")

  cov <- res$coverage$per_clock[[1]]$Horvath1
  expect_false(cov$normalizes)
  expect_equal(cov$norm_imputed_partial, 0L)
  expect_equal(res$provenance$normalized, character(0))
  # a clean run records no scoring failures: list(), not an absent field
  expect_equal(res$provenance$scoring_failures, list())
})

# declining is the absence of a panel, not a special case downstream
test_that("a declined clock asks for no normalization panel", {
  skip_on_cran()
  expect_equal(length(clock_norm_cpgs("Horvath1", FALSE)), 0L)
  # accepting asks for exactly the declared gold panel, whatever its size
  expect_equal(length(clock_norm_cpgs("Horvath1", TRUE)), length(GOLD))

  # the 21k gold panel never reaches the required CpG set
  DNAm <- random_betas(clock_scoring_cpgs("Knight"), n = 4L)
  res <- calc_clocks(DNAm, "Knight")
  expect_false(res$coverage$per_clock[[1]]$Knight$normalizes)
  expect_false("Knight" %in% colnames(res$coverage$sample_miss$norm))

  # and sim_DNAm builds over the same resolved decision
  sim <- sim_DNAm("Horvath1", n = 3L)
  expect_equal(ncol(sim$DNAm), length(clock_scoring_cpgs("Horvath1")))
})

# the norm half of sample_miss is keyed by clock, like the score half
test_that("a clock with no norm panel keeps its entry rather than losing it", {
  skip_on_cran()
  spec <- mc_spec(c("Hannum", "DunedinPACE"))
  DNAm <- random_betas(panels_union(spec$panels), n = 4L)
  miss <- score_cohort(DNAm, spec, mc_cohort(DNAm, spec))$coverage$sample_miss

  expect_equal(names(miss$norm), spec$sequence)
  expect_equal(names(miss$score), spec$sequence)
  # present but empty for the clock that does not normalize
  expect_null(miss$norm[["Hannum"]])
  expect_equal(length(miss$norm[["DunedinPACE"]]), nrow(DNAm))
})

# normalized arithmetic is in test-fixtures-parity.R. this file covers the record half.
test_that("a normalized run says on the record that it normalized", {
  skip_on_cran()
  skip_if_not_installed("betanorm")
  # record-only: gates off because a thinned background is deliberately short
  res <- calc_clocks(
    methylation_betas(background = 1000L),
    "Horvath1",
    normalize = c(Horvath1 = TRUE),
    min_clocks_coverage = 0,
    min_samples_coverage = 0
  )

  # the record must be able to say which Horvath1 it holds
  expect_equal(res$provenance$normalized, "Horvath1")
  cov <- res$coverage$per_clock[[1]]$Horvath1
  expect_true(cov$normalizes)
  expect_true("Horvath1" %in% colnames(res$coverage$sample_miss$norm))

  # constitutive normalization is recorded even though nobody asked for it
  pace <- calc_clocks(
    random_betas(names(clock_norm_target("DunedinPACE")), n = 3L),
    "DunedinPACE",
    min_samples_coverage = 0
  )
  expect_equal(pace$provenance$normalized, "DunedinPACE")
})

# unfit BMIQ sample: NA score + notes entry (coverage still full)
test_that("a sample BMIQ cannot fit is on the record, not a bare NA", {
  skip_on_cran()
  skip_if_not_installed("betanorm")
  # the failure is a property of the unfittable sample, not of the width
  DNAm <- methylation_betas(background = 1000L)
  DNAm[2, ] <- 0.5

  res <- suppressWarnings(
    calc_clocks(
      DNAm,
      "Horvath1",
      normalize = c(Horvath1 = TRUE),
      min_clocks_coverage = 0,
      min_samples_coverage = 0
    )
  )
  got <- res$scores[, "Horvath1"]

  expect_equal(res$provenance$scoring_failures$Horvath1, rownames(DNAm)[2])
  expect_true(is.na(got[[2]]))
  expect_false(anyNA(got[-2]))

  # coverage stays full -- notes is what distinguishes the NA
  cov <- res$coverage$per_clock[[1]]$Horvath1
  expect_equal(cov$score_used, cov$score_needed)
})

# absent background CpGs are dropped from the fit, never filled from the target
test_that("BMIQ drops absent background CpGs rather than filling them", {
  skip_on_cran()
  skip_if_not_installed("betanorm")
  full <- methylation_betas()
  # drop background-only probes: the scoring panel stays whole, so no gate fires
  norm_only <- setdiff(names(GOLD), clock_scoring_cpgs("Horvath1"))
  dropped <- norm_only[seq_len(2000L)]
  thin <- full[, setdiff(names(GOLD), dropped), drop = FALSE]

  res <- calc_clocks(thin, "Horvath1", normalize = c(Horvath1 = TRUE))
  cov <- res$coverage$per_clock[[1]]$Horvath1
  expect_equal(cov$norm_present, length(GOLD) - 2000L)
  # the record says dropped, not filled
  expect_equal(cov$norm_dropped, 2000L)
  expect_equal(cov$norm_imputed_full, 0L)
  expect_false(anyNA(res$scores[, "Horvath1"]))

  # and the fit really did ignore them: filling from the target moves the score
  filled <- full
  filled[, dropped] <- rep(as.numeric(GOLD[dropped]), each = nrow(full))
  expect_false(isTRUE(all.equal(
    as.numeric(res$scores[, "Horvath1"]),
    horvath1_score(bmiq_calibrated(filled[, names(GOLD)]))
  )))
})
