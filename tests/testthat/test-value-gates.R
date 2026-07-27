# The value gates over the one col_stats() sweep: +/-Inf stops, values outside
# [0, 1] warn. Asserted through calc_clocks() rather than the kernel wherever
# possible, because what matters is that a user handing over bad data hears
# about it.

gate_betas <- function(n = 8L) {
  spec <- mc_spec("Hannum")
  panel <- panels_union(spec$panels)
  list(DNAm = random_betas(panel, n = n), panel = panel)
}

# collect every warning instead of just the first -- the two range flags are
# independent and a matrix can trip both
warnings_of <- function(expr) {
  seen <- character(0)
  withCallingHandlers(
    expr,
    warning = function(cond) {
      seen <<- c(seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  seen
}

test_that("an infinite value stops scoring, whatever else is in the matrix", {
  b <- gate_betas()

  # no NA anywhere: `anyNA()` does not see an Inf, so a scan gated on it
  # would skip this entirely
  inf_only <- b$DNAm
  inf_only[1, b$panel[1]] <- Inf
  expect_error(calc_clocks(inf_only, "Hannum"))

  # ... and the same Inf must still stop when an unrelated NA is present.
  # These two disagreed once: the NA-free matrix scored -Inf while the other
  # silently filled the Inf as if it were missing.
  plus_na <- inf_only
  plus_na[2, b$panel[2]] <- NA
  expect_error(calc_clocks(plus_na, "Hannum"))

  neg <- b$DNAm
  neg[3, b$panel[5]] <- -Inf
  expect_error(calc_clocks(neg, "Hannum"))
})

test_that("the abort names where the infinite value is", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[3, b$panel[5]] <- Inf

  # wording is normally not asserted (CLAUDE.md, "Test altitude"), but
  # reporting the position IS the feature here -- a bare "there is an Inf"
  # over an 866k-column matrix is not actionable
  err <- tryCatch(calc_clocks(DNAm, "Hannum"), error = function(e) e)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_true(grepl(rownames(DNAm)[3], msg, fixed = TRUE))
  expect_true(grepl(b$panel[5], msg, fixed = TRUE))
})

test_that("the kernel bails on Inf and reports nothing else usable", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[3, b$panel[5]] <- Inf

  scan <- col_stats(DNAm[, b$panel[1:8], drop = FALSE])
  # position is 1-based and relative to the block handed to the kernel
  expect_equal(scan$inf_at, c(3L, 5L))
  # stats must not be consumed when inf_at is set
  expect_null(scan$stats)
})

test_that("each range flag warns on its own, and both can fire", {
  b <- gate_betas()

  low <- b$DNAm
  low[4, b$panel[7]] <- -0.2
  expect_equal(length(warnings_of(calc_clocks(low, "Hannum"))), 1L)

  high <- b$DNAm
  high[4, b$panel[7]] <- 1.4
  expect_equal(length(warnings_of(calc_clocks(high, "Hannum"))), 1L)

  both <- b$DNAm
  both[4, b$panel[7]] <- -0.2
  both[5, b$panel[9]] <- 1.4
  expect_equal(length(warnings_of(calc_clocks(both, "Hannum"))), 2L)

  # an M-value matrix (betas through log2(p / (1 - p))) spans both sides
  m <- log2(b$DNAm / (1 - b$DNAm))
  expect_equal(length(warnings_of(res <- calc_clocks(m, "Hannum"))), 2L)
  # warnings, not a refusal -- the caller still gets their scores
  expect_equal(nrow(res$scores), nrow(m))
})

test_that("ordinary betas pass both gates in silence", {
  b <- gate_betas()
  expect_no_warning(calc_clocks(b$DNAm, "Hannum"))

  # NAs are missing, not bad: they fill and say nothing about range
  with_na <- b$DNAm
  with_na[1:3, b$panel[1]] <- NA
  expect_no_warning(filled <- calc_clocks(with_na, "Hannum"))
  expect_equal(
    filled$coverage$per_clock[["Hannum"]]$score_imputed_partial,
    3L
  )
})

test_that("an all-missing column classifies rather than erroring", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[, b$panel[2]] <- NA

  # no observations means no cohort mean exists, which is an ordinary case,
  # not a failure: the column leaves usable_cols and the clock's vendored ref
  # or the drop policy takes over. The 0/0 that would come of averaging it is
  # avoided by classification, so nothing here should warn or stop.
  expect_no_warning(res <- calc_clocks(DNAm, "Hannum"))
  expect_true(all(is.finite(res$scores[, "Hannum"])))

  mna <- scan_missing_cpgs(DNAm, b$panel)
  expect_true(b$panel[2] %in% mna$all_na_cols)
  expect_false(b$panel[2] %in% mna$usable_cols)
  expect_false(b$panel[2] %in% names(mna$col_mean))
  expect_true(all(is.finite(mna$col_mean)))
})

test_that("col_stats sums and counts observed entries in one sweep", {
  b <- gate_betas(n = 10L)
  DNAm <- b$DNAm
  DNAm[1:4, b$panel[1]] <- NA
  DNAm[, b$panel[2]] <- NA

  scan <- col_stats(DNAm[, b$panel[1:3], drop = FALSE])
  expect_null(scan$inf_at)
  expect_false(scan$any_lt0)
  expect_false(scan$any_gt1)

  st <- scan$stats
  expect_equal(rownames(st), c("sum", "n_obs"))
  val <- function(row, col) unname(st[row, col])

  # partial NA: observed count drops, the mean comes off the same sweep
  expect_equal(val("n_obs", 1), 6)
  expect_equal(
    val("sum", 1) / val("n_obs", 1),
    mean(DNAm[5:10, b$panel[1]])
  )

  # all NA: nothing observed, so no mean is defined for it
  expect_equal(val("n_obs", 2), 0)
  expect_equal(val("sum", 2), 0)

  # untouched column
  expect_equal(val("n_obs", 3), 10)
  expect_equal(val("sum", 3), sum(DNAm[, b$panel[3]]))
})
