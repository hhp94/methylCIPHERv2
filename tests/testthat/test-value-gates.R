# value gates: Inf stops, out-of-[0,1] warns (via calc_clocks where possible)

gate_betas <- function(n = 8L) {
  spec <- mc_spec("Hannum")
  panel <- panels_union(spec$panels)
  list(DNAm = random_betas(panel, n = n), panel = panel)
}

# collect every warning (the two range flags are independent)
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

  # no NA anywhere (anyNA does not see Inf)
  inf_only <- b$DNAm
  inf_only[1, b$panel[1]] <- Inf
  expect_error(calc_clocks(inf_only, "Hannum"))

  # inf must still stop when an unrelated NA is present
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

  # position in the message is the feature (asserted on purpose)
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

  # all-NA column is ordinary -- classified, not averaged
  expect_no_warning(res <- calc_clocks(DNAm, "Hannum"))
  expect_true(all(is.finite(res$scores[, "Hannum"])))

  mna <- scan_missing_cpgs(DNAm, b$panel)
  expect_true(b$panel[2] %in% mna$all_na_cols)
  expect_false(b$panel[2] %in% mna$usable_cols)
  expect_false(b$panel[2] %in% names(mna$col_mean))
  expect_true(all(is.finite(mna$col_mean)))
})

test_that("a column that overflows its own sum stops, and names no sample", {
  b <- gate_betas()
  DNAm <- b$DNAm
  # overflow: column reported, row is NA (no invented position)
  DNAm[, b$panel[4]] <- 1e308

  scan <- col_stats(DNAm[, b$panel[1:8], drop = FALSE])
  expect_true(is.na(scan$inf_at[[1L]]))
  expect_equal(scan$inf_at[[2L]], 4L)
  expect_null(scan$stats)

  # naming the column is the feature, same as for Inf above
  err <- tryCatch(calc_clocks(DNAm, "Hannum"), error = function(e) e)
  expect_s3_class(err, "error")
  expect_true(grepl(b$panel[4], conditionMessage(err), fixed = TRUE))
})

test_that("a sample with nothing on the scoring panel stops, off-panel or not", {
  b <- gate_betas()
  extra <- paste0("cg_offpanel_", seq_len(20L))
  DNAm <- cbind(b$DNAm, random_betas(extra, n = nrow(b$DNAm)))

  # all panel CpGs partial-NA cohort-wide would fill this sample from others
  DNAm[1, b$panel] <- NA
  expect_error(calc_clocks(DNAm, "Hannum"))

  # off-panel observations do not make a sample scoreable
  expect_true(all(is.finite(DNAm[1, extra])))
})

test_that("col_stats counts observed entries per row as well as per column", {
  b <- gate_betas(n = 6L)
  DNAm <- b$DNAm
  DNAm[2, b$panel[1:4]] <- NA
  DNAm[3, ] <- NA

  scan <- col_stats(DNAm[, b$panel[1:4], drop = FALSE])
  expect_equal(scan$row_obs, c(4L, 0L, 0L, 4L, 4L, 4L))
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
