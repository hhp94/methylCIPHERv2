# value gates: overflow stops, Inf and out-of-[0,1] warn (via calc_clocks
# where possible)

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

test_that("an infinite value scores as missing, and warns that it did", {
  b <- gate_betas()

  # no NA anywhere (anyNA does not see Inf)
  inf_only <- b$DNAm
  inf_only[1, b$panel[1]] <- Inf
  expect_equal(length(warnings_of(res <- calc_clocks(inf_only, "Hannum"))), 1L)
  expect_true(all(is.finite(res$scores[, "Hannum"])))
  # imputed exactly like an NA in the same cell would have been
  expect_equal(res$coverage$per_clock[["Hannum"]]$score_imputed_partial, 1L)

  neg <- b$DNAm
  neg[3, b$panel[5]] <- -Inf
  expect_equal(length(warnings_of(calc_clocks(neg, "Hannum"))), 1L)
})

test_that("an Inf and an NA in the same cell score identically", {
  b <- gate_betas()
  as_inf <- as_na <- b$DNAm
  as_inf[3, b$panel[5]] <- Inf
  as_na[3, b$panel[5]] <- NA

  inf_res <- suppressWarnings(calc_clocks(as_inf, "Hannum"))
  na_res <- calc_clocks(as_na, "Hannum")

  expect_equal(inf_res$scores, na_res$scores)
  expect_equal(inf_res$coverage$sample_miss, na_res$coverage$sample_miss)
})

test_that("a wholly infinite column classifies absent, like an all-NA one", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[, b$panel[2]] <- Inf

  mna <- suppressWarnings(scan_missing_cpgs(DNAm, b$panel, b$panel))
  expect_true(b$panel[2] %in% mna$all_na_cols)
  expect_false(b$panel[2] %in% mna$usable_cols)
  expect_true(all(is.finite(mna$col_mean)))
})

test_that("the kernel counts Inf as missing, not as observed", {
  b <- gate_betas(n = 6L)
  DNAm <- b$DNAm
  DNAm[3, b$panel[5]] <- Inf

  scan <- col_stats(DNAm[, b$panel[1:8], drop = FALSE])
  expect_true(scan$any_inf)
  expect_null(scan$overflow_col)
  # the Inf row observed 7 of the 8 columns, and the mean skips it
  expect_equal(scan$row_obs, c(8L, 8L, 7L, 8L, 8L, 8L))
  expect_equal(unname(scan$stats["n_obs", 5]), 5)
  expect_equal(
    unname(scan$stats["sum", 5]),
    sum(DNAm[-3, b$panel[5]])
  )
  # an Inf is missing, so it says nothing about range
  expect_false(scan$any_gt1)
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

test_that("per-sample fill counts land on the samples that were filled", {
  b <- gate_betas(n = 6L)
  DNAm <- b$DNAm
  # 2 CpGs filled for sample 1, 1 for sample 4, none for anyone else
  DNAm[1, b$panel[1:2]] <- NA
  DNAm[4, b$panel[1]] <- NA

  res <- calc_clocks(DNAm, "Hannum")
  miss <- res$coverage$sample_miss$score[, "Hannum"]

  expect_equal(unname(miss), c(2L, 0L, 0L, 1L, 0L, 0L))
  expect_equal(names(miss), rownames(DNAm))
  # the per-clock record aggregates the same counts
  expect_equal(res$coverage$per_clock[["Hannum"]]$score_imputed_partial, 3L)
})

test_that("an all-missing column classifies rather than erroring", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[, b$panel[2]] <- NA

  # all-NA column is ordinary -- classified, not averaged
  expect_no_warning(res <- calc_clocks(DNAm, "Hannum"))
  expect_true(all(is.finite(res$scores[, "Hannum"])))

  mna <- scan_missing_cpgs(DNAm, b$panel, b$panel)
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
  expect_equal(scan$overflow_col, 4L)
  expect_null(scan$stats)

  # naming the column is the feature
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

test_that("a dead row is judged on the scoring panel, not the norm one", {
  spec <- mc_spec("DunedinPACE")
  score <- panels_union(spec$panels, "score")
  norm_only <- setdiff(panels_union(spec$panels), score)
  expect_true(length(norm_only) > 0L)

  DNAm <- random_betas(panels_union(spec$panels), n = 4L)

  # a full normalization background does not make a sample scoreable
  dead <- DNAm
  dead[1, score] <- NA
  expect_true(all(is.finite(dead[1, norm_only])))
  expect_error(mc_cohort(dead, spec))

  # the converse is not fatal -- a thin background only warns
  thin <- DNAm
  thin[1, norm_only] <- NA
  expect_no_error(mc_cohort(thin, spec))
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
  expect_null(scan$overflow_col)
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

# cols: scan a subset of DNAm's columns without materializing a slice

test_that("scanning by cols agrees with scanning a pre-subset matrix", {
  b <- gate_betas(n = 8L)
  DNAm <- b$DNAm
  DNAm[1:3, b$panel[1]] <- NA
  DNAm[, b$panel[2]] <- NA

  # deliberately not in column order
  sel <- b$panel[c(5, 1, 9, 2)]
  direct <- col_stats(DNAm, match(sel, colnames(DNAm)))
  sliced <- col_stats(DNAm[, sel, drop = FALSE])

  expect_equal(direct$stats, sliced$stats)
  expect_equal(direct$row_obs, sliced$row_obs)
  expect_equal(direct$any_lt0, sliced$any_lt0)
  expect_equal(direct$any_gt1, sliced$any_gt1)
  expect_null(direct$overflow_col)

  # results are ordered by cols, not by position in DNAm
  expect_equal(unname(direct$stats["n_obs", 2]), 5)
  expect_equal(unname(direct$stats["sum", 2]), sum(DNAm[4:8, b$panel[1]]))
  expect_equal(unname(direct$stats["n_obs", 4]), 0)
})

test_that("cols narrows the sweep -- unlisted columns are not scanned", {
  b <- gate_betas(n = 6L)
  DNAm <- b$DNAm
  DNAm[2, b$panel[1:4]] <- NA
  # out-of-range value the narrow scan must not see
  DNAm[, b$panel[20]] <- -0.5

  narrow <- col_stats(DNAm, match(b$panel[1:4], colnames(DNAm)))
  expect_equal(narrow$row_obs, c(4L, 0L, 4L, 4L, 4L, 4L))
  expect_false(narrow$any_lt0)

  # widen over the same matrix and the flag fires
  wide <- col_stats(DNAm, match(b$panel[c(1:4, 20)], colnames(DNAm)))
  expect_true(wide$any_lt0)
  expect_equal(wide$row_obs, c(5L, 1L, 5L, 5L, 5L, 5L))
})

test_that("overflow_col names a position within cols, not a column of DNAm", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[, b$panel[5]] <- 1e308

  # panel[5] is the 2nd entry of cols
  scan <- col_stats(DNAm, match(b$panel[c(1, 5, 9)], colnames(DNAm)))
  expect_equal(scan$overflow_col, 2L)
  expect_null(scan$stats)
  expect_null(scan$row_obs)
})

test_that("cols defaults to every column of DNAm", {
  b <- gate_betas(n = 5L)
  sub <- b$DNAm[, b$panel[1:6], drop = FALSE]
  expect_equal(col_stats(sub), col_stats(sub, seq_len(ncol(sub))))
})

test_that("cols must be in-range column indices", {
  b <- gate_betas(n = 4L)
  DNAm <- b$DNAm
  expect_error(col_stats(DNAm, 0L))
  expect_error(col_stats(DNAm, -1L))
  expect_error(col_stats(DNAm, ncol(DNAm) + 1L))
  expect_error(col_stats(DNAm, NA_integer_))
})

test_that("an empty cols scans nothing and observes nothing", {
  b <- gate_betas(n = 4L)
  scan <- col_stats(b$DNAm, integer(0))
  expect_equal(dim(scan$stats), c(2L, 0L))
  expect_equal(scan$row_obs, integer(4))
  expect_null(scan$overflow_col)
})
