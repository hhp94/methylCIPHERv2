# value gates: overflow stops. Inf and out-of-[0,1] warn.
# col_stats() kernel blocks below are gated. they pin index base, mask width, overflow bail.

GATE_PANEL <- mc_spec("Hannum")$needed_union

gate_betas <- function(n = 8L) {
  list(DNAm = random_betas(GATE_PANEL, n = n), panel = GATE_PANEL)
}

test_that("an infinite value scores as missing, and warns that it did", {
  b <- gate_betas()

  # no NA anywhere (anyNA does not see Inf)
  as_inf <- as_na <- b$DNAm
  as_inf[3, b$panel[5]] <- Inf
  as_na[3, b$panel[5]] <- NA

  expect_equal(
    length(capture_warnings(inf_res <- calc_clocks(as_inf, "Hannum"))),
    1L
  )
  expect_true(all(is.finite(inf_res$scores[, "Hannum"])))

  # imputed exactly like an NA in the same cell would have been
  na_res <- calc_clocks(as_na, "Hannum")
  expect_equal(inf_res$scores, na_res$scores)
  expect_equal(inf_res$coverage$sample_miss, na_res$coverage$sample_miss)
  expect_equal(
    inf_res$coverage$per_clock[[1]][["Hannum"]]$score_imputed_partial,
    1L
  )

  neg <- b$DNAm
  neg[3, b$panel[5]] <- -Inf
  expect_equal(length(capture_warnings(calc_clocks(neg, "Hannum"))), 1L)
})

test_that("each range flag warns on its own, and both can fire", {
  b <- gate_betas()

  low <- b$DNAm
  low[4, b$panel[7]] <- -0.2
  expect_equal(length(capture_warnings(calc_clocks(low, "Hannum"))), 1L)

  high <- b$DNAm
  high[4, b$panel[7]] <- 1.4
  expect_equal(length(capture_warnings(calc_clocks(high, "Hannum"))), 1L)

  # an M-value matrix (betas through log2(p / (1 - p))) spans both sides
  m <- log2(b$DNAm / (1 - b$DNAm))
  expect_equal(length(capture_warnings(res <- calc_clocks(m, "Hannum"))), 2L)
  # warnings, not a refusal -- the caller still gets their scores
  expect_equal(nrow(res$scores), nrow(m))
})

test_that("ordinary betas pass both gates in silence, and NA fills are counted", {
  b <- gate_betas(n = 6L)
  expect_no_warning(calc_clocks(b$DNAm, "Hannum"))

  DNAm <- b$DNAm
  # nas are missing, not bad: 2 CpGs filled for sample 1, 1 for sample 4
  DNAm[1, b$panel[1:2]] <- NA
  DNAm[4, b$panel[1]] <- NA
  expect_no_warning(res <- calc_clocks(DNAm, "Hannum"))

  miss <- res$coverage$sample_miss$score[, "Hannum"]
  expect_equal(unname(miss), c(2L, 0L, 0L, 1L, 0L, 0L))
  expect_equal(names(miss), rownames(DNAm))
  # the record counts the other axis: 2 distinct CpGs, not the 3 filled cells
  expect_equal(
    res$coverage$per_clock[[1]][["Hannum"]]$score_imputed_partial,
    2L
  )
})

test_that("a column that overflows its own sum stops, and names the column", {
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[, b$panel[4]] <- 1e308

  # naming the column is the feature, so the content is asserted, not the wording
  err <- tryCatch(calc_clocks(DNAm, "Hannum"), error = function(e) e)
  expect_s3_class(err, "error")
  expect_true(grepl(b$panel[4], conditionMessage(err), fixed = TRUE))
})

test_that("a sample with nothing on the scoring panel stops, off-panel or not", {
  b <- gate_betas()
  extra <- paste0("cg_offpanel_", seq_len(20L))
  DNAm <- cbind(b$DNAm, random_betas(extra, n = nrow(b$DNAm)))

  # every panel CpG is partial-NA cohort-wide. fill cannot rescue this row.
  DNAm[1, b$panel] <- NA
  expect_error(calc_clocks(DNAm, "Hannum"))
})

test_that("the moments span columns the column stats and gates never see", {
  skip_on_cran()
  # zhang2019EN z-scores each sample over the whole matrix.
  panel <- suppressMessages(panels_union(mc_spec("Zhang2019EN")$panels))
  DNAm <- random_betas(panel, n = 4L)
  off <- paste0("cg_offpanel_", seq_len(50L))
  wide <- cbind(DNAm, random_betas(off, n = 4L))
  wide[, off] <- 0.9

  # off-panel columns are in the moments, so every sample's score moves ...
  suppressMessages({
    with_off <- calc_clocks(wide, "Zhang2019EN")$scores
    panel_only <- calc_clocks(DNAm, "Zhang2019EN")$scores
  })
  expect_true(all(abs(with_off - panel_only) > 1e-6))

  # ... but not in the column stats, so no value gate ever looks at them
  wide[1, off[[1]]] <- -0.5
  expect_no_warning(suppressMessages(calc_clocks(wide, "Zhang2019EN")))

  # off-panel Inf is skipped by the moments, not poison.
  as_inf <- as_na <- wide
  as_inf[2, off[[2]]] <- Inf
  as_na[2, off[[2]]] <- NA
  suppressWarnings(suppressMessages({
    inf_res <- calc_clocks(as_inf, "Zhang2019EN")
    na_res <- calc_clocks(as_na, "Zhang2019EN")
  }))
  expect_equal(inf_res$scores, na_res$scores)
  expect_true(all(is.finite(inf_res$scores[, "Zhang2019EN"])))
})

test_that("an unscored sample is NA, which is not a non-finite score", {
  skip_on_cran()
  expect_no_warning(check_score_values(list(Hannum = matrix(c(1, NA_real_)))))
  expect_warning(check_score_values(list(Hannum = matrix(c(1, NaN)))))
  expect_warning(check_score_values(list(Hannum = matrix(c(1, Inf)))))
})

test_that("an all-NA or wholly infinite column classifies rather than erroring", {
  skip_on_cran()
  b <- gate_betas()

  # all-NA column is ordinary -- classified, not averaged
  na_col <- b$DNAm
  na_col[, b$panel[2]] <- NA
  expect_no_warning(res <- calc_clocks(na_col, "Hannum"))
  expect_true(all(is.finite(res$scores[, "Hannum"])))

  # a wholly infinite column lands in the same bucket
  inf_col <- b$DNAm
  inf_col[, b$panel[2]] <- Inf
  mna <- suppressWarnings(scan_missing_cpgs(inf_col, b$panel, b$panel))
  expect_true(b$panel[2] %in% mna$all_na_cols)
  expect_false(b$panel[2] %in% mna$usable_cols)
  expect_false(b$panel[2] %in% names(mna$col_mean))
  expect_true(all(is.finite(mna$col_mean)))
})

test_that("a dead row is judged on the scoring panel, not the norm one", {
  skip_on_cran()
  spec <- mc_spec("DunedinPACE")
  score <- panels_union(spec$panels, "score")
  norm_only <- setdiff(panels_union(spec$panels), score)
  expect_true(length(norm_only) > 0L)

  DNAm <- random_betas(panels_union(spec$panels), n = 4L)

  # a full normalization background does not make a sample scoreable
  dead <- DNAm
  dead[1, score] <- NA
  expect_error(mc_cohort(dead, spec))

  # the converse is not fatal -- a thin background only warns
  thin <- DNAm
  thin[1, norm_only] <- NA
  expect_no_error(mc_cohort(thin, spec))
})

test_that("the kernel counts Inf as missing, not as observed", {
  skip_on_cran()
  b <- gate_betas(n = 6L)
  DNAm <- b$DNAm
  DNAm[3, b$panel[5]] <- Inf

  scan <- col_stats(DNAm[, b$panel[1:8], drop = FALSE])
  expect_true(scan$any_inf)
  expect_null(scan$overflow_col)
  # the Inf row observed 7 of the 8 columns, and the column count skips it
  expect_equal(scan$row_obs, c(8L, 8L, 7L, 8L, 8L, 8L))
  expect_equal(unname(scan$stats["n_obs", 5]), 5)
  # an Inf is missing, so it says nothing about range
  expect_equal(scan$max_val, 1)
  expect_true(is.na(scan$max_col))
})

# cols: scan a subset of DNAm's columns without materializing a slice

test_that("cols selects positionally, and results are ordered by cols", {
  skip_on_cran()
  b <- gate_betas(n = 8L)
  DNAm <- b$DNAm
  DNAm[1:3, b$panel[1]] <- NA
  DNAm[, b$panel[2]] <- NA
  # out-of-range value only the wider sweep may see
  DNAm[, b$panel[20]] <- -0.5

  # deliberately not in column order
  sel <- b$panel[c(5, 1, 9, 2)]
  direct <- col_stats(DNAm, match(sel, colnames(DNAm)))
  sliced <- col_stats(DNAm[, sel, drop = FALSE])
  expect_equal(direct$stats, sliced$stats)
  expect_equal(direct$row_obs, sliced$row_obs)
  expect_null(direct$overflow_col)

  # results are ordered by cols, not by position in DNAm
  expect_equal(unname(direct$stats["n_obs", 2]), 5)
  expect_equal(unname(direct$stats["sum", 2]), sum(DNAm[4:8, b$panel[1]]))
  expect_equal(unname(direct$stats["n_obs", 4]), 0)

  # unlisted columns are not scanned, so the range moves only when they join
  expect_equal(direct$min_val, 0)
  wide <- col_stats(DNAm, match(c(sel, b$panel[20]), colnames(DNAm)))
  expect_equal(wide$min_val, -0.5)
  # position within cols, like overflow_col -- panel[20] is the 5th entry
  expect_equal(wide$min_col, 5L)

  # cols defaults to every column of the matrix it is handed
  sub <- DNAm[, sel, drop = FALSE]
  expect_equal(col_stats(sub), col_stats(sub, seq_len(ncol(sub))))
})

test_that("the overflow bail nulls every output, at a position within cols", {
  skip_on_cran()
  b <- gate_betas()
  DNAm <- b$DNAm
  DNAm[, b$panel[5]] <- 1e308

  # panel[5] is the 2nd entry of cols
  sel <- match(b$panel[c(1, 5, 9)], colnames(DNAm))
  scan <- col_stats(DNAm, sel, list(a = seq_len(20L)))
  expect_equal(scan$overflow_col, 2L)
  expect_null(scan$stats)
  expect_null(scan$row_obs)
  expect_null(scan$row_moment_obs)

  # cols is bounds-checked ahead of the kernel
  expect_error(col_stats(b$DNAm, 0L))
  expect_error(col_stats(b$DNAm, ncol(b$DNAm) + 1L))
  expect_error(col_stats(b$DNAm, NA_integer_))

  # an empty cols scans nothing and observes nothing
  empty <- col_stats(b$DNAm, integer(0))
  expect_equal(dim(empty$stats), c(2L, 0L))
  expect_equal(empty$row_obs, integer(8))
})

test_that("R and the kernel agree on the mask width, and the catalog fits it", {
  skip_on_cran()
  b <- gate_betas(n = 3L)
  one <- lapply(seq_len(MAX_MOMENT_SETS), function(i) i)
  expect_no_error(col_stats(b$DNAm, NULL, one))
  # one past the width is a stop, not a silently truncated mask
  expect_error(col_stats(b$DNAm, NULL, c(one, list(1L))))
  expect_error(check_moment_sets(c(one, list(1L)), 40L))

  # validator keeps what the kernel accepts and rejects what would be fatal there.
  expect_null(check_moment_sets(NULL, 40L))
  expect_equal(check_moment_sets(list(a = c(1, 2, 3)), 40L)[["a"]], 1:3)
  # a ref meeting no measured column is a data fact, not a usage error
  expect_equal(check_moment_sets(list(a = integer(0)), 40L)[["a"]], integer(0))
  expect_error(check_moment_sets(list(a = "x"), 40L))
  expect_error(check_moment_sets(list(), 40L))
  expect_error(check_moment_sets(list(a = 0L), 40L))
  expect_error(check_moment_sets(list(a = 41L), 40L))

  # mask width is a ceiling on domains. census keys only.
  keys <- vapply(
    names(mc_catalog),
    function(id) clock_moment_key(id) %||% NA_character_,
    ""
  )
  expect_lte(length(unique(stats::na.omit(keys))), MAX_MOMENT_SETS)
})

# moment domains: each set carries its own counter.

# per-row golden over a column subset, finite entries only
domain_golden <- function(DNAm, cols, f) {
  sub <- DNAm[, cols, drop = FALSE]
  unname(apply(sub, 1, function(v) f(v[is.finite(v)])))
}

# obs/mean/sd triple one set's kernel outputs must match.
expect_domain_moments <- function(got, DNAm, cols, nm) {
  expect_equal(
    got$row_moment_obs[, nm],
    as.integer(rowSums(is.finite(DNAm[, cols, drop = FALSE]))),
    info = nm
  )
  expect_equal(got$row_mean[, nm], domain_golden(DNAm, cols, mean), info = nm)
  expect_equal(
    sqrt(got$row_m2[, nm] / (got$row_moment_obs[, nm] - 1)),
    domain_golden(DNAm, cols, stats::sd),
    info = nm
  )
}

test_that("overlapping domains are each counted in full in one sweep", {
  skip_on_cran()
  b <- gate_betas(n = 6L)
  DNAm <- b$DNAm
  DNAm[1, b$panel[1]] <- NA
  DNAm[2, b$panel[12]] <- NA

  # a partial overlap and the whole matrix, all in one call
  cols <- list(
    a = b$panel[1:10],
    b = b$panel[6:15],
    all = colnames(DNAm)
  )
  idx <- lapply(cols, function(cc) match(cc, colnames(DNAm)))
  sel <- match(b$panel[1:3], colnames(DNAm))
  got <- col_stats(DNAm, sel, idx)

  for (nm in names(cols)) {
    expect_domain_moments(got, DNAm, cols[[nm]], nm)
  }
  # the moment outputs are nr x K, named by their sets, and absent without them
  expect_equal(dim(got$row_mean), c(6L, 3L))
  expect_equal(colnames(got$row_moment_obs), names(cols))
  expect_null(col_stats(DNAm, sel)$row_moment_obs)

  # a moment column is counted once. widening moments never moves the column half or row_obs.
  flat <- col_stats(DNAm, sel)
  expect_equal(got$stats, flat$stats)
  expect_equal(got$row_obs, flat$row_obs)
  disjoint <- col_stats(DNAm, match(b$panel[20:25], colnames(DNAm)), idx)
  expect_equal(disjoint$row_mean, got$row_mean)
  expect_equal(disjoint$row_moment_obs, got$row_moment_obs)
})

test_that("scan_missing_cpgs banks one moment entry per declared domain", {
  skip_on_cran()
  b <- gate_betas(n = 7L)
  DNAm <- b$DNAm
  sel <- b$panel[1:5]
  measured <- b$panel[20:30]
  ref <- c(measured, "cg_not_measured_at_all")
  DNAm[1, sel[1]] <- NA
  DNAm[2, setdiff(colnames(DNAm), sel)[1]] <- Inf

  mna <- suppressWarnings(scan_missing_cpgs(
    DNAm,
    sel,
    sel,
    moment_domains = list(full = NULL, ref = ref)
  ))
  expect_equal(names(mna$sample_moments), c("full", "ref"))
  # a NULL domain is every column
  expect_equal(
    mna$sample_moments$full$mean,
    domain_golden(DNAm, colnames(DNAm), mean)
  )
  expect_equal(
    mna$sample_moments$full$sd,
    domain_golden(DNAm, colnames(DNAm), stats::sd)
  )
  # a declared domain resolves against what was measured. unmeasured CpGs drop out.
  expect_equal(mna$sample_moments$ref$mean, domain_golden(DNAm, measured, mean))
  expect_equal(
    mna$sample_moments$ref$sd,
    domain_golden(DNAm, measured, stats::sd)
  )

  # no declared domain banks no moments at all
  expect_null(scan_missing_cpgs(b$DNAm, sel, sel)$sample_moments)
  expect_null(
    scan_missing_cpgs(b$DNAm, sel, sel, moment_domains = list())$sample_moments
  )
})

test_that("a domain too thin to describe a sample reports NA, not zero spread", {
  skip_on_cran()
  b <- gate_betas(n = 5L)
  DNAm <- b$DNAm
  thin <- b$panel[1:3]
  # sample 2 observes nothing in the domain, sample 3 observes exactly one
  DNAm[2, thin] <- NA
  DNAm[3, thin[2:3]] <- NA
  score <- b$panel[10:20]

  mna <- scan_missing_cpgs(
    DNAm,
    score,
    score,
    moment_domains = list(thin = thin, one = thin[1], none = character(0))
  )
  got <- mna$sample_moments$thin
  # n = 0 gives neither. n = 1 gives a mean but no spread
  expect_true(is.na(got$mean[2]))
  expect_true(is.na(got$sd[2]))
  expect_equal(got$mean[3], DNAm[3, thin[1]])
  expect_true(is.na(got$sd[3]))
  # the rest of the column is untouched by the guard
  expect_false(anyNA(got$sd[c(1, 4, 5)]))

  # a one-column domain can never carry a spread, and never reads as 0
  expect_true(all(is.na(mna$sample_moments$one$sd)))
  # an empty domain is a data fact, reported as NA on both moments
  expect_true(all(is.na(mna$sample_moments$none$mean)))
  expect_true(all(is.na(mna$sample_moments$none$sd)))
})

# pheno Age units gate: per row, warn only, both sides independent.

test_that("the Age units gate judges each surviving row on its own", {
  skip_on_cran()
  age_pheno <- function(age) {
    data.frame(ID = paste0("s", seq_along(age)), Age = age)
  }

  # the whole reason this is per row: a cohort statistic cannot see one bad row
  ph <- age_pheno(c(45, 52, 600, 38, 61, 47, 55, 39))
  expect_warning(flagged <- warn_age_units(ph, "ID", ph$ID))
  expect_equal(flagged, "s3")

  # both bounds fire independently on one pheno
  two <- age_pheno(c(45, 600, -38, 38))
  expect_equal(
    length(capture_warnings(flagged <- warn_age_units(two, "ID", two$ID))),
    2L
  )
  expect_equal(flagged, c("s2", "s3"))

  # -0.5/-1 are a legitimate pre-birth convention and 122 is a real maximum
  ok <- age_pheno(c(-0.5, -1, 122, 0, 45, NA, 118))
  expect_no_warning(flagged <- warn_age_units(ok, "ID", ok$ID))
  expect_equal(flagged, character(0))

  # the out-of-range row is not in this run's sample set, so it is not judged
  expect_no_warning(warn_age_units(ph, "ID", c("s1", "s2")))
})

test_that("calc_accel reaches the gate through its own pheno merge", {
  b <- gate_betas(n = 4L)
  res <- suppressMessages(calc_clocks(b$DNAm, "Hannum"))
  bad <- data.frame(ID = rownames(b$DNAm), Age = c(45, 52, 600, 38))
  expect_warning(calc_accel(res, type = "diff", data = bad))
})
