# report(): DNAm-input QC + score QC. Assert the structured mc_report it
# produces, not the PDF pixels or message wording (see "Test altitude").

# linear, bundled clocks that need no covariates -- keeps the fixtures small
REPORT_CLOCKS <- c("Horvath1", "Hannum")

# quiet the file-written note and QC warnings; we assert on the returned object.
# default the output to a session tempfile so tests never write into the package
# dir (report()'s real default is the working directory).
quiet_report <- function(...) {
  args <- list(...)
  if (is.null(args[["file"]])) {
    args[["file"]] <- tempfile(fileext = ".pdf")
  }
  suppressMessages(suppressWarnings(do.call(report, args)))
}

test_that("report(DNAm) builds a DNAm QC section over a clean matrix", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 6L)
  f <- withr::local_tempfile(fileext = ".pdf")

  rep <- quiet_report(sim$DNAm, clocks = REPORT_CLOCKS, file = f, ask = FALSE)

  expect_s3_class(rep, "mc_report")
  expect_true(rep$meta$has_dnam)
  expect_false(rep$meta$has_score)

  d <- rep$dnam
  expect_true(d$format$ok)
  expect_identical(d$format$n_probes, ncol(sim$DNAm))
  # every scoring panel is fully present -> coverage 1 for each assessed clock
  expect_true(all(d$coverage$coverage == 1))
  expect_setequal(intersect(REPORT_CLOCKS, d$coverage$clock_id), REPORT_CLOCKS)
  expect_true(file.exists(f) && file.info(f)$size > 0)
})

test_that("report(DNAm) counts out-of-range beta values and their probes", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 6L)
  DNAm <- sim$DNAm
  DNAm[2L, 3L] <- 1.5
  DNAm[4L, 5L] <- -0.2

  b <- quiet_report(DNAm, clocks = REPORT_CLOCKS, ask = FALSE)$dnam$beta
  expect_identical(b$n_above1, 1L)
  expect_identical(b$n_below0, 1L)
  expect_identical(b$n_out_of_range, 2L)
  expect_setequal(b$out_of_range_probes, colnames(DNAm)[c(3L, 5L)])
})

test_that("report(DNAm) separates partial-NA from fully-missing probes", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 6L)
  DNAm <- sim$DNAm
  DNAm[1:2, 4L] <- NA # partial
  DNAm[, 7L] <- NA # fully missing

  ms <- quiet_report(DNAm, clocks = REPORT_CLOCKS, ask = FALSE)$dnam$missing
  expect_true(colnames(DNAm)[4L] %in% ms$partial_probes)
  expect_true(colnames(DNAm)[7L] %in% ms$full_na_probes)
  expect_false(colnames(DNAm)[7L] %in% ms$partial_probes)
})

test_that("report(DNAm) reports absent probes as reduced coverage", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 6L)
  drop <- clock_scoring_cpgs("Hannum")[1:5]
  DNAm <- sim$DNAm[, setdiff(colnames(sim$DNAm), drop), drop = FALSE]

  d <- quiet_report(DNAm, clocks = REPORT_CLOCKS, ask = FALSE)$dnam
  expect_identical(d$missing$n_absent_from_matrix, 5L)
  hannum <- d$coverage[d$coverage$clock_id == "Hannum", ]
  expect_identical(hannum$absent, 5L)
  expect_true(hannum$coverage < 1)
})

test_that("detect_array uses probe count and EPICv2 replicate suffixes", {
  # EPICv2 replicate suffixes win regardless of count
  v2_cols <- c(paste0("cg", 1:60), paste0("cg", 1:20, "_TC21"))
  m_v2 <- matrix(0.5, nrow = 1L, ncol = length(v2_cols), dimnames = list("s1", v2_cols))
  expect_identical(detect_array(m_v2)$guess, "EPICv2")

  # a full 450K-sized cg set classifies as 450K
  cols_450 <- paste0("cg", sprintf("%08d", seq_len(MC_ARRAY_SIZES[["450K"]])))
  m_450 <- matrix(0.5, nrow = 1L, ncol = length(cols_450), dimnames = list("s1", cols_450))
  expect_identical(detect_array(m_450)$guess, "450K")

  # a small subset is not claimed as any full array
  expect_identical(detect_array(m_v2[, 1:60, drop = FALSE])$guess, "subset/unknown")
})

test_that("report_check_dnam flags a transposed matrix and bad structure", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 4L)
  expect_true(report_check_dnam(t(sim$DNAm))$transposed)
  expect_false(report_check_dnam(sim$DNAm)$transposed)

  no_rows <- sim$DNAm
  rownames(no_rows) <- NULL
  expect_error(report_check_dnam(no_rows))
  expect_error(report_check_dnam(as.data.frame(sim$DNAm)))
})

test_that("report(DNAm) builds a per-sample QC table and an overall verdict", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 8L)
  rep <- quiet_report(sim$DNAm, clocks = REPORT_CLOCKS, ask = FALSE)

  sm <- rep$dnam$samples
  expect_identical(nrow(sm$table), nrow(sim$DNAm))
  expect_true(all(c("na_frac", "mean_beta", "mid_frac", "flags") %in% names(sm$table)))

  vd <- rep$meta$verdict
  expect_true(vd$overall %in% c("PASS", "WARN", "FAIL"))
  expect_true(all(c("format", "beta", "coverage", "samples") %in% names(vd$sections)))
})

test_that("a minority of low-coverage clocks is partial (WARN), not FAIL", {
  clocks <- c("Horvath1", "Hannum", "PhenoAge", "DunedinPACE")
  sim <- sim_DNAm(clocks, n = 6L)
  # strip every Hannum probe -> that one clock is ~0 coverage, the rest are full
  D <- sim$DNAm[, setdiff(colnames(sim$DNAm), clock_scoring_cpgs("Hannum")), drop = FALSE]
  v <- quiet_report(D, clocks = clocks, ask = FALSE)$meta$verdict
  expect_identical(v$sections[["coverage"]], "WARN")
})

test_that("report(DNAm) detects a duplicated sample", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 6L)
  DNAm <- sim$DNAm
  DNAm[5L, ] <- DNAm[4L, ] # exact duplicate

  dups <- quiet_report(DNAm, clocks = REPORT_CLOCKS, ask = FALSE)$dnam$samples$duplicates
  expect_identical(nrow(dups), 1L)
  expect_setequal(c(dups$a, dups$b), rownames(DNAm)[c(4L, 5L)])
})

test_that("a cohort-wide non-bimodal distribution is flagged once, not per sample", {
  # uniform betas (sim_DNAm) are genuinely not bimodal: cohort flag fires, but
  # no single sample is a cohort outlier of the others
  sim <- sim_DNAm(REPORT_CLOCKS, n = 8L)
  sm <- quiet_report(sim$DNAm, clocks = REPORT_CLOCKS, ask = FALSE)$dnam$samples
  expect_false(sm$cohort_bimodal_ok)
  # the cohort signal must not collapse into a per-sample flag on *every* sample
  # (the original absolute-threshold bug); relative MAD outliers may be a few
  expect_true(sum(grepl("distro-outlier", sm$table$flags)) < nrow(sm$table))
})

test_that("a genuinely bimodal cohort is not flagged as non-bimodal", {
  # real methylation is bimodal with substantial intermediate CpGs; the check
  # must key on the density trough, not the middle-band fraction
  cpgs <- clock_cpgs(resolve_clocks_sequence(resolve_clocks(
    c("Horvath1", "Hannum", "PhenoAge")
  )))
  bimodal <- function(n) {
    u <- runif(n)
    ifelse(u < 0.45, rbeta(n, 2, 18), ifelse(u < 0.9, rbeta(n, 18, 2), runif(n, 0.3, 0.7)))
  }
  m <- matrix(
    0,
    8L,
    length(cpgs),
    dimnames = list(paste0("s", seq_len(8L)), cpgs)
  )
  for (i in seq_len(8L)) m[i, ] <- bimodal(length(cpgs))

  sm <- quiet_report(m, clocks = c("Horvath1", "Hannum", "PhenoAge"), ask = FALSE)$dnam$samples
  expect_true(sm$cohort_bimodal_ok)
  expect_false(any(grepl("distro-outlier", sm$table$flags)))
})

test_that("report(result) flags samples with NA scores", {
  # a probe absent for a sample cannot make a whole clock NA on its own, so
  # drive an NA score via a fully-missing covariate for one sample
  sim <- sim_DNAm("DNAmGrip_wAge", n = 6L, Female = TRUE)
  sim$pheno$Age <- c(NA, 50, 55, 60, 65, 70)
  res <- suppressWarnings(
    calc_clocks(sim$DNAm, "DNAmGrip_wAge", pheno = sim$pheno)
  )
  na <- quiet_report(res)$score$na_samples
  expect_true(nrow(na) >= 1L)
  expect_true(sim$pheno$ID[1] %in% na$sample_id)
})

test_that("report(result) computes a clock correlation matrix for >1 clock", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 10L)
  res <- calc_clocks(sim$DNAm, REPORT_CLOCKS, pheno = sim$pheno)
  cm <- quiet_report(res)$score$correlations
  expect_true(is.matrix(cm))
  expect_identical(dim(cm), c(ncol(res$scores), ncol(res$scores)))
})

test_that("out-of-range flags fire when the reference carries expect_lo/expect_hi", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 8L)
  res <- calc_clocks(sim$DNAm, REPORT_CLOCKS, pheno = sim$pheno)
  # an impossible window forces every finite score out of range
  ref <- data.frame(
    clock_id = colnames(res$scores),
    expect_lo = 1e6,
    expect_hi = 1e6 + 1,
    stringsAsFactors = FALSE
  )
  rf <- quiet_report(result = res, score_reference = ref)$score$range_flags
  expect_s3_class(rf, "data.frame")
  expect_true(all(rf$n_out_of_range > 0L))
})

test_that("report(result) builds a score section, one row per score column", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 8L)
  res <- calc_clocks(sim$DNAm, REPORT_CLOCKS, pheno = sim$pheno)
  f <- withr::local_tempfile(fileext = ".pdf")

  rep <- quiet_report(res, file = f)
  expect_true(rep$meta$has_score)
  expect_false(rep$meta$has_dnam)
  expect_identical(nrow(rep$score$summary), ncol(res$scores))
  expect_setequal(rep$score$summary$clock_id, colnames(res$scores))
  expect_null(rep$score$flags)
  expect_true(file.exists(f) && file.info(f)$size > 0)
})

test_that("report scores against a supplied reference and flags deviations", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 8L)
  res <- calc_clocks(sim$DNAm, REPORT_CLOCKS, pheno = sim$pheno)
  ref <- data.frame(
    clock_id = colnames(res$scores),
    ref_mean = 0,
    ref_sd = 1,
    stringsAsFactors = FALSE
  )

  rep <- quiet_report(result = res, score_reference = ref)
  expect_s3_class(rep$score$flags, "data.frame")
  # age clocks (mean ~ decades) vs ref_mean 0, sd 1 -> large |z|, all flagged
  expect_identical(nrow(rep$score$flags), ncol(res$scores))
})

test_that("report(result, pheno) checks age associations against the reference", {
  skip_if(is.null(mc_clock_reference()))
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 30L, Age = TRUE)
  res <- calc_clocks(sim$DNAm, c("Horvath1", "Hannum"), pheno = sim$pheno)
  age <- sim$pheno$Age
  # overwrite scores with known age relationships (deterministic, unseeded-safe):
  # Horvath1 anti-correlated (broken age clock), Hannum tracks age (working)
  res$scores[, "Horvath1"] <- -age + rnorm(length(age), 0, 2)
  res$scores[, "Hannum"] <- age + rnorm(length(age), 0, 2)

  a <- quiet_report(result = res, pheno = sim$pheno, pheno_id = "ID")$score$associations
  expect_s3_class(a$table, "data.frame")
  expect_true(all(c("clock", "obs_age_r", "exp_age_r", "exp_lo", "exp_hi", "outside")
  %in% names(a$table)))
  hor <- a$table[a$table$clock == "Horvath1", ]
  han <- a$table[a$table$clock == "Hannum", ]
  expect_true(hor$outside && hor$wrong_sign) # broken age clock is caught
  expect_false(han$outside) # working age clock is within range
})

test_that("age associations are absent without Age in pheno", {
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 8L)
  res <- calc_clocks(sim$DNAm, c("Horvath1", "Hannum"), pheno = sim$pheno)
  expect_null(quiet_report(result = res)$score$associations)
})

test_that("report(DNAm, result) produces both sections", {
  sim <- sim_DNAm(REPORT_CLOCKS, n = 6L)
  res <- calc_clocks(sim$DNAm, REPORT_CLOCKS, pheno = sim$pheno)
  f <- withr::local_tempfile(fileext = ".pdf")

  rep <- quiet_report(sim$DNAm, res, clocks = REPORT_CLOCKS, file = f, ask = FALSE)
  expect_true(rep$meta$has_dnam)
  expect_true(rep$meta$has_score)
})

test_that("report() with no usable input errors", {
  expect_error(report())
})

test_that("a numeric data.frame DNAm is coerced; non-numeric columns abort", {
  m <- random_betas(clock_scoring_cpgs("Hannum"), n = 4L)
  df <- as.data.frame(m)
  # report() and calc_clocks() both coerce and run
  expect_no_error(quiet_report(df, clocks = "Hannum", ask = FALSE))
  expect_no_error(suppressMessages(calc_clocks(df, "Hannum")))
  # a non-numeric column cannot be a beta matrix -> refuse
  bad <- df
  bad$label <- "x"
  expect_error(suppressMessages(calc_clocks(bad, "Hannum")))
})

test_that("report survives degenerate inputs (all-NA aborts cleanly; n=1 runs)", {
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 4L)
  allna <- sim$DNAm
  allna[] <- NA_real_
  # clean diagnostic, not the cli "Cannot pluralize" crash (distinct failure mode)
  err <- tryCatch(
    quiet_report(allna, clocks = c("Horvath1", "Hannum"), ask = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("pluralize", err))
  expect_match(err, "observed CpGs")

  # a single sample must not crash the renderer (hist / heatmap guards)
  s1 <- sim_DNAm(c("Horvath1", "Hannum"), n = 1L, Age = TRUE)
  r1 <- suppressWarnings(calc_clocks(s1$DNAm, c("Horvath1", "Hannum"), pheno = s1$pheno))
  f <- withr::local_tempfile(fileext = ".pdf")
  expect_no_error(quiet_report(result = r1, pheno = s1$pheno, file = f))
})
