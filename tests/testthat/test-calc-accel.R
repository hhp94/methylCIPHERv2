# finalizer output: as.data.frame(mc_result) and calc_accel()

# grimageV1 needs Age + Female. hannum needs neither.
ACCEL_CLOCKS <- c("GrimAgeV1", "Hannum")

accel_fixture <- function(n = 12L) {
  DNAm <- random_betas(clock_cpgs(ACCEL_CLOCKS), n = n)
  pheno <- mc_pheno(
    rownames(DNAm),
    Age = stats::rnorm(n, 45, 8),
    Female = rep(c(0, 1), length.out = n)
  )
  list(
    DNAm = DNAm,
    pheno = pheno,
    res = calc_clocks(DNAm, ACCEL_CLOCKS, pheno = pheno)
  )
}

# one draw, shared by every test below that only reads the record
FX <- accel_fixture()

test_that("as.data.frame long and wide carry the same scores", {
  fx <- FX
  wide <- as.data.frame(fx$res, long = FALSE)
  long <- as.data.frame(fx$res)
  clocks <- colnames(as.matrix(fx$res))

  # id first on both shapes. one batch, so no batch column
  expect_equal(names(wide), c("ID", clocks))
  expect_equal(names(long), c("ID", "clock_id", "score"))

  # wide is as.matrix() with an id column bolted on
  expect_equal(
    unname(as.matrix(wide[, clocks, drop = FALSE])),
    unname(as.matrix(fx$res))
  )
  # long is the same cells, clock-major
  expect_equal(long$score, as.vector(as.matrix(fx$res)))
  expect_equal(nrow(long), nrow(wide) * length(clocks))
  expect_equal(long$ID, rep(fx$pheno$ID, times = length(clocks)))

  # keyed by the id column only -- no row names
  expect_equal(attr(wide, "row.names"), seq_len(nrow(wide)))
})

test_that("accel is the residual of the fitted model, default rhs or supplied", {
  fx <- FX
  scores <- as.matrix(fx$res)
  one <- calc_accel(fx$res, long = FALSE)

  # the default rhs is ~ Age, fitted per clock
  for (id in colnames(scores)) {
    d <- data.frame(y = scores[, id], Age = fx$pheno$Age)
    expect_equal(
      one[[paste0(id, "_Age_accel")]],
      unname(residuals(lm(y ~ Age, data = d)))
    )
  }

  # a formula replaces it, and moves the answer
  two <- calc_accel(fx$res, ~ Age + Female, data = fx$pheno, long = FALSE)
  d <- data.frame(
    y = scores[, "GrimAgeV1"],
    Age = fx$pheno$Age,
    Female = fx$pheno$Female
  )
  expect_equal(
    two$GrimAgeV1_Age_Female_accel,
    unname(residuals(lm(y ~ Age + Female, data = d)))
  )
  expect_false(isTRUE(all.equal(
    one$GrimAgeV1_Age_accel,
    two$GrimAgeV1_Age_Female_accel
  )))
})

test_that("diff is the raw difference, and a formula residualizes it", {
  fx <- FX
  half <- calc_clocks(
    fx$DNAm[1:6, , drop = FALSE],
    ACCEL_CLOCKS,
    pheno = fx$pheno[1:6, ]
  )

  d_full <- calc_accel(fx$res, type = "diff", long = FALSE)
  expect_equal(
    d_full$GrimAgeV1_diff,
    unname(as.matrix(fx$res)[, "GrimAgeV1"] - fx$pheno$Age)
  )
  # the per-sample property: dropping samples does not move a diff
  expect_equal(
    calc_accel(half, type = "diff", long = FALSE)$GrimAgeV1_diff,
    d_full$GrimAgeV1_diff[1:6]
  )

  # with an rhs, diff is that difference residualized
  b <- calc_accel(fx$res, ~Female, type = "diff", long = FALSE)
  d <- data.frame(
    y = as.matrix(fx$res)[, "GrimAgeV1"] - fx$pheno$Age,
    Female = fx$pheno$Female
  )
  expect_equal(
    b$GrimAgeV1_Female_diff,
    unname(residuals(lm(y ~ Female, data = d)))
  )
})

test_that("an NA covariate drops that sample, scoped to the formula's variables", {
  fx <- FX
  ph <- fx$pheno
  ph$Age[c(2L, 5L)] <- NA_real_
  expect_warning(res <- calc_clocks(fx$DNAm, ACCEL_CLOCKS, pheno = ph))

  # the pheno NA count is reported, not silent
  expect_warning(acc <- calc_accel(res, long = FALSE))
  expect_true(all(is.na(acc$GrimAgeV1_Age_accel[c(2L, 5L)])))
  expect_equal(sum(!is.na(acc$GrimAgeV1_Age_accel)), nrow(ph) - 2L)

  # and the fit is over the surviving samples only
  keep <- !is.na(ph$Age)
  d <- data.frame(y = as.matrix(res)[, "GrimAgeV1"], Age = ph$Age)[keep, ]
  expect_equal(
    acc$GrimAgeV1_Age_accel[keep],
    unname(residuals(lm(y ~ Age, data = d)))
  )

  # a covariate the formula does not name drops nothing
  ph2 <- fx$pheno
  ph2$Female[c(1L, 3L)] <- NA_real_
  expect_warning(res2 <- calc_clocks(fx$DNAm, ACCEL_CLOCKS, pheno = ph2))
  expect_equal(
    sum(is.na(calc_accel(res2, ~Age, long = FALSE)$Hannum_Age_accel)),
    0L
  )
  expect_warning(both <- calc_accel(res2, ~ Age + Female, long = FALSE))
  expect_equal(which(is.na(both$Hannum_Age_Female_accel)), c(1L, 3L))
})

test_that("data supplies what the record lacks and may never change what it has", {
  fx <- FX

  # the modal workflow: the same Age back, plus the new column
  expect_no_warning(calc_accel(fx$res, ~ Age + Female, data = fx$pheno))
  # storage is not disagreement: integer against double is the same Age
  as_int <- fx$pheno
  as_int$Age <- round(as_int$Age)
  ok <- fx$pheno
  ok$Age <- as.integer(round(ok$Age))
  int_res <- calc_clocks(fx$DNAm, ACCEL_CLOCKS, pheno = as_int)
  expect_no_warning(calc_accel(int_res, ~ Age + Female, data = ok))
  # but a type change is
  chr <- fx$pheno
  chr$Age <- as.character(chr$Age)
  expect_error(calc_accel(fx$res, ~Age, data = chr))

  # a genuinely different Age is a different pheno, not a supplement
  bad <- fx$pheno
  bad$Age[[1L]] <- bad$Age[[1L]] + 5
  expect_error(calc_accel(fx$res, ~ Age + Female, data = bad))
  # and the id column must be there to join on
  expect_error(calc_accel(fx$res, ~Age, data = fx$pheno[, "Age", drop = FALSE]))

  # a record carrying no covariate at all needs data, and refuses without it
  n <- 8L
  DNAm <- random_betas(clock_scoring_cpgs("Hannum"), n = n)
  res <- calc_clocks(DNAm, "Hannum")
  expect_error(calc_accel(res))
  # diff needs Age too, whatever the rhs says
  expect_error(calc_accel(res, type = "diff"))

  ph <- mc_pheno(rownames(DNAm), Age = stats::rnorm(n, 45, 8))
  expect_no_error(calc_accel(res, data = ph))
  expect_error(calc_accel(res, ~ Age + Smoking, data = ph))
  # a two-sided formula is not an rhs
  expect_error(calc_accel(res, score ~ Age, data = ph))
})

test_that("a sample data has no row for is reported, not silently NA-filled", {
  n <- 8L
  DNAm <- random_betas(clock_scoring_cpgs("Hannum"), n = n)
  res <- calc_clocks(DNAm, "Hannum")

  # ids that match nothing at all.
  wrong <- mc_pheno(
    sub("^sample", "Sample", rownames(DNAm)),
    Age = stats::rnorm(n, 45, 8)
  )
  msgs <- capture_warnings(
    value <- calc_accel(res, data = wrong, long = FALSE)
  )
  # pinned because three warnings fire here and only one is this behaviour
  expect_true(any(grepl("no row for", msgs)))
  expect_true(all(is.na(value$Hannum_Age_accel)))

  # a complete data says nothing
  right <- mc_pheno(rownames(DNAm), Age = stats::rnorm(n, 45, 8))
  expect_no_warning(calc_accel(res, data = right))
})

test_that("mc_batch_id is a formula variable and is reserved against data", {
  n <- 10L
  m1 <- sim_DNAm("Hannum", n = n, suffix = "T1")$DNAm
  m2 <- sim_DNAm("Hannum", n = n, suffix = "T2")$DNAm
  both <- rbind(calc_clocks(m1, "Hannum"), calc_clocks(m2, "Hannum"))
  ph <- mc_pheno(both$provenance$sample_id, Age = stats::rnorm(2L * n, 45, 8))

  # the record's own label, without the caller rebuilding it from provenance
  acc <- calc_accel(both, ~ Age + mc_batch_id, data = ph, long = FALSE)
  d <- data.frame(
    y = as.matrix(both)[, "Hannum"],
    Age = ph$Age,
    b = both$provenance$mc_batch_id
  )
  expect_equal(
    acc$Hannum_Age_mc_batch_id_accel,
    unname(residuals(lm(y ~ Age + b, data = d)))
  )

  # the name belongs to the record
  clash <- ph
  clash$mc_batch_id <- "mine"
  expect_error(calc_accel(both, ~Age, data = clash))
})

test_that("accel long and wide agree, and accel_id names the spec", {
  fx <- FX
  wide <- calc_accel(fx$res, long = FALSE)
  long <- calc_accel(fx$res)
  clocks <- colnames(as.matrix(fx$res))

  # accel_id sits beside clock_id: the two together are what a pivot keys on
  expect_equal(names(long), c("ID", "clock_id", "accel_id", "accel"))
  expect_equal(names(wide), c("ID", paste0(clocks, "_Age_accel")))
  expect_equal(nrow(long), nrow(wide) * length(clocks))
  expect_equal(long$accel, as.vector(as.matrix(wide[, -1L, drop = FALSE])))

  # one label per (type, formula), derived from the call.
  labels_of <- function(...) unique(calc_accel(fx$res, ...)$accel_id)
  expect_equal(labels_of(~ Age + Female), "Age_Female_accel")
  expect_equal(labels_of(type = "diff"), "diff")
  stacked <- rbind(calc_accel(fx$res, ~Age), calc_accel(fx$res, type = "diff"))
  expect_equal(anyDuplicated(stacked[, c("ID", "clock_id", "accel_id")]), 0L)
})

test_that("too few samples to fit gives NA and one warning, not an error", {
  skip_on_cran()
  fx <- accel_fixture(n = 2L)
  clocks <- paste0(colnames(as.matrix(fx$res)), "_Age_accel")
  expect_warning(acc <- calc_accel(fx$res, long = FALSE))

  # every column is present and every one of them is NA
  expect_equal(names(acc), c("ID", clocks))
  expect_true(all(is.na(as.matrix(acc[, clocks, drop = FALSE]))))
})

test_that("clocks with different missingness patterns each fit their own rows", {
  skip_on_cran()
  cl <- c("DNAmFitAge", "Hannum")
  n <- 12L
  DNAm <- random_betas(clock_cpgs(cl), n = n)
  ph <- mc_pheno(
    rownames(DNAm),
    Age = stats::rnorm(n, 45, 8),
    Female = rep(c(0, 1), length.out = n)
  )
  ph$Female[c(2L, 7L)] <- NA_real_
  expect_warning(res <- calc_clocks(DNAm, cl, pheno = ph))

  # ~ Age drops no rows. groups are the score NAs alone.
  scores <- as.matrix(res)
  expect_true(all(is.na(scores[c(2L, 7L), "DNAmFitAge"])))
  expect_false(anyNA(scores[, "Hannum"]))

  acc <- calc_accel(res, ~Age, long = FALSE)
  # the short group fits over the 10 rows it scored, and stays NA elsewhere
  fit <- !is.na(scores[, "DNAmFitAge"])
  d2 <- data.frame(y = scores[, "DNAmFitAge"], Age = ph$Age)[fit, ]
  expect_equal(
    acc$DNAmFitAge_Age_accel[fit],
    unname(residuals(lm(y ~ Age, data = d2)))
  )
  expect_true(all(is.na(acc$DNAmFitAge_Age_accel[!fit])))
})

test_that("a cohort-mean fill across batches is reported, never injected", {
  skip_on_cran()
  n <- 10L
  mk <- function(tag, na_frac) {
    m <- sim_DNAm("Hannum", n = n, suffix = tag)$DNAm
    if (na_frac > 0) {
      m[sample.int(length(m), floor(length(m) * na_frac))] <- NA_real_
    }
    m
  }
  ph <- function(x) {
    mc_pheno(x$provenance$sample_id, Age = stats::rnorm(nrow(x$scores), 45, 8))
  }

  # filled in both batches -> the offset is real, so say so
  dirty <- rbind(
    calc_clocks(mk("A", 0.05), "Hannum"),
    calc_clocks(mk("B", 0.05), "Hannum")
  )
  expect_message(calc_accel(dirty, ~Age, data = ph(dirty)))
  # naming it in the rhs is the fix, so the note stops
  expect_no_message(calc_accel(dirty, ~ Age + mc_batch_id, data = ph(dirty)))

  # nothing was filled -> the batches are numerically irrelevant, so stay quiet
  clean <- rbind(
    calc_clocks(mk("C", 0), "Hannum"),
    calc_clocks(mk("D", 0), "Hannum")
  )
  expect_no_message(calc_accel(clean, ~Age, data = ph(clean)))
})
