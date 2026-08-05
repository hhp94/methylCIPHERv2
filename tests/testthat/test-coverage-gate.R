# upfront coverage gate: unscorable clock stops before any scoring

# a panel sharing no CpG with `ids`.
foreign_panel <- function(ids) {
  out <- setdiff(
    clock_scoring_cpgs("PedBE"),
    unlist(lapply(ids, clock_scoring_cpgs))
  )
  # a donor swallowed by the request would make every caller vacuous
  stopifnot(length(out) > 0L)
  out
}

# the leading `frac` of a panel, which is what a coverage ratio is measured on
thin_panel <- function(cpgs, frac) cpgs[seq_len(round(frac * length(cpgs)))]

test_that("under-covered clocks stop instead of scoring", {
  cpgs <- clock_scoring_cpgs("Hannum")
  keep <- thin_panel(cpgs, 0.5)
  expect_error(calc_clocks(random_betas(keep, n = 4L), "Hannum"))

  # a failing clock stops the whole call, not just its own column
  bad <- thin_panel(clock_scoring_cpgs("PedBE"), 0.3)
  expect_error(calc_clocks(
    random_betas(c(cpgs, bad), n = 4L),
    c("Hannum", "PedBE")
  ))

  # columns now clear -- every row is still half-imputed, so the row gate warns
  expect_warning(
    res <- calc_clocks(
      random_betas(keep, n = 4L),
      "Hannum",
      min_clocks_coverage = 0.4
    )
  )
  expect_true(all(is.finite(res$scores[, "Hannum"])))
})

test_that("zero observed CpGs stops even at min_clocks_coverage = 0", {
  # omit policy: fully absent panel is unscoreable at any threshold
  DNAm <- random_betas(foreign_panel("Hannum"), n = 4L)
  expect_error(calc_clocks(DNAm, "Hannum", min_clocks_coverage = 0))
})

test_that("min_clocks_coverage = 0 really is off for a vendor-filled clock", {
  # vendor_mean: threshold 0 means no gate (parity runs there)
  DNAm <- random_betas(foreign_panel("DNAmCRP"), n = 4L)

  res <- suppressWarnings(calc_clocks(
    DNAm,
    "DNAmCRP",
    min_clocks_coverage = 0,
    min_samples_coverage = 0
  ))
  expect_true(all(is.finite(res$scores[, "DNAmCRP"])))
  # every CpG came from the vendored ref, so every sample gets the same score
  expect_equal(length(unique(res$scores[, "DNAmCRP"])), 1L)

  # the gate is still a gate one notch up
  expect_error(calc_clocks(DNAm, "DNAmCRP", min_clocks_coverage = 0.01))
})

test_that("the gate names a clock the caller is allowed to request", {
  skip_on_cran()
  routed <- sex_routed_members()
  member <- names(routed$alias)[[1]]
  alias <- routed$alias[[member]]

  # thin matrix on the alias must not print any member id.
  DNAm <- random_betas(thin_panel(clock_scoring_cpgs(member), 0.3), n = 4L)
  pheno <- mc_pheno(
    rownames(DNAm),
    Age = mc_ages(4L),
    Female = c(1L, 0L, 1L, 0L)
  )
  msg <- conditionMessage(tryCatch(
    calc_clocks(DNAm, alias, pheno = pheno),
    error = function(e) e
  ))
  # it must be the coverage gate talking, not a missing-pheno abort
  expect_true(grepl("min_clocks_coverage", msg, fixed = TRUE))
  expect_true(grepl(alias, msg, fixed = TRUE))
  for (nm in names(routed$alias)) {
    expect_false(grepl(nm, msg, fixed = TRUE))
  }
})

test_that("a sparse normalization panel warns but still scores (does not stop)", {
  skip_if_not_installed("betanorm")
  gold <- names(clock_norm_target("DunedinPACE"))
  model <- clock_scoring_cpgs("DunedinPACE")

  keep <- union(model, thin_panel(gold, 0.5))
  # two distinct warnings: thin background (column) and per-sample imputation (row)
  expect_warning(
    expect_warning(
      res <- calc_clocks(random_betas(keep, n = 4L), "DunedinPACE")
    )
  )
  expect_true(all(is.finite(res$scores[, "DunedinPACE"])))
})

test_that("the warn band sits above min_clocks_coverage and is silent at full coverage", {
  cpgs <- clock_scoring_cpgs("Hannum")
  keep <- thin_panel(cpgs, 0.78)
  DNAm <- random_betas(keep, n = 4L)

  # clearing the default floor by under 10% warns instead of stopping
  expect_warning(res <- calc_clocks(DNAm, "Hannum"))
  expect_true(all(is.finite(res$scores[, "Hannum"])))

  # the band moves with the floor, and never fires on a full panel
  expect_silent(calc_clocks(DNAm, "Hannum", min_clocks_coverage = 0.5))
  expect_silent(calc_clocks(random_betas(cpgs, n = 4L), "Hannum"))
})

# row gate reads every branch's coverage record, not just Dunedin
test_that("an under-covered sample warns on an ordinary linear clock", {
  cpgs <- clock_scoring_cpgs("Hannum")
  DNAm <- random_betas(cpgs, n = 4L)
  DNAm[1, seq_len(round(0.5 * ncol(DNAm)))] <- NA

  expect_warning(res <- calc_clocks(DNAm, "Hannum"))
  expect_true(all(is.finite(res$scores[, "Hannum"])))

  expect_silent(calc_clocks(DNAm, "Hannum", min_samples_coverage = 0))
})
