# clocks_coverage() / samples_coverage() formatters over a finished record.

report_pheno <- function(ids, female) {
  data.frame(
    ID = ids,
    Female = as.integer(female),
    Age = seq(45, 70, length.out = length(ids)),
    stringsAsFactors = FALSE
  )
}

test_that("clocks_coverage marks members routing_target and the alias returned", {
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 4L)
  pheno <- report_pheno(rownames(DNAm), c(1, 1, 0, 0))

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)
  cc <- clocks_coverage(res)

  expect_s3_class(cc, "data.frame")
  expect_true(all(
    c("clock_id", "role", "policy", "score_needed", "score_used") %in% names(cc)
  ))

  # the alias is a returned column with an all-NA panel row
  alias <- cc[cc$clock_id == "DNAmGrip_wAge", ]
  expect_equal(alias$role, "returned")
  expect_true(is.na(alias$score_needed))

  # the two members are routing targets carrying the per-sex denominators
  members <- cc[
    cc$clock_id %in% c("DNAmGrip_wAge_Female", "DNAmGrip_wAge_Male"),
  ]
  expect_equal(sort(unique(members$role)), "routing_target")
  expect_setequal(members$score_needed, c(length(fem), length(mal)))
})

test_that("clocks_coverage reports score_used = present + imputed_full per row", {
  id <- "DNAmGait_noAge_Female"
  coef <- clock_coefs(id)
  medians <- clock_impute_ref(id)
  full <- union(names(coef), names(clock_coefs("DNAmGait_noAge_Male")))
  DNAm <- random_betas(full, n = 4L)
  drop <- intersect(names(coef), names(medians))[1:5]
  DNAm <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]
  pheno <- report_pheno(rownames(DNAm), rep(1L, 4L))

  cc <- clocks_coverage(calc_clocks(DNAm, "DNAmGait_noAge", pheno = pheno))
  row <- cc[cc$clock_id == id, ]
  expect_equal(row$score_imputed_full, 5L)
  expect_equal(row$score_used, row$score_present + row$score_imputed_full)
  expect_equal(row$missing_cpgs[[1]], drop)
})

test_that("samples_coverage is long with a per-sex denominator for an alias", {
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- report_pheno(rownames(DNAm), female)

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)
  sc <- samples_coverage(res)

  expect_equal(
    names(sc),
    c("id", "clock_id", "panel", "n_observed", "n_needed", "coverage")
  )
  alias <- sc[sc$clock_id == "DNAmGrip_wAge", ]
  expect_equal(nrow(alias), 6L)
  # n_needed is the scoring model of the sex that scored each sample
  expect_equal(
    alias$n_needed,
    ifelse(female == 1, length(fem), length(mal))
  )
  # full coverage: every declared CpG present
  expect_true(all(alias$coverage == 1))
  expect_equal(alias$n_observed, alias$n_needed)
})

test_that("samples_coverage coverage is literally row_coverage() for a partial fill", {
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- report_pheno(rownames(DNAm), female)
  fem_only <- setdiff(names(fem), names(mal))[1]
  DNAm[1, fem_only] <- NA_real_ # one female leans on a cohort mean

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)
  sc <- samples_coverage(res)
  alias <- sc[sc$clock_id == "DNAmGrip_wAge", ]

  # sample 1 is one CpG short of its female panel, the rest are full
  expect_equal(alias$n_observed[1], length(fem) - 1L)
  expect_equal(alias$coverage[1], (length(fem) - 1L) / length(fem))
  expect_true(all(alias$coverage[-1] == 1))
})

test_that("samples_coverage gives a normalizing clock a score and a norm row", {
  skip_if_not_installed("betanorm")
  norm_panel <- names(clock_norm_target("DunedinPACE"))
  score_panel <- clock_scoring_cpgs("DunedinPACE")
  DNAm <- random_betas(norm_panel, n = 4L)

  res <- calc_clocks(DNAm, "DunedinPACE", min_samples_coverage = 0)
  sc <- samples_coverage(res)
  dp <- sc[sc$clock_id == "DunedinPACE", ]

  expect_setequal(unique(dp$panel), c("score", "norm"))
  needed <- tapply(dp$n_needed, dp$panel, unique)
  expect_equal(as.integer(needed[["score"]]), length(score_panel))
  # the norm background is the larger, separate panel
  expect_gt(needed[["norm"]], needed[["score"]])
})
