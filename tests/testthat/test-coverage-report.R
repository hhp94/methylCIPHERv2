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

test_that("samples_coverage emits no rows for a clock that reads no CpGs", {
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 6L)
  female <- c(1, 1, 1, 0, 0, 0)
  pheno <- report_pheno(rownames(DNAm), female)

  res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno)
  sc <- samples_coverage(res)

  expect_equal(
    names(sc),
    c("id", "clock_id", "panel", "n_observed", "n_needed", "coverage", "batch")
  )
  # the alias is a returned column but reads no CpGs, so it reports nothing
  expect_false("DNAmGrip_wAge" %in% sc$clock_id)

  # rows are the members' (the clocks that read betas). masked sex rows drop,
  # so each sample appears once under the model that scored it
  expect_setequal(
    unique(sc$clock_id),
    c("DNAmGrip_wAge_Female", "DNAmGrip_wAge_Male")
  )
  expect_equal(nrow(sc), 6L)
  expect_setequal(sc$id, rownames(DNAm))
  expect_false(anyNA(sc$coverage))
  expect_equal(
    sc$clock_id[match(rownames(DNAm), sc$id)],
    ifelse(female == 1, "DNAmGrip_wAge_Female", "DNAmGrip_wAge_Male")
  )
  expect_equal(
    sc$n_needed[match(rownames(DNAm), sc$id)],
    ifelse(female == 1, length(fem), length(mal))
  )
})

test_that("a sample no model scored gets no samples_coverage row", {
  fem <- clock_coefs("DNAmGrip_wAge_Female")
  mal <- clock_coefs("DNAmGrip_wAge_Male")
  DNAm <- random_betas(union(names(fem), names(mal)), n = 4L)
  # the last sample's sex is unknown, so neither member claims it
  pheno <- report_pheno(rownames(DNAm), c(1, 1, 0, NA))

  # the unknown sex is worth a warning, and it is not what this test asserts
  res <- NULL
  expect_warning(res <- calc_clocks(DNAm, "DNAmGrip_wAge", pheno = pheno))
  sc <- samples_coverage(res)

  unscored <- rownames(DNAm)[4]
  expect_true(is.na(res$scores[unscored, "DNAmGrip_wAge"]))
  expect_false(unscored %in% sc$id)
  expect_setequal(sc$id, rownames(DNAm)[1:3])
})

test_that("samples_coverage coverage is literally row_coverage() for a partial fill", {
  panel <- clock_scoring_cpgs("Hannum")
  DNAm <- random_betas(panel, n = 6L)
  DNAm[1, panel[[1]]] <- NA_real_ # one sample leans on a cohort mean

  sc <- samples_coverage(calc_clocks(DNAm, "Hannum", min_samples_coverage = 0))
  hannum <- sc[sc$clock_id == "Hannum", ]

  expect_equal(nrow(hannum), 6L)
  # sample 1 is one CpG short of the panel, the rest are full
  expect_equal(hannum$n_observed[1], length(panel) - 1L)
  expect_equal(hannum$coverage[1], (length(panel) - 1L) / length(panel))
  expect_true(all(hannum$coverage[-1] == 1))
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
