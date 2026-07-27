# result-object verbs: as.matrix / as.data.frame / augment, and the catalog
# verbs codebook / bibliography. Assert what they produce.

test_that("as.matrix / as.data.frame extract scores; as.data.frame carries no pheno", {
  sim <- sim_DNAm("DNAmGrip_wAge", n = 4L, Age = TRUE, Female = TRUE)
  res <- suppressWarnings(calc_clocks(sim$DNAm, "DNAmGrip_wAge", pheno = sim$pheno))

  m <- as.matrix(res)
  expect_true(is.matrix(m))
  expect_identical(dim(m), dim(res$scores))

  df <- as.data.frame(res)
  expect_s3_class(df, "data.frame")
  # id column + one column per clock, and nothing else -- no pheno leak
  expect_identical(names(df), c("ID", "DNAmGrip_wAge"))
  expect_false(any(c("Age", "Female") %in% names(df)))
  expect_setequal(df$ID, rownames(res$scores))
})

test_that("augment surfaces the run's covariates and joins external data by id", {
  sim <- sim_DNAm("DNAmGrip_wAge", n = 4L, Age = TRUE, Female = TRUE)
  res <- suppressWarnings(calc_clocks(sim$DNAm, "DNAmGrip_wAge", pheno = sim$pheno))

  a <- augment(res)
  expect_true(all(c("Age", "Female") %in% names(a)))

  extra <- data.frame(ID = sim$pheno$ID, BMI = c(20, 25, 30, 22))
  a2 <- augment(res, data = extra)
  expect_true("BMI" %in% names(a2))
  expect_equal(a2$BMI[match(extra$ID, a2$ID)], extra$BMI)

  expect_error(augment(res, data = data.frame(x = 1))) # missing id column
})

test_that("codebook describes the requested clocks", {
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 3L)
  res <- calc_clocks(sim$DNAm, c("Horvath1", "Hannum"))

  cb <- codebook(res)
  expect_setequal(cb$clock_id, c("Horvath1", "Hannum"))
  expect_true(all(c("group_id", "n_cpgs", "computation", "pmid", "reference") %in% names(cb)))
  expect_identical(
    cb$n_cpgs[cb$clock_id == "Hannum"],
    length(clock_scoring_cpgs("Hannum"))
  )
  expect_gte(nrow(codebook("all")), 60L)
})

test_that("bibliography returns unique references with PubMed links and BibTeX", {
  b <- bibliography(c("Horvath1", "Hannum"))
  expect_true(all(c("reference", "citation", "pmid", "url") %in% names(b)))
  expect_true(all(grepl("pubmed", b$url)))
  expect_gte(nrow(b), 1L)

  bt <- utils::capture.output(bibliography("Horvath1", format = "bibtex"))
  expect_true(any(grepl("@article", bt)))
})
