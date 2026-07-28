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

test_that("augment refuses duplicate ids and warns on column collisions", {
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 5L, Age = TRUE, Female = TRUE)
  res <- suppressWarnings(calc_clocks(sim$DNAm, c("Horvath1", "Hannum"), pheno = sim$pheno))

  dup <- data.frame(ID = rep(sim$pheno$ID[1L], 2L), X = 1:2)
  expect_error(augment(res, data = dup)) # would fan out the join

  clash <- data.frame(ID = sim$pheno$ID, Horvath1 = 1)
  expect_warning(augment(res, data = clash)) # column already on the table
})

test_that("augment(adjust=) appends residual columns orthogonal to the covariates", {
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 12L, Age = TRUE, Female = TRUE)
  res <- calc_clocks(sim$DNAm, c("Horvath1", "Hannum"), pheno = sim$pheno)

  # these clocks need no covariates, so Age is not on the record -> needs data=
  expect_error(augment(res, adjust = "Age"))

  a <- augment(res, data = sim$pheno, adjust = c("Age", "Female"))
  expect_true(all(c("Horvath1_resid", "Hannum_resid") %in% names(a)))
  # a residual is orthogonal to its adjustment set
  expect_lt(abs(stats::cor(a$Horvath1_resid, a$Age)), 1e-6)
})

test_that("augment accepts a formula for adjust and rejects other types clearly", {
  sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 10L, Age = TRUE, Female = TRUE)
  res <- calc_clocks(sim$DNAm, c("Horvath1", "Hannum"), pheno = sim$pheno)
  a <- augment(res, data = sim$pheno, adjust = ~ Age + Female)
  expect_true(all(c("Horvath1_resid", "Hannum_resid") %in% names(a)))
  expect_error(augment(res, data = sim$pheno, adjust = 1)) # not names / formula
})

test_that("codebook / bibliography reject non-clock input clearly", {
  expect_error(codebook(NULL))
  expect_error(codebook(123))
  expect_error(bibliography(NULL))
})

test_that("augment(adjust=) uses record covariates when a clock required them", {
  sim <- sim_DNAm("DNAmGrip_wAge", n = 8L, Age = TRUE, Female = TRUE)
  res <- suppressWarnings(calc_clocks(sim$DNAm, "DNAmGrip_wAge", pheno = sim$pheno))
  a <- augment(res, adjust = c("Age", "Female")) # no data= needed
  expect_true("DNAmGrip_wAge_resid" %in% names(a))
})

test_that("codebook reports panel size for sex-routed aliases (not 0)", {
  cb <- codebook("DNAmFitAge")
  expect_true(all(cb$n_cpgs > 0L))
})

test_that("bibliography enriches from clocks.bib and emits full BibTeX", {
  b <- bibliography(c("Horvath1", "Hannum"))
  expect_true(all(c("reference", "citation", "title", "journal", "year", "doi", "pmid", "url")
  %in% names(b)))
  expect_true(all(grepl("pubmed", b$url)))
  # both are in clocks.bib, so titles/journals resolve
  expect_false(any(is.na(b$title)))
  expect_false(any(is.na(b$journal)))

  bt <- utils::capture.output(bibliography("Horvath1", format = "bibtex"))
  expect_true(any(grepl("@article", bt)))
  expect_true(any(grepl("title = \\{", bt))) # full entry, not just a stub
})
