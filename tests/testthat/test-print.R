# the shared printer grammar (R/print.R). loose on purpose: the class tag and
# the return value are the contract, the counts and layout are not

test_that("every mc_* printer runs and returns its input invisibly", {
  sim <- sim_DNAm("Hannum", n = 3L)
  res <- calc_clocks(sim[["DNAm"]], "Hannum", pheno = sim[["pheno"]])
  cit <- cite_clocks("Hannum")

  for (x in list(sim, res, cit)) {
    got <- suppressMessages(capture.output(out <- print(x)))
    expect_equal(out, x)
  }

  # the cat-based printers write to stdout. mc_citation goes through cli, whose
  # output is on the message stream, so only its return value is checked above
  expect_output(print(sim), "<mc_sim>", fixed = TRUE)
  expect_output(print(res), "<mc_result>", fixed = TRUE)
})

test_that("a run with no pheno argument still carries the id column", {
  DNAm <- random_betas(clock_scoring_cpgs("Hannum"), n = 3L)
  res <- calc_clocks(DNAm, "Hannum")

  expect_equal(names(res[["pheno"]]), "ID")
  expect_equal(res[["pheno"]][["ID"]], rownames(DNAm))
  expect_output(print(res), "pheno", fixed = TRUE)
})
