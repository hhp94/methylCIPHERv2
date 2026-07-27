# crash-smoke: every bundled callable clock scores on random betas without error

bundled_smoke_clocks <- function() {
  # callable pool only (routed members exercise via their alias)
  ids <- resolve_clocks("all")
  ids[!vapply(ids, clock_is_external, logical(1))]
}

# betanorm soft dep for QN clocks.
betanorm_installed <- requireNamespace("betanorm", quietly = TRUE)

for (id in bundled_smoke_clocks()) {
  local({
    clock_id <- id

    # does this clock normalize in the default configuration?
    needs_betanorm <- isTRUE(resolve_normalize(NULL, clock_id)[[clock_id]])
    test_that(paste0("sim_DNAm smoke: ", clock_id), {
      if (needs_betanorm) {
        skip_if_not(betanorm_installed, "betanorm not installed")
      }
      sim <- sim_DNAm(clock_id, n = 4L, Age = TRUE, Female = TRUE)
      expect_no_error(
        calc_clocks(sim$DNAm, clock_id, pheno = sim$pheno)
      )
    })
  })
}
