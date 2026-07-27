# MiAge optimizer branch (parity owns the numeric golden).

p <- miage_params("MiAge")

# betas the model generates exactly at the given divisions (zero residual).
miage_model_betas <- function(divisions, cpgs) {
  m <- t(vapply(
    divisions,
    function(n) p$c[cpgs] + p$b[cpgs]^(n - 1) * p$d[cpgs],
    numeric(length(cpgs))
  ))
  dimnames(m) <- list(paste0("sample", seq_along(divisions)), cpgs)
  m
}

test_that("MiAge recovers the divisions that generated the betas", {
  panel <- clock_scoring_cpgs("MiAge")
  divisions <- c(300, 620, 900, 1400)
  DNAm <- miage_model_betas(divisions, panel)

  got <- calc_clocks(DNAm, "MiAge")$scores[, "MiAge"]
  expect_equal(unname(got), divisions, tolerance = 1e-6)
})

test_that("absent MiAge CpGs drop out of the objective", {
  panel <- clock_scoring_cpgs("MiAge")
  divisions <- c(450, 1100)
  DNAm <- miage_model_betas(divisions, panel)
  kept <- setdiff(panel, panel[1:40])

  res <- calc_clocks(DNAm[, kept, drop = FALSE], "MiAge")
  expect_equal(unname(res$scores[, "MiAge"]), divisions, tolerance = 1e-6)

  cov <- res$coverage$per_clock[["MiAge"]]
  expect_equal(cov$score_dropped, 40L)
  expect_equal(cov$score_present, length(kept))
  expect_equal(cov$score_imputed_full, 0L)
})
