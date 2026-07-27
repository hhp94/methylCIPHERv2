# normalize-then-linear: the clock's declared scheme applied to its own declared
# panel and target, then the shared linear engine over the calibrated betas.

score_normalized <- function(
  id,
  cpgs,
  DNAm,
  partial_cache = NULL,
  pheno = NULL,
  packs = NULL
) {
  # declined, or a clock that declares a scheme but was not asked to apply it
  if (!isTRUE(cpgs$normalizes)) {
    return(linear_score(cpgs, DNAm, partial_cache, pheno, packs))
  }
  require_betanorm(id)

  scheme <- clock_norm_scheme(id)
  target <- clock_norm_target(id)
  # scoring CpGs are a subset of the background panel, so one pass over the
  # present background covers both
  obs <- observed_panel(cpgs$norm_present, DNAm, partial_cache)

  calibrated <- switch(
    scheme,
    bmiq = bmiq_panel(obs, target),
    cli::cli_abort(
      "No normalization branch for scheme {.val {scheme}} (clock
       {.val {id}}).",
      call = NULL
    )
  )

  linear_score(
    cpgs,
    DNAm,
    partial_cache,
    pheno,
    packs,
    observed = list(
      cols = cpgs$score_present,
      values = calibrated[, cpgs$score_present, drop = FALSE]
    )
  )
}

# BMIQ calibration onto the vendored gold standard, at Horvath's fixed settings.
# A fully absent probe is dropped rather than filled from the target: BMIQ fits
# each sample's own mixture from the panel, so target-drawn values pull that fit
# toward the gold standard and shrink the correction being computed.
bmiq_panel <- function(obs, target) {
  fit <- betanorm::bmiq_calibration(
    obs$values,
    goldstandard.beta = target[obs$cols],
    nfit = ncol(obs$values),
    verbose = FALSE,
    on.sample.error = "continue",
    failed.sample = "NA"
  )
  fit$calibrated
}
