# normalize-then-linear: the clock's declared scheme applied to its own declared
# panel and target, then the shared linear engine over the calibrated betas.

score_normalized <- function(id, cpgs, block, results) {
  # declined, or a clock that declares a scheme but was not asked to apply it
  if (!cpgs[["normalizes"]]) {
    return(linear_score(cpgs, block))
  }
  require_betanorm(id)

  scheme <- clock_norm_scheme(id)
  target <- clock_norm_target(id)
  # scoring CpGs are a subset of the background panel, so one pass over the
  # present background covers both
  obs <- observed_panel(cpgs[["norm_present"]], block)

  # routing totality is score_type()'s job -- a scheme with no arm here stops
  # there, before any DNAm is read, not part-way through the scoring loop
  calibrated <- switch(
    scheme,
    bmiq = bmiq_panel(obs, target, id, block),
    cli::cli_abort(
      "No normalization branch for scheme {.val {scheme}} (clock
       {.val {id}}).",
      call = NULL
    )
  )

  linear_score(
    cpgs,
    block,
    observed = list(
      cols = cpgs[["score_present"]],
      values = calibrated[, cpgs[["score_present"]], drop = FALSE]
    )
  )
}

# BMIQ calibration onto the vendored gold standard, at Horvath's fixed settings.
# A fully absent probe is dropped rather than filled from the target: BMIQ fits
# each sample's own mixture from the panel, so target-drawn values pull that fit
# toward the gold standard and shrink the correction being computed.
#
# A sample whose mixture will not fit is calibrated to NA, so its score is NA
# while its coverage -- computed before any scoring -- still reads full. Say so
# and record it, or the record cannot tell that NA from a missing covariate.
bmiq_panel <- function(obs, target, id, block) {
  fit <- betanorm::bmiq_calibration(
    obs[["values"]],
    goldstandard.beta = target[obs[["cols"]]],
    nfit = ncol(obs[["values"]]),
    verbose = FALSE,
    on.sample.error = "continue",
    failed.sample = "NA"
  )

  failed <- block[["sample_id"]][!fit[["success"]]]
  if (length(failed)) {
    note_scoring_failure(block, id, failed)
    cli::cli_warn(
      c(
        "BMIQ calibration failed for {length(failed)} sample{?s} scoring
         {.val {id}}.",
        "!" = "Scored NA: {.val {failed}}.",
        "i" = "Also recorded in {.code $provenance$scoring_failures}."
      ),
      call = NULL
    )
  }
  fit[["calibrated"]]
}
