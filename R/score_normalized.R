# normalize-then-linear over the clock's declared panel and target

score_normalized <- function(id, cpgs, block, results) {
  # declined or not requested -- score raw
  if (!cpgs[["normalizes"]]) {
    return(linear_score(cpgs, block))
  }
  require_betanorm(id)

  scheme <- clock_norm_scheme(id)
  target <- clock_norm_target(id)
  # scoring CpGs are a subset of the background panel
  obs <- observed_panel(cpgs[["norm_present"]], block)

  # unimplemented scheme already stopped in score_type()
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

# BMIQ onto vendored gold (absent probes dropped, unfit samples -> NA + notes)
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
