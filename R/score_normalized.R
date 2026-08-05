# normalize-then-linear over the clock's declared panel and target

# betanorm is a soft dep (every normalizing branch needs it)
require_betanorm <- function(id) {
  if (!requireNamespace("betanorm", quietly = TRUE)) {
    stop(
      sprintf(
        paste0(
          "%s needs the betanorm package for normalization. ",
          "Install it from GitHub: pak::pak(\"hhp94/betanorm\")."
        ),
        id
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

score_normalized <- function(id, cpgs, block, results) {
  # declined or not requested -- score raw
  if (!cpgs[["normalizes"]]) {
    return(linear_score(cpgs, block))
  }
  require_betanorm(id)

  scheme <- clock_norm_scheme(id)
  target <- clock_norm_target(id)
  # scoring CpGs are a subset of the background panel
  obs <- observed_panel(
    cpgs[["norm_present"]],
    cpgs[["norm_present_idx"]],
    block
  )

  # unimplemented scheme already stopped in score_type()
  calibrated <- switch(
    scheme,
    bmiq = bmiq_panel(obs, target, id, block),
    stop(
      sprintf(
        "No normalization branch for scheme %s (clock %s).",
        scheme,
        id
      ),
      call. = FALSE
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

# bmiq onto vendored gold (absent probes dropped, unfit samples -> NA + notes)
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
    warning(
      sprintf(
        paste0(
          "BMIQ calibration failed for %d sample(s) scoring %s. ",
          "Scored NA: %s. Also recorded in $provenance$scoring_failures."
        ),
        length(failed),
        id,
        paste(failed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  fit[["calibrated"]]
}
