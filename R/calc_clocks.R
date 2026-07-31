# public front door

#' Stub
#'
#' Rcpp needs some roxygen2 stub
#'
#' @param DNAm x
#' @param clocks x
#' @param pheno x
#' @param pheno_id x
#' @param min_clocks_coverage x
#' @param min_samples_coverage x
#' @param normalize x
#' @param ext_data x
#' @param ask x
#'
#' @returns x
#' @export
calc_clocks <- function(
  DNAm,
  clocks,
  pheno = NULL,
  pheno_id = "ID",
  min_clocks_coverage = 0.75,
  min_samples_coverage = 0.75,
  normalize = NULL,
  ext_data = NULL,
  ask = TRUE
) {
  checkmate::assert_number(min_samples_coverage, lower = 0, upper = 1)

  spec <- mc_spec(clocks, pheno_id, normalize, ext_data, ask)
  facts <- mc_cohort(DNAm, spec, pheno, min_clocks_coverage)
  scored <- score_cohort(DNAm, spec, facts)
  # shared with refinalize_clocks() -- a no-op when pending is empty
  scores <- finalize_cross_sample(scored[["scores"]], scored[["pending"]])

  # per-sample coverage gate (warn only, after scoring)
  check_row_coverage(scored[["coverage"]], min_samples_coverage)
  # value gate on output columns. nan/inf land here, not in the panel scan
  check_score_values(scores[spec[["output_ids"]]])

  construct_mc_result(
    scores,
    scored[["coverage"]],
    spec[["output_ids"]],
    spec[["clock_ids"]],
    facts[["sample_id"]],
    pheno = facts[["pheno"]],
    pheno_id = spec[["pheno_id"]],
    covariates_used = spec[["covariates"]],
    normalized = names(spec[["normalize"]])[spec[["normalize"]]],
    scoring_failures = scored[["notes"]],
    # kept, not discarded, so a bound record can re-finalize exactly
    pending = scored[["pending"]]
  )
}
