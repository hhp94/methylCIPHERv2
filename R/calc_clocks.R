# public front door

#' Epigenetic Clock Scores
#'
#' Scores CpG-based epigenetic clocks on a matrix of methylation beta values.
#'
#' @inheritParams mc-params
#' @param pheno_id A string. The name of the column in `pheno` that holds the
#'   sample ids. Default is `"ID"`.
#' @param min_clocks_coverage A number between 0 and 1. The smallest fraction
#'   of a clock's CpGs that must be present for that clock to score. Default is
#'   `0.75`.
#' @param min_samples_coverage A number between 0 and 1. The smallest fraction
#'   of a clock's CpGs that must be present for a sample to score without a
#'   warning. Default is `0.75`.
#'
#' @inheritSection mc-params The assets directory
#'
#' @details
#' [list_clocks()] and [list_clock_tags()] show every value `clocks` accepts.
#'
#' `normalize` turns on the schemes that a clock declares as optional. It
#' cannot turn off a scheme that is part of the clock. The `normalize` column
#' of `list_clocks(all_columns = TRUE)` gives the scheme each clock uses.
#'
#' The two coverage arguments differ. `calc_clocks()` stops when a clock has
#' too few CpGs present, so every clock in the returned scores passed
#' `min_clocks_coverage`. A clock just above that floor, and a clock whose
#' normalization panel falls below it, each raise a warning and still score.
#' A sample with too few CpGs present raises a warning and still scores. Pass
#' the returned value to [clocks_coverage()] or [samples_coverage()] to see
#' the counts.
#'
#' `calc_clocks()` narrows `pheno` before it stores it. The returned value
#' keeps the id column and the covariates that the clocks need, and drops the
#' other columns.
#'
#' @returns An `mc_result` object. It holds the scores, the narrowed `pheno`,
#'   the coverage counts, and the provenance of the run.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 20)
#'
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#' res
#'
#' # pheno is narrowed to the id column and the covariates the clocks need
#' pheno <- data.frame(ID = rownames(sim[["DNAm"]]), Age = runif(20, 20, 80))
#' res <- calc_clocks(sim[["DNAm"]], clocks, pheno = pheno)
#' head(res[["pheno"]])
#'
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
  # the boundary: everything downstream assumes these are already validated
  checkmate::assert_string(pheno_id)
  checkmate::assert_number(min_clocks_coverage, lower = 0, upper = 1)
  checkmate::assert_number(min_samples_coverage, lower = 0, upper = 1)

  spec <- mc_spec(clocks, pheno_id, normalize, ext_data, ask)
  facts <- mc_cohort(DNAm, spec, pheno, min_clocks_coverage)
  scored <- score_cohort(DNAm, spec, facts)
  # shared with refinalize_clocks() -- a no-op when pending is empty
  scores <- finalize_cross_sample(scored[["scores"]], scored[["pending"]])

  # per-sample coverage gate (warn only, after scoring)
  check_row_coverage(scored[["coverage"]], min_samples_coverage)
  # value gate on output columns. nan/inf land here.
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
    min_clocks_coverage = min_clocks_coverage,
    min_samples_coverage = min_samples_coverage,
    scoring_failures = scored[["notes"]],
    # kept, not discarded, so a bound record can re-finalize exactly
    pending = scored[["pending"]]
  )
}
