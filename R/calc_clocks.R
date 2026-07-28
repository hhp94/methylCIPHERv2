# stop for an unroutable catalog entry (names all four routing facts)
unroutable <- function(p) {
  stop(
    sprintf(
      paste0(
        "No scoring path for clock %s ",
        "(group %s, weights_format %s, computation_type %s, ",
        "normalization %s). This is a package bug -- please report it."
      ),
      p,
      clock_group_id(p),
      clock_weights_format(p),
      clock_type(p),
      clock_norm_scheme(p)
    ),
    call. = FALSE
  )
}

# scorer tag for calc_clocks() dispatch
score_type <- function(p) {
  # package-minted aliases route on kind first
  if (identical(clock_kind(p), "sex_routed_alias")) {
    return("sex_routed")
  }

  ct <- clock_type(p)
  wf <- clock_weights_format(p)
  gid <- clock_group_id(p)

  if (clock_is_external(p)) {
    if (identical(gid, "SystemsAge")) {
      return("pack_systemsage")
    }
    if (wf == "cpg_coefficient" && ct %in% c("linear", "linear_transformed")) {
      return("pack_linear")
    }
    unroutable(p)
  }

  # group-specific tags first
  gtag <- switch(
    gid,
    Dunedin = "Dunedin",
    Zhang2019 = "Zhang2019",
    GrimAge = switch(
      ct,
      linear = "linear",
      linear_transformed = "GrimAge",
      NULL
    ),
    DNAmFitAge = switch(
      ct,
      linear = "linear",
      linear_transformed = "DNAmFitAge",
      NULL
    ),
    PhysAge = switch(ct, linear_transformed = "PhysAge", NULL),
    EpiTOC2 = switch(ct, reference_code_required = "EpiTOC2", NULL),
    MiAge = switch(ct, reference_code_required = "MiAge", NULL),
    CellDRIFT = switch(ct, reference_code_required = "linear", NULL),
    NULL
  )
  if (!is.null(gtag)) {
    return(gtag)
  }

  # declared scheme (except Dunedin) -> normalize-then-linear (else stop)
  scheme <- clock_norm_scheme(p)
  if (scheme %in% NORM_SCHEMES) {
    if (!scheme %in% NORM_SCHEMES_ROUTED) {
      unroutable(p)
    }
    return("normalized")
  }

  if (wf == "cpg_coefficient" && ct %in% c("linear", "linear_transformed")) {
    return("linear")
  }
  unroutable(p)
}

# pack groups use score_pack_group()
PACK_SCORE_TYPES <- c("pack_linear", "pack_systemsage")

is_pack_scored <- function(p) {
  score_type(p) %in% PACK_SCORE_TYPES
}

# external pack groups needed for a compute sequence
pack_groups_needed <- function(clock_sequence) {
  unique(unlist(lapply(clock_sequence, function(p) {
    if (is_pack_scored(p)) clock_group_id(p) else NULL
  })))
}


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
  # single-pass: finalize still runs here so both front ends share it
  scores <- finalize_cross_sample(scored[["scores"]], scored[["pending"]])

  # per-sample coverage gate (warn only, after scoring)
  check_row_coverage(scored[["coverage"]], min_samples_coverage)

  construct_mc_result(
    scores,
    scored[["coverage"]],
    spec[["output_ids"]],
    spec[["clock_ids"]],
    facts[["sample_id"]],
    pheno = facts[["pheno"]],
    pheno_id = pheno_id,
    covariates_used = spec[["covariates"]],
    normalized = names(spec[["normalize"]])[spec[["normalize"]]],
    scoring_failures = scored[["notes"]]
  )
}

# n x length(ids) per-sample miss matrix (NULL entry -> NA)
miss_matrix <- function(miss_list, ids, sample_id) {
  m <- matrix(
    NA_integer_,
    nrow = length(sample_id),
    ncol = length(ids),
    dimnames = list(sample_id, ids)
  )
  for (id in ids) {
    v <- miss_list[[id]]
    if (!is.null(v)) {
      m[, id] <- v
    }
  }
  m
}

# stack scorer outputs into an mc_result record
construct_mc_result <- function(
  results,
  coverage,
  output_ids,
  requested_ids,
  sample_id,
  pheno,
  pheno_id,
  covariates_used,
  normalized,
  scoring_failures = list()
) {
  scores <- do.call(cbind, results[output_ids])
  dimnames(scores) <- list(sample_id, output_ids)

  # sample_miss: score for every column, norm for normalizers only
  per_clock <- coverage[["per_clock"]]
  # normalizers from the record's normalizes flag
  norm_ids <- output_ids[vapply(
    output_ids,
    function(id) isTRUE(per_clock[[id]][["normalizes"]]),
    logical(1L)
  )]
  sample_miss <- list(
    score = miss_matrix(coverage[["sample_miss"]][["score"]], output_ids, sample_id),
    norm = miss_matrix(coverage[["sample_miss"]][["norm"]], norm_ids, sample_id)
  )

  structure(
    list(
      scores = scores,
      pheno = pheno,
      coverage = list(per_clock = per_clock, sample_miss = sample_miss),
      provenance = list(
        sample_id = sample_id,
        pheno_id = pheno_id,
        clocks = output_ids,
        requested = requested_ids,
        dependencies = setdiff(output_ids, requested_ids),
        covariates_used = covariates_used,
        # which clocks were actually normalized
        normalized = normalized,
        # clock id -> sample ids the scorer could not fit
        scoring_failures = scoring_failures
      )
    ),
    class = "mc_result"
  )
}
