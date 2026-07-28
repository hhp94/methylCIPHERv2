# unroutable catalog entry. Reads the four routing facts off the entry itself so
# no call site can report a subset of them.
unroutable <- function(p) {
  cli::cli_abort(
    c(
      "No scoring path for clock {.val {p}}.",
      "*" = "group {.val {clock_group_id(p)}},
             weights_format {.val {clock_weights_format(p)}},
             computation_type {.val {clock_type(p)}},
             normalization {.val {clock_norm_scheme(p)}}",
      "i" = "This is a package bug -- please report it."
    ),
    call = NULL
  )
}

# scorer tag for calc_clocks() dispatch
score_type <- function(p) {
  # package-minted aliases route on `kind` before weights_format / computation_type
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

  # a declared normalization scheme routes to normalize-then-linear; the
  # Dunedin group keeps its own tag above. Catalog-only, never the request --
  # the branch itself reads `cpgs$normalizes` to decide whether to apply it.
  #
  # A scheme score_normalized() does not implement stops here rather than
  # falling through: the fall-through lands on `linear`, which would score the
  # clock raw and silently drop a normalization step the catalog declared.
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
#' @param from x
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
  from = NULL,
  ask = TRUE
) {
  checkmate::assert_number(min_samples_coverage, lower = 0, upper = 1)

  spec <- mc_spec(clocks, pheno_id, normalize, from, ask)
  facts <- mc_cohort(DNAm, spec, pheno, min_clocks_coverage)
  scored <- score_cohort(DNAm, spec, facts)
  # one block, so assembly is the identity -- but the reduction still runs
  # here rather than in the loop, so both front ends share this step
  scores <- finalize_cross_sample(scored[["scores"]], scored[["pending"]])

  # row-gate every clock actually computed (including routed members). Warn
  # only, so it reads assembled counts and runs after scoring.
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

  # sample_miss is per panel role: score for every column, norm for normalizers only
  per_clock <- coverage[["per_clock"]]
  # normalizers from the record's `normalizes` flag (aliases are NULL / dropped)
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
        # which clocks were actually normalized -- absent means scored raw
        normalized = normalized,
        # clock id -> sample ids the scorer could not fit, so an NA score is
        # readable as a failure rather than as absent input. Empty is the
        # normal case; coverage cannot carry this (it precedes scoring).
        scoring_failures = scoring_failures
      )
    ),
    class = "mc_result"
  )
}
