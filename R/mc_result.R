# mc_result record: constructor + methods for class "mc_result"

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

  # per_clock, sample_miss and samples_coverage() span one set: clocks that
  # read CpGs. pure composites count nothing of their own, so they are in none
  per_clock <- coverage[["per_clock"]]
  record_ids <- names(per_clock)[
    !vapply(per_clock, is.null, logical(1L))
  ]
  # normalizers from the record's normalizes flag
  norm_ids <- record_ids[vapply(
    record_ids,
    function(id) isTRUE(per_clock[[id]][["normalizes"]]),
    logical(1L)
  )]
  sample_miss <- list(
    score = miss_matrix(
      coverage[["sample_miss"]][["score"]],
      record_ids,
      sample_id
    ),
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

# print mc_result with labelled components
#' @export
print.mc_result <- function(x, n = 6, p = 6, ...) {
  scores <- x[["scores"]]
  pheno <- x[["pheno"]]
  nr <- nrow(scores)
  nc <- ncol(scores)

  cat(sprintf("<mc_result> %d sample(s) x %d clock(s)\n", nr, nc))

  cat("\n<$pheno>")
  if (is.null(pheno)) {
    cat(" NULL\n")
  } else {
    pn <- min(n, nrow(pheno))
    cat(sprintf(" [%d of %d row(s)]\n", pn, nrow(pheno)))
    print(pheno[seq_len(pn), , drop = FALSE])
    if (pn < nrow(pheno)) {
      cat("...\n")
    }
  }

  ni <- min(n, nr)
  pi <- min(p, nc)
  cat(sprintf(
    "\n<$scores> [%d of %d row(s), %d of %d clock(s)]\n",
    ni,
    nr,
    pi,
    nc
  ))
  print(scores[seq_len(ni), seq_len(pi), drop = FALSE])
  if (ni < nr || pi < nc) {
    cat("...\n")
  }

  invisible(x)
}

# naked scores (coverage and provenance stay on the record)
#' @export
as.matrix.mc_result <- function(x, ...) {
  x[["scores"]]
}

# rbind refused -- re-run calc_clocks() on the combined DNAm
#' @export
rbind.mc_result <- function(..., deparse.level = 1) {
  cli::cli_abort(
    c(
      "Can't combine {.cls mc_result} records with {.fn rbind}.",
      "i" = "Please re-run {.fn calc_clocks} on the combined DNAm matrix so
             batch-dependent clocks see all samples at once."
    ),
    call = NULL
  )
}
