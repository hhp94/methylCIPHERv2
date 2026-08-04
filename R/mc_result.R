# mc_result record: constructor + methods for class "mc_result"

# batch label is a hash of the pheno id column.
batch_hash <- function(ids) {
  # hash the id set, not the sequence or r representation
  key <- paste0(
    sort(unname(as.character(ids)), method = "radix"),
    collapse = "\r"
  )
  digest::digest(key, algo = "xxhash64", serialize = FALSE)
}

# the record's distinct batch labels.
batch_labels <- function(x) {
  unique(x[["provenance"]][[MC_BATCH]])
}

# true when the record spans more than one batch.
is_multi_batch <- function(batch) {
  length(unique(batch)) > 1L
}

# how many batches the record spans. the per-sample provenance vector is
# authoritative -- it is the vector that fills mc_batch_id -- and per_clock
# must agree with it. a disagreement means the record was assembled wrong, so
# it stops rather than silently picking one of the two counts.
n_batches <- function(x) {
  n <- length(unique(x[["provenance"]][[MC_BATCH]]))
  n_cov <- length(x[["coverage"]][["per_clock"]])
  if (n != n_cov) {
    stop(
      sprintf(
        paste0(
          "mc_result batch count disagrees: provenance has %d, coverage has ",
          "%d. This is a package bug -- please report it."
        ),
        n,
        n_cov
      ),
      call. = FALSE
    )
  }
  n
}

# drop the batch column when it is a single repeated hash.
drop_single_batch <- function(df, batch) {
  if (is_multi_batch(batch)) {
    return(df)
  }
  # no-op where the column was never built.
  df[[MC_BATCH]] <- NULL
  df
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
  min_clocks_coverage,
  min_samples_coverage,
  scoring_failures = list(),
  pending = list()
) {
  scores <- do.call(cbind, results[output_ids])
  dimnames(scores) <- list(sample_id, output_ids)

  # derived, never passed in: the id column is the batch's whole identity
  batch <- batch_hash(pheno[[pheno_id]])

  # coverage spans clocks that read cpgs.
  per_clock <- coverage[["per_clock"]]
  record_ids <- covered_ids(per_clock)
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
      coverage = list(
        # batch -> clock -> record. one fill regime per batch
        per_clock = stats::setNames(list(per_clock), batch),
        sample_miss = sample_miss
      ),
      provenance = list(
        sample_id = sample_id,
        # per-sample, aligned to sample_id (never a key -- see rbind)
        mc_batch_id = rep(batch, length(sample_id)),
        pheno_id = pheno_id,
        clocks = output_ids,
        requested = requested_ids,
        dependencies = setdiff(output_ids, requested_ids),
        covariates_used = covariates_used,
        # which clocks were actually normalized
        normalized = normalized,
        # the gates this batch was scored under, keyed by batch like per_clock
        min_clocks_coverage = stats::setNames(min_clocks_coverage, batch),
        min_samples_coverage = stats::setNames(min_samples_coverage, batch),
        # clock id -> sample ids the scorer could not fit
        scoring_failures = scoring_failures,
        # retained per-sample intermediates, so a bind can re-finalize exactly
        pending = pending
      )
    ),
    class = "mc_result"
  )
}

# scores then pheno, in the shared printer grammar (R/print.R)
#' Print Method For An mc_result Object
#'
#' Prints the scores and the `pheno` table for an `mc_result` object.
#'
#' @inheritParams mc-params
#' @param n A single whole number. The number of sample rows to print.
#'   Default is `6`.
#' @param p A single whole number. The number of clock columns to print for
#'   the scores table. Default is `6`.
#' @param ... Not used.
#'
#' @details
#' The output lists batch labels only when `x` spans more than one
#' `mc_batch_id`.
#'
#' @returns An `mc_result` object. Returns `x`, invisibly, after printing it.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 20)
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#' print(res, n = 3, p = 2)
#'
#' @export
print.mc_result <- function(x, n = 6, p = 6, ...) {
  scores <- x[["scores"]]
  pheno <- x[["pheno"]]

  cat(
    fmt_header("mc_result", nrow(scores), "sample", ncol(scores), "clock"),
    "\n",
    sep = ""
  )
  print_block(
    "scores",
    scores,
    min(n, nrow(scores)),
    min(p, ncol(scores)),
    "clock"
  )
  # always present -- the id column at minimum (see resolve_pheno)
  print_block(
    "pheno",
    pheno,
    min(n, nrow(pheno)),
    ncol(pheno),
    "column",
    cut_cols = FALSE
  )

  # multi-batch only, same test as the exit frames.
  labels <- batch_labels(x)
  if (length(labels) > 1L) {
    cat(
      "\n",
      fmt_section("provenance", plural_count(length(labels), "batch", "es")),
      "\n",
      paste(labels, collapse = ", "),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

# naked scores (coverage and provenance stay on the record)
#' Matrix Method For An mc_result Object
#'
#' Converts the [calc_clocks()] output to a matrix containing just the clocks.
#'
#' @inheritParams mc-params
#' @param ... Not used.
#'
#' @details
#' This function recalculates any clock that depends on sample-wise
#' information, such as a z-score, from all the available samples when `x`
#' holds more than one batch. This is the same calculation as
#' [refinalize_clocks()].
#'
#' @returns A numeric matrix. The scores, with samples in the rows and
#'   clocks in the columns.
#'
#' @seealso
#' - [as.data.frame.mc_result()] for the scores as a data.frame.
#' - [rbind.mc_result()] for two runs combined into one object.
#' - [refinalize_clocks()] for a cross-sample score recomputed after a bind.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 10)
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#' as.matrix(res)
#'
#' @export
as.matrix.mc_result <- function(x, ...) {
  check_mc_result(x)
  # a finalizer: the matrix cannot carry `pending`, so resolve it on the way out
  finalized(x)[["scores"]]
}
