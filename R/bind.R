# rbind over mc_result records. record batching, refuse differing caller choices.

# one provenance field, flattened over every record
prov <- function(recs, field) {
  unlist(lapply(recs, function(r) r[["provenance"]][[field]]))
}

# gates
check_bind_inputs <- function(recs) {
  bad <- which(!vapply(recs, inherits, logical(1L), "mc_result"))
  if (length(bad)) {
    cli::cli_abort(
      c(
        "{.fn rbind} needs every argument to be an {.cls mc_result}.",
        # qty() on a length-1 numeric reads the value, not the count
        "x" = "Argument{cli::qty(length(bad))}{?s} {.val {capped_vals(bad)}}
               {cli::qty(length(bad))}{?is/are} a different class.",
        "i" = "{.fn calc_clocks} returns an {.cls mc_result}."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# a sample id may belong to exactly one batch
gate_disjoint_ids <- function(recs) {
  ids <- prov(recs, "sample_id")
  dup <- unique(ids[duplicated(ids)])
  if (length(dup)) {
    cli::cli_abort(
      c(
        "{length(dup)} sample id{?s} appear{?s/} in more than one
         {.cls mc_result}:",
        capped_bullets(dup, val_lines),
        "i" = "Give each batch its own sample ids before scoring.",
        "i" = "For example,
               {.code rownames(DNAm) <- paste0(rownames(DNAm), '_T1')}.",
        "i" = "The returned value must not count one sample twice."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# one provenance field, compared as a set against the first record
gate_same_set <- function(recs, field, what, hint) {
  ref <- recs[[1L]][["provenance"]][[field]]
  for (i in seq_along(recs)[-1L]) {
    got <- recs[[i]][["provenance"]][[field]]
    if (setequal(ref, got)) {
      next
    }
    diff <- c(setdiff(ref, got), setdiff(got, ref))
    cli::cli_abort(
      c(
        "Argument {i} has different {what} from argument 1:",
        capped_bullets(diff, val_lines),
        "i" = hint
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# pheno columns follow pheno_id and the clock sequence.
gate_same_pheno_id <- function(recs) {
  ref_id <- recs[[1L]][["provenance"]][["pheno_id"]]
  ids <- vapply(
    recs,
    function(r) r[["provenance"]][["pheno_id"]],
    character(1L)
  )
  if (!all(ids == ref_id)) {
    cli::cli_abort(
      c(
        "{.fn calc_clocks} scored these records with different
         {.arg pheno_id} columns: {.val {capped_vals(unique(ids))}}.",
        "i" = "Call {.fn calc_clocks} again with one {.arg pheno_id} column
               for every batch."
      ),
      call = NULL
    )
  }

  invisible(NULL)
}

run_bind_gates <- function(recs) {
  gate_disjoint_ids(recs)
  gate_same_set(
    recs,
    "clocks",
    "score columns",
    "Bind only records that {.fn calc_clocks} scored from the same
     {.arg clocks} request."
  )
  gate_same_pheno_id(recs)
  gate_same_set(
    recs,
    "normalized",
    "normalized clocks",
    "Use one {.arg normalize} setting for every batch. A clock that is
     normalized in one batch, and not in another, gives two different
     columns."
  )
  invisible(NULL)
}

# batch labels are a hash of the id set gate 1 made disjoint.

# assembly
stack_by_cols <- function(recs, get, cols) {
  # rbind carries the rownames through, including the zero-column norm panel
  do.call(rbind, lapply(recs, function(r) get(r)[, cols, drop = FALSE]))
}

# a clock-keyed provenance list, folded over records by `combine`
bind_by_key <- function(recs, field, combine) {
  each <- lapply(recs, function(r) r[["provenance"]][[field]])
  ids <- sort(unique(unlist(lapply(each, names))))
  if (!length(ids)) {
    return(list())
  }
  stats::setNames(
    lapply(ids, function(id) combine(lapply(each, function(e) e[[id]]))),
    ids
  )
}

# say when a bound record holds one reduction per batch, not one for the cohort
say_pending <- function(x) {
  ids <- names(x[["provenance"]][["pending"]])
  n_batch <- n_batches(x)
  if (!length(ids) || n_batch < 2L) {
    return(invisible(NULL))
  }
  cli::cli_inform(c(
    "!" = "The returned value holds {n_batch} separate values for
           {cli::qty(ids)}column{?s} {.val {ids}}, one value per batch.",
    "i" = "{cli::qty(ids)}{?This/These} score{?s} {cli::qty(ids)}{?is/are}
           computed from all the samples together, and {.fn calc_clocks}
           scored each batch on its own.",
    "i" = "Call {.fn refinalize_clocks} to compute {cli::qty(ids)}{?it/them}
           again from all {nrow(x[['scores']])} samples."
  ))
  invisible(NULL)
}

# bind into one labelled union. not a re-run
#' Combined Batches Of Scores
#'
#' Stacks two or more outputs from [calc_clocks()] runs into one object of
#' multiple batches.
#'
#' @param ... Two or more `mc_result` objects.
#' @param deparse.level A single whole number. Not used by this method.
#'   Default is `1`.
#'
#' @details
#' Each input must use disjoint sample ids, the same scored clocks, the
#' same `pheno_id`, and the same `normalize` setting. `rbind()` stops when
#' any of those differ between inputs.
#'
#' The combined value gets one `mc_batch_id` label for each input. A clock
#' that depends on sample-wise information, such as a z-score, keeps the
#' value each input calculated on its own samples. Call
#' [refinalize_clocks()] to calculate it again from every sample in the
#' combined value.
#'
#' @returns An `mc_result` object. It holds the stacked scores, `pheno`, and
#'   `coverage` of every input, under one `mc_batch_id` label for each.
#'
#' @seealso
#' - [refinalize_clocks()] for a cross-sample score recomputed after a bind.
#' - [as.data.frame.mc_result()] for the scores as a data.frame.
#' - [as.matrix.mc_result()] for the scores as a numeric matrix.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim1 <- sim_DNAm(clocks, n = 10)
#' sim2 <- sim_DNAm(clocks, n = 10, suffix = "b")
#'
#' res1 <- calc_clocks(sim1[["DNAm"]], clocks)
#' res2 <- calc_clocks(sim2[["DNAm"]], clocks)
#'
#' combined <- rbind(res1, res2)
#' combined
#'
#' @export
rbind.mc_result <- function(..., deparse.level = 1) {
  # names are dropped, not read. labels are derived
  args <- unname(list(...))
  args <- args[!vapply(args, is.null, logical(1L))]
  if (!length(args)) {
    return(NULL)
  }
  check_bind_inputs(args)
  run_bind_gates(args)

  # every gated field is the same across records, so record 1 speaks for all
  ref <- args[[1L]][["provenance"]]
  clocks <- ref[["clocks"]]

  scores <- stack_by_cols(args, function(r) r[["scores"]], clocks)
  # no rownames reset. resolve_pheno leaves automatic ones.
  pheno <- stack_by_cols(
    args,
    function(r) r[["pheno"]],
    names(args[[1L]][["pheno"]])
  )

  miss_panel <- function(panel) {
    stack_by_cols(
      args,
      function(r) r[["coverage"]][["sample_miss"]][[panel]],
      colnames(args[[1L]][["coverage"]][["sample_miss"]][[panel]])
    )
  }

  # honest union of what was asked for anywhere (columns already match)
  requested <- unique(prov(args, "requested"))

  out <- structure(
    list(
      scores = scores,
      pheno = pheno,
      coverage = list(
        # batch -> clock -> record, each under the label it was born with
        per_clock = unlist(
          lapply(args, function(r) r[["coverage"]][["per_clock"]]),
          recursive = FALSE
        ),
        sample_miss = list(
          score = miss_panel("score"),
          norm = miss_panel("norm")
        )
      ),
      provenance = list(
        sample_id = rownames(scores),
        mc_batch_id = prov(args, MC_BATCH),
        pheno_id = ref[["pheno_id"]],
        clocks = clocks,
        requested = requested,
        dependencies = setdiff(clocks, requested),
        covariates_used = ref[["covariates_used"]],
        normalized = ref[["normalized"]],
        # kept per batch. never reconciled.
        min_clocks_coverage = prov(args, "min_clocks_coverage"),
        min_samples_coverage = prov(args, "min_samples_coverage"),
        # clock -> the sample ids it failed on anywhere
        scoring_failures = bind_by_key(args, "scoring_failures", function(v) {
          unique(unlist(v))
        }),
        # intermediates stack by row, like the scores they will become
        pending = bind_by_key(args, "pending", function(v) do.call(rbind, v))
      )
    ),
    class = "mc_result"
  )
  say_pending(out)
  out
}

# recompute cohort-reducing clocks over every sample. leaves pending so calls compose.
#' Scores Recomputed From All Samples
#'
#' Recalculates every clock that depends on sample-wise information, such as
#' a z-score, using the samples now in `x`.
#'
#' @inheritParams mc-params
#'
#' @details
#' A clock of that kind is calculated once from every sample in a
#' `calc_clocks()` call, not one sample at a time. After `rbind()` combines
#' several such values, each score still holds the value its own input
#' calculated. `refinalize_clocks()` calculates it again from every sample in
#' `x`.
#'
#' `refinalize_clocks()` changes nothing, and returns `x` unchanged, when
#' `x` holds no clock of that kind.
#'
#' @returns An `mc_result` object. The same as `x`, with any score that is
#'   computed from all its samples together computed again from every
#'   sample in `x`.
#'
#' @seealso
#' - [rbind.mc_result()] for two runs combined into one object.
#' - [as.data.frame.mc_result()] for the scores as a data.frame.
#' - [as.matrix.mc_result()] for the scores as a numeric matrix.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim1 <- sim_DNAm(clocks, n = 10)
#' sim2 <- sim_DNAm(clocks, n = 10, suffix = "b")
#'
#' res1 <- calc_clocks(sim1[["DNAm"]], clocks)
#' res2 <- calc_clocks(sim2[["DNAm"]], clocks)
#' combined <- rbind(res1, res2)
#'
#' # a no-op here, because neither Horvath1 nor Hannum is scored from all
#' # samples together
#' refinalize_clocks(combined)
#'
#' @export
refinalize_clocks <- function(x) {
  check_mc_result(x)
  pending <- x[["provenance"]][["pending"]]
  if (!length(pending)) {
    cli::cli_inform(
      c(
        "i" = "{.fn refinalize_clocks} changes only clocks that are scored
               from all the samples together. {.arg x} has none."
      )
    )
    return(invisible(x))
  }

  done <- finalize_cross_sample(list(), pending)
  ids <- intersect(names(done), colnames(x[["scores"]]))
  for (id in ids) {
    col <- done[[id]]
    # match by name, not row order
    x[["scores"]][, id] <- col[rownames(x[["scores"]]), 1L]
  }
  cli::cli_inform(c(
    "v" = "{cli::qty(ids)}{?Column/Columns} {.val {ids}}
           {cli::qty(ids)}{?is/are} now computed from all
           {nrow(x[['scores']])} samples."
  ))
  invisible(x)
}
