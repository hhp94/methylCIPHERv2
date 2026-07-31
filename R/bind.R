# rbind over mc_result records. record what batching forced, refuse what the caller chose differently.

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
        "x" = "Argument{cli::qty(length(bad))}{?s} {.val {bad}}
               {cli::qty(length(bad))}{?is/are} not."
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
        "{length(dup)} sample id{?s} appear{?s/} in more than one record:",
        capped_bullets(dup),
        "i" = "Give each batch its own ids before scoring -- e.g.
               {.code rownames(DNAm) <- paste0(rownames(DNAm), '_T1')} -- so a
               bound record cannot double-count a sample."
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
        "Record {i} has different {what} from record 1.",
        capped_bullets(diff),
        "i" = hint
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# pheno columns need no gate. they follow pheno_id and the clock sequence
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
        "Records were scored against different {.arg pheno_id} columns:
         {.val {unique(ids)}}.",
        "i" = "Re-score with one id column so the bound pheno has one identity."
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
    "Bind records scored from the same {.arg clocks} request."
  )
  gate_same_pheno_id(recs)
  gate_same_set(
    recs,
    "normalized",
    "normalized clocks",
    "Use one {.arg normalize} setting across every batch -- the same clock
     normalized in one batch and not another is not one column."
  )
  invisible(NULL)
}

# no gate on batch labels: they are a hash of the ids gate 1 just made disjoint

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
  n_batch <- length(x[["coverage"]][["per_clock"]])
  if (!length(ids) || n_batch < 2L) {
    return(invisible(NULL))
  }
  cli::cli_inform(c(
    "!" = "Cohort-reducing {cli::qty(ids)}column{?s} {.val {ids}}
           {cli::qty(ids)}{?was/were} reduced within each batch, so this record
           holds {n_batch} separate reductions rather than one.",
    "i" = "Call {.fn refinalize_clocks} to recompute {cli::qty(ids)}{?it/them}
           against all {nrow(x[['scores']])} samples."
  ))
  invisible(NULL)
}

# bind into one labelled union. not a re-run
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
  # no rownames reset. resolve_pheno leaves automatic ones, and rbind renumbers them
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
        batch = prov(args, "batch"),
        pheno_id = ref[["pheno_id"]],
        clocks = clocks,
        requested = requested,
        dependencies = setdiff(clocks, requested),
        covariates_used = ref[["covariates_used"]],
        normalized = ref[["normalized"]],
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

# recompute cohort-reducing clocks over every sample. never automatic, and it leaves pending so calls compose
#' @export
refinalize_clocks <- function(x) {
  check_mc_result(x)
  pending <- x[["provenance"]][["pending"]]
  if (!length(pending)) {
    cli::cli_inform(
      c(
        "i" = "No cohort-reducing clocks on this record; nothing to re-finalize."
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
    "v" = "Re-finalized {cli::qty(ids)}{?column/columns} {.val {ids}} against
           all {nrow(x[['scores']])} samples."
  ))
  invisible(x)
}
