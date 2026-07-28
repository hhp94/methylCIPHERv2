# The seam. calc_clocks() -- and any chunked front end -- composes these three,
# split by what each depends on: `spec` is decided before any DNAm is read,
# `facts` is everything only definable against the whole cohort, and
# score_cohort() sees a matrix plus those two. `facts` is the only channel
# between the front ends, so score_cohort() never learns which one called it.

# data-independent: resolved once, whatever the front end
mc_spec <- function(
  clocks,
  pheno_id = "ID",
  normalize = NULL,
  from = NULL,
  ask = TRUE
) {
  checkmate::assert_string(pheno_id)
  clock_ids <- resolve_clocks(clocks)
  clock_sequence <- resolve_clocks_sequence(clock_ids)
  normalize <- resolve_normalize(normalize, clock_sequence)
  # routed members are internal machinery: scored, but never a score column
  output_ids <- drop_routed_members(c(
    clock_ids,
    setdiff(clock_sequence, clock_ids)
  ))
  note_full_panel_clocks(clock_sequence)

  # the one covariate union: gates the pheno check, narrows the carried
  # pheno, and is stamped as provenance
  covariates <- unique(unlist(lapply(
    clock_sequence,
    clock_covariates_required
  )))
  if (is.null(covariates)) {
    covariates <- character(0)
  }

  packs <- load_mc_assets(pack_groups_needed(clock_sequence), from, ask)

  list(
    clock_ids = clock_ids,
    sequence = clock_sequence,
    output_ids = output_ids,
    normalize = normalize,
    covariates = covariates,
    pheno_id = pheno_id,
    packs = packs,
    panels = clock_panels(clock_sequence, packs, normalize),
    # sample-axis split, declared by the catalog. `cross_sample` names the
    # clocks whose reduction is still inside their branch, so scoring a row
    # subset does not reproduce their whole-cohort rows -- the list Phase 3
    # empties by emitting per-sample intermediates and finalizing after
    # assembly. Derived here so no caller carries a clock list of its own.
    cross_sample = split_cross_sample(clock_sequence)[["cross_sample"]]
  )
}

# cohort-set facts, plus the gates that must fire before anything is scored.
# One scan today; a chunked front end accumulates the same fields over blocks.
mc_cohort <- function(DNAm, spec, pheno = NULL, min_clocks_coverage = 0.75) {
  if (length(spec[["covariates"]]) && is.null(pheno)) {
    cli::cli_abort(
      c(
        "These clocks need {cli::qty(spec[['covariates']])} pheno column{?s}
         {.field {spec[['covariates']]}}, but {.arg pheno} is missing.",
        "i" = "Pass a pheno table with {cli::qty(spec[['covariates']])}
               {?that/those} column{?s}."
      ),
      call = NULL
    )
  }

  # sample ids are the DNAm rownames -- mandatory, enforced by check_DNAm().
  check_DNAm(DNAm)
  sample_id <- rownames(DNAm)

  check_pheno(
    pheno,
    ID = spec[["pheno_id"]],
    extra_columns = spec[["covariates"]],
    sample_id = sample_id
  )
  pheno <- resolve_pheno(DNAm, pheno, spec[["pheno_id"]], spec[["covariates"]])

  mna <- scan_missing_cpgs(DNAm, panels_union(spec[["panels"]]))
  cpg_list <- resolve_cpgs(mna[["usable_cols"]], spec[["panels"]])
  check_coverage(cpg_list, min_clocks_coverage)

  list(
    sample_id = sample_id,
    pheno = pheno,
    usable_cols = mna[["usable_cols"]],
    cpg_list = cpg_list,
    # the names are the column classification, not just a lookup of fill
    # values -- a block must never re-derive which columns are cohort-partial.
    # Every partial column is by construction a usable panel column, so this is
    # already exactly the set the panels can reach.
    partial_fill = mna[["col_mean"]]
  )
}

# pheno rows follow facts$sample_id; narrow to the rows in hand
block_pheno <- function(DNAm, facts) {
  if (is.null(facts[["pheno"]])) {
    return(NULL)
  }
  facts[["pheno"]][match(rownames(DNAm), facts[["sample_id"]]), , drop = FALSE]
}

# Everything a branch may read about the rows in hand, built once per block.
# Branches take `(id, cpgs, block, results)` and nothing else, so a new cohort
# fact reaches all of them by being added here rather than by re-threading it
# through every signature between here and the leaf that wants it.
#
# `DNAm` is narrowed to the panel union up front. Every branch otherwise gathers
# its own panel out of the caller's full-width array, once per clock, and that
# strided read grows with the array while the panel does not. `DNAm_full` is the
# caller's matrix untouched, read only by a clock whose recipe z-scores each
# sample over every probe it was handed (Zhang2019) -- narrowing that input
# would move its moments and its score.
mc_block <- function(DNAm, spec, facts) {
  narrowed <- DNAm[, facts[["usable_cols"]], drop = FALSE]
  list(
    DNAm = narrowed,
    DNAm_full = DNAm,
    partial_cache = build_partial_cache(narrowed, facts[["partial_fill"]]),
    pheno = block_pheno(DNAm, facts),
    packs = spec[["packs"]],
    usable = facts[["usable_cols"]],
    sample_id = rownames(DNAm),
    # write-only collector for scoring-time failures; see new_notes()
    notes = new_notes()
  )
}

# a matrix plus the cohort facts -> score matrices, per-sample intermediates
# for the cohort-reducing clocks, the coverage fragment for those rows, and any
# scoring-time failures. Returns no record: assembly is the front end's job.
score_cohort <- function(DNAm, spec, facts) {
  clock_sequence <- spec[["sequence"]]
  cpg_list <- facts[["cpg_list"]]
  block <- mc_block(DNAm, spec, facts)

  # coverage/QC needs no score: compute it once, keyed by clock id
  coverage <- compute_coverage(clock_sequence, cpg_list, block)

  # scoring loop returns score matrices only -- coverage was computed above
  results <- vector("list", length(clock_sequence))
  names(results) <- clock_sequence
  # per-sample intermediates for the cohort-reducing clocks; their slot in
  # `results` stays empty until finalize_cross_sample() fills it
  pending <- list()

  is_pack <- vapply(clock_sequence, is_pack_scored, logical(1))
  if (any(is_pack)) {
    pack_ids <- clock_sequence[is_pack]
    pgroups <- vapply(pack_ids, clock_group_id, character(1))
    for (g in unique(pgroups)) {
      grp <- score_pack_group(pack_ids[pgroups == g], block)
      results[names(grp)] <- grp
    }
  }

  # one line per tag: every branch reads the same `(id, cpgs, block, results)`,
  # so the switch dispatches rather than re-typing arguments. `linear_score()`
  # is the shared engine, not only a branch -- it keeps the engine signature
  # because other branches call it too.
  for (p in clock_sequence[!is_pack]) {
    cpgs <- cpg_list[["per_clock"]][[p]]
    out <- switch(
      score_type(p),
      linear = linear_score(cpgs, block),
      GrimAge = score_GrimAge(p, cpgs, block, results),
      DNAmFitAge = score_DNAmFitAge(p, cpgs, block, results),
      PhysAge = physage_raws(p, cpgs, block, results),
      Dunedin = score_Dunedin(p, cpgs, block, results),
      normalized = score_normalized(p, cpgs, block, results),
      EpiTOC2 = score_EpiTOC2(p, cpgs, block, results),
      MiAge = score_MiAge(p, cpgs, block, results),
      Zhang2019 = score_Zhang2019(p, cpgs, block, results),
      sex_routed = score_sex_routed(p, cpgs, block, results),
      cli::cli_abort(
        "No dispatch branch for score_type {.val {score_type(p)}}
         (clock {.val {p}}).",
        call = NULL
      )
    )
    # a cohort-reducing clock yields its per-sample intermediate here, not a
    # score. Which clocks those are is the catalog's declaration, never a
    # clock list and never the branch's own choice.
    if (p %in% spec[["cross_sample"]]) {
      pending[[p]] <- out
    } else {
      results[[p]] <- out
    }
  }

  list(
    scores = results,
    coverage = coverage,
    pending = pending,
    # per-clock sample ids the branch could not score, for these rows only
    notes = collect_notes(block[["notes"]])
  )
}

# The one cohort-reduction point, called unconditionally by every front end:
# a single-pass run finalizes its own block, a chunked run finalizes the
# assembled intermediates. `pending` is empty for all but the cohort-reducing
# clocks, so this is a no-op for a typical request and nothing downstream
# branches on whether it was chunked.
finalize_cross_sample <- function(scores, pending) {
  for (p in names(pending)) {
    scores[[p]] <- switch(
      score_type(p),
      PhysAge = finalize_PhysAge(p, pending[[p]]),
      cli::cli_abort(
        "No finalize branch for score_type {.val {score_type(p)}}
         (clock {.val {p}}).",
        call = NULL
      )
    )
  }
  scores
}
