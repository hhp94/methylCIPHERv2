# scoring seam: mc_spec (data-independent) + mc_cohort (cohort facts) + score_cohort

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
  # routed members are scored but never a score column
  output_ids <- drop_routed_members(c(
    clock_ids,
    setdiff(clock_sequence, clock_ids)
  ))
  note_full_panel_clocks(clock_sequence)

  # covariate union for pheno check, carried pheno, and provenance
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
    # clocks whose reduction is still inside the branch (catalog-declared)
    cross_sample = split_cross_sample(clock_sequence)[["cross_sample"]]
  )
}

# cohort-set facts + pre-score gates (chunked front end accumulates these)
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

  # sample ids are the DNAm rownames
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
    # partial_fill names are the column classification (do not re-derive)
    partial_fill = mna[["col_mean"]]
  )
}

# pheno rows follow facts$sample_id, narrowed to the rows in hand
block_pheno <- function(DNAm, facts) {
  if (is.null(facts[["pheno"]])) {
    return(NULL)
  }
  facts[["pheno"]][match(rownames(DNAm), facts[["sample_id"]]), , drop = FALSE]
}

# per-block view: narrowed DNAm, full DNAm, partial cache, pheno, notes
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
    # write-only collector for scoring-time failures
    notes = new_notes()
  )
}

# score one block: scores, coverage, pending intermediates, notes
score_cohort <- function(DNAm, spec, facts) {
  clock_sequence <- spec[["sequence"]]
  cpg_list <- facts[["cpg_list"]]
  block <- mc_block(DNAm, spec, facts)

  # coverage before scoring, keyed by clock id
  coverage <- compute_coverage(clock_sequence, cpg_list, block)

  results <- vector("list", length(clock_sequence))
  names(results) <- clock_sequence
  # per-sample intermediates for cohort-reducing clocks
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

  # branch dispatch: every scorer takes (id, cpgs, block, results)
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
    # cohort-reducing clocks yield intermediates into pending
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
    # per-clock sample ids the branch could not score
    notes = collect_notes(block[["notes"]])
  )
}

# cohort reduction after assembly (no-op when pending is empty)
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
