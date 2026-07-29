# scoring seam: mc_spec (data-independent) + mc_cohort (cohort facts) + score_cohort

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

# data-independent: resolved once, whatever the front end
mc_spec <- function(
  clocks,
  pheno_id = "ID",
  normalize = NULL,
  ext_data = NULL,
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

  packs <- load_mc_assets(pack_groups_needed(clock_sequence), ext_data, ask)

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

  mna <- scan_missing_cpgs(
    DNAm,
    panels_union(spec[["panels"]]),
    panels_union(spec[["panels"]], "score")
  )
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

# scoring-time failures per clock (coverage cannot see these)
new_notes <- function() {
  new.env(parent = emptyenv())
}

# record that `sample_id` could not be scored for clock `id`
note_scoring_failure <- function(block, id, sample_id) {
  notes <- block[["notes"]]
  if (!length(sample_id)) {
    return(invisible(NULL))
  }
  notes[[id]] <- union(notes[[id]], sample_id)
  invisible(NULL)
}

# collector -> plain named list (empty when nothing failed)
collect_notes <- function(notes) {
  ids <- sort(names(notes))
  if (!length(ids)) {
    return(list())
  }
  stats::setNames(lapply(ids, function(id) notes[[id]]), ids)
}

# per-block view: DNAm + its usable column index, partial cache, pheno, notes
mc_block <- function(DNAm, spec, facts) {
  usable <- facts[["usable_cols"]]
  # name -> column position, resolved once instead of per panel gather
  usable_idx <- match(usable, colnames(DNAm))
  if (anyNA(usable_idx)) {
    stop(
      sprintf(
        paste0(
          "mc_block: %d usable CpG(s) are not columns of this block. ",
          "This is a package bug -- please report it."
        ),
        sum(is.na(usable_idx))
      ),
      call. = FALSE
    )
  }

  block <- list(
    DNAm = DNAm,
    DNAm_full = DNAm,
    pheno = block_pheno(DNAm, facts),
    packs = spec[["packs"]],
    usable = usable,
    usable_idx = usable_idx,
    sample_id = rownames(DNAm),
    # write-only collector for scoring-time failures
    notes = new_notes()
  )
  fill <- facts[["partial_fill"]]
  block[["partial_cache"]] <- build_partial_cache(
    DNAm,
    block_cols(names(fill), block),
    fill
  )
  block
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

  # one score_type() per clock, reused by the pack filter and the dispatch
  types <- vapply(clock_sequence, score_type, character(1))
  is_pack <- unname(types) %in% PACK_SCORE_TYPES
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
    ty <- types[[p]]
    out <- switch(
      ty,
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
      stop(
        sprintf("No dispatch branch for score_type %s (clock %s).", ty, p),
        call. = FALSE
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
    ty <- score_type(p)
    scores[[p]] <- switch(
      ty,
      PhysAge = finalize_PhysAge(p, pending[[p]]),
      stop(
        sprintf("No finalize branch for score_type %s (clock %s).", ty, p),
        call. = FALSE
      )
    )
  }
  scores
}
