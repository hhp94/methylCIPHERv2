# cohort-gated golden parity (needs MC_PARITY=1 and a staged cohort)

# cohort fixture access
meta_clone_path <- function(...) {
  testthat::test_path("..", "..", "data-raw", "methylCIPHER-meta", ...)
}

# registry cohorts (scripts/cohorts.py). paths derive from the id.
PARITY_COHORTS <- c("cohort_EPICv1", "cohort_450K")

cohort_beta_db <- function(cohort) {
  meta_clone_path("fixtures", cohort, "beta.duckdb")
}

# samples x CpGs from the tall beta table
read_betas <- function(con, sql) {
  raw <- DBI::dbGetQuery(con, sql)
  mat <- as.matrix(raw[, setdiff(names(raw), "cpg"), drop = FALSE])
  rownames(mat) <- raw$cpg
  rm(raw)
  t(mat)
}

cohort_betas <- function(con, cpgs) {
  read_betas(
    con,
    sprintf(
      "SELECT * FROM beta WHERE cpg IN (%s)",
      paste0("'", cpgs, "'", collapse = ",")
    )
  )
}

# whole array for clocks whose recipe moments over the full input
cohort_betas_full <- function(con) {
  read_betas(con, "SELECT * FROM beta")
}

# cohort pheno: id, Tissue, Age, Female
cohort_pheno <- function(cohort) {
  ph <- utils::read.csv(
    meta_clone_path("fixtures", cohort, "pheno.csv"),
    stringsAsFactors = FALSE
  )
  data.frame(
    ID = ph$id,
    Age = ph$Age,
    Female = as.integer(ph$Female),
    stringsAsFactors = FALSE
  )
}

# golden scores from the clock's fixture block for this cohort
expected_scores <- function(id, cohort) {
  rel <- clock_fixture(id, cohort)[["expected"]]
  if (is.null(rel)) {
    stop(
      "Clock '",
      id,
      "' has no fixture expected path for ",
      cohort,
      call. = FALSE
    )
  }
  utils::read.csv(gzfile(meta_clone_path(rel)), stringsAsFactors = FALSE)
}

# both max_abs and max_rel must clear. only abs tol varies by block (units).
PARITY_REL_TOL <- 1e-10

# packs: abs tol 1e-6 (scale ~3e6).
PARITY_ABS_TOL <- c(core = 1e-10, fitage = 1e-10, packs = 1e-6, horvath = 1e-10)

# snapshot of the normalized-horvath residual (not an agreement target). keyed clock@cohort.
#   Horvath1@cohort_450K  max_abs 1.137877e-01  max_rel 1.926282e-03
# a pair with no entry here fails. measure new pairs. do not default.
HORVATH_NORM_TOL <- list(
  "Horvath1@cohort_450K" = c(abs = 1.2e-1, rel = 2.0e-3)
)

# horvath-online oracle clocks (declared, not listed)
is_horvath_online <- function(id) {
  any(vapply(
    clock_fixtures(id) %||% list(),
    function(fx) identical(as.character(fx[["oracle"]]), "horvath_online"),
    logical(1)
  ))
}

# horvath-online clocks with an expressible scheme. today only Horvath1 (bmiq).
is_normalized_horvath <- function(id) {
  is_horvath_online(id) && clock_norm_scheme(id) %in% NORM_SCHEMES
}

# which block a clock is graded in. Derived from the catalog, never a clock list.
parity_block <- function(id) {
  if (is_horvath_online(id)) {
    "horvath"
  } else if (is_pack_scored(id)) {
    # the relaxed tolerance follows the pack scoring path, not externality
    "packs"
  } else if (identical(clock_group_id(id), "DNAmFitAge")) {
    "fitage"
  } else {
    "core"
  }
}

# max, not median (median hides bad samples)
rel_diff <- function(got, want) {
  max(abs(got - want) / pmax(abs(want), .Machine$double.eps))
}

expect_parity <- function(got, id, cohort) {
  exp <- expected_scores(id, cohort)
  aligned <- as.numeric(got[exp$sample_id])
  testthat::expect_false(
    anyNA(aligned),
    label = paste0(id, "/", cohort, ": scored samples missing for fixture ids")
  )
  testthat::expect_lt(
    max(abs(aligned - exp$value)),
    PARITY_ABS_TOL[[parity_block(id)]],
    label = sprintf("%s/%s max_abs_diff", id, cohort)
  )
  testthat::expect_lt(
    rel_diff(aligned, exp$value),
    PARITY_REL_TOL,
    label = sprintf("%s/%s max_rel_diff", id, cohort)
  )
}

# parity tier flag (gates duckdb, pack scan, and per-test skips)
parity_on <- nzchar(Sys.getenv("MC_PARITY"))

# cached external packs (empty when tier is off)
cached_pack_groups <- if (parity_on) {
  Filter(function(g) length(mc_staged_files(g)) > 0L, mc_external_groups())
} else {
  character(0)
}

# skip external clocks whose pack is not cached
skip_if_no_pack <- function(clock_id) {
  if (!clock_is_external(clock_id)) {
    return(invisible())
  }
  gid <- clock_group_id(clock_id)
  testthat::skip_if_not(
    gid %in% cached_pack_groups,
    paste0("external pack for '", gid, "' not cached")
  )
}

# one read-only duckdb connection per staged cohort, for this file
cohort_cons <- list()
if (
  parity_on &&
    requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)
) {
  # duckdb extensions in a throwaway temp dir.
  withr::local_options(
    list(
      duckdb.extension_directory = withr::local_tempdir(
        .local_envir = testthat::teardown_env()
      )
    ),
    .local_envir = testthat::teardown_env()
  )
  for (cohort in PARITY_COHORTS) {
    if (!file.exists(cohort_beta_db(cohort))) {
      next
    }
    con <- DBI::dbConnect(
      duckdb::duckdb(),
      cohort_beta_db(cohort),
      read_only = TRUE
    )
    cohort_cons[[cohort]] <- con
    local({
      cc <- con
      withr::defer(
        try(DBI::dbDisconnect(cc, shutdown = TRUE), silent = TRUE),
        envir = testthat::teardown_env()
      )
    })
  }
}

# cohorts actually on disk. generators key on this. unstaged tiers emit nothing.
staged_cohorts <- PARITY_COHORTS[PARITY_COHORTS %in% names(cohort_cons)]
unstaged_cohorts <- if (parity_on) {
  setdiff(PARITY_COHORTS, staged_cohorts)
} else {
  character(0)
}

# wang@cohort_450K: no sex-chromosome probes. skip. epicv1 passes.
WANG_450K_GAP <- paste0(
  "cohort_450K has no sex-chromosome probes, so the panel is 0% present. ",
  "The fixture expects the oracle's empty-panel 0; we decline to score."
)

# known gaps (clock- or clock@cohort-keyed).
KNOWN_PARITY_GAPS <- c(
  "DNAmSex_Wang_ChrX@cohort_450K" = WANG_450K_GAP,
  "DNAmSex_Wang_ChrY@cohort_450K" = WANG_450K_GAP
)

# horvath block skipped: oracle filled absent probes server-side.
HORVATH_ONLINE_GAP <- paste0(
  "horvath_online oracle -- server-side fill of absent probes is unpublished. ",
  "Pairs with no absent probes already match to ~1e-8"
)

# group-keyed gaps (separate map: group ids share namespace with clock ids)
KNOWN_PARITY_GAP_GROUPS <- character(0)

parity_gap <- function(id, cohort) {
  fx <- clock_fixture(id, cohort)
  if (identical(as.character(fx[["oracle"]] %||% NA), "horvath_online")) {
    return(HORVATH_ONLINE_GAP)
  }
  key <- paste0(id, "@", cohort)
  if (key %in% names(KNOWN_PARITY_GAPS)) {
    return(KNOWN_PARITY_GAPS[[key]])
  }
  if (id %in% names(KNOWN_PARITY_GAPS)) {
    return(KNOWN_PARITY_GAPS[[id]])
  }
  gid <- clock_group_id(id)
  if (!is.null(gid) && gid %in% names(KNOWN_PARITY_GAP_GROUPS)) {
    return(KNOWN_PARITY_GAP_GROUPS[[gid]])
  }
  NULL
}

# (clock, cohort) pairs with an upstream fixture. unstaged cohorts are not targets.
parity_targets <- function(block) {
  if (!length(staged_cohorts)) {
    return(list())
  }
  out <- list()
  for (id in names(mc_catalog)) {
    if (!identical(parity_block(id), block)) {
      next
    }
    for (fx in clock_fixtures(id) %||% list()) {
      cohort <- as.character(fx[["cohort"]])
      if (!cohort %in% staged_cohorts) {
        next
      }
      out[[length(out) + 1L]] <- list(id = id, cohort = cohort)
    }
  }
  out
}

run_parity_target <- function(clock_id, cohort) {
  skip_if_no_pack(clock_id)
  gap <- parity_gap(clock_id, cohort)
  if (!is.null(gap)) {
    skip(paste0("known parity gap -- ", gap))
  }
  # routed members scored as their alias's dependency
  routed <- sex_routed_members()$alias
  request <- if (clock_id %in% names(routed)) {
    routed[[clock_id]]
  } else {
    clock_id
  }
  # packs carry their group's scoring panel -- resolve before the union
  seq_ids <- resolve_clocks_sequence(resolve_clocks(request))
  packs <- load_mc_assets(pack_groups_needed(seq_ids), NULL, FALSE)
  DNAm <- if (any(vapply(seq_ids, clock_needs_full_panel, logical(1)))) {
    cohort_betas_full(cohort_cons[[cohort]])
  } else {
    # same union as clock_cpgs() (panels alone would drop moment refs).
    cohort_betas(cohort_cons[[cohort]], sequence_cpgs(seq_ids, packs))
  }
  # parity gates numbers, not coverage policy
  res <- calc_clocks(
    DNAm,
    request,
    pheno = cohort_pheno(cohort),
    ext_data = packs,
    min_clocks_coverage = 0,
    min_samples_coverage = 0
  )
  # routed member scores land on the alias column for that sex's samples
  expect_parity(res$scores[, request], clock_id, cohort)
}

# what the tier could not run, said once each rather than once per target
if (!parity_on) {
  test_that("cohort parity tier", {
    skip("parity tier off (set MC_PARITY=1, e.g. via dev test_parity())")
  })
}
for (cohort_i in unstaged_cohorts) {
  local({
    cohort <- cohort_i
    test_that(paste0("cohort parity: ", cohort), {
      skip(paste0(cohort, " fixture not staged"))
    })
  })
}

# census: every catalog clock declares a fixture per cohort. needs the flag, not duckdb.
if (parity_on) {
  test_that("every clock declares a fixture for every registry cohort", {
    ids <- names(mc_catalog)
    declared <- lapply(stats::setNames(ids, ids), function(id) {
      unique(vapply(
        clock_fixtures(id) %||% list(),
        function(fx) as.character(fx[["cohort"]]),
        character(1)
      ))
    })

    # sex-routed aliases declare no fixture (members do)
    aliases <- unique(unlist(sex_routed_members()$alias))
    expect_setequal(names(Filter(function(x) !length(x), declared)), aliases)

    # every other clock: both cohorts. Partial coverage is a gap, not a pass.
    rest <- setdiff(ids, aliases)
    incomplete <- rest[
      !vapply(
        rest,
        function(id) setequal(declared[[id]], PARITY_COHORTS),
        logical(1)
      )
    ]
    expect_equal(incomplete, character(0))

    # a cohort outside the registry would generate targets that can never run
    expect_setequal(unique(unlist(declared)), PARITY_COHORTS)
  })
}

# the four blocks and their labels. membership is derived by parity_block().
PARITY_BLOCK_LABELS <- c(
  core = "parity", # bundled clocks outside the DNAmFitAge family
  fitage = "parity (fitage)", # DNAmFitAge family at core tolerances
  packs = "parity (packs)", # external packs (relaxed abs tol only)
  horvath = "parity (horvath online)" # online oracles (skipped wholesale)
)

for (block_i in names(PARITY_BLOCK_LABELS)) {
  for (target_i in parity_targets(block_i)) {
    local({
      label <- PARITY_BLOCK_LABELS[[block_i]]
      clock_id <- target_i$id
      cohort <- target_i$cohort
      test_that(paste0(label, ": ", clock_id, " @ ", cohort), {
        run_parity_target(clock_id, cohort)
      })
    })
  }
}

# both PhysAge composites in one call, per cohort
for (cohort_i in staged_cohorts) {
  local({
    cohort <- cohort_i
    test_that(
      paste0("PhysAge composites match the author fixtures @ ", cohort),
      {
        members <- mc_groups[["PhysAge"]]$members
        cpgs <- unique(unlist(lapply(members, clock_scoring_cpgs)))
        DNAm <- cohort_betas(cohort_cons[[cohort]], cpgs)
        res <- calc_clocks(
          DNAm,
          c("DNAmPhysAge", "DNAmPhysAge_years"),
          pheno = cohort_pheno(cohort),
          min_clocks_coverage = 0,
          min_samples_coverage = 0
        )
        expect_parity(res$scores[, "DNAmPhysAge"], "DNAmPhysAge", cohort)
        expect_parity(
          res$scores[, "DNAmPhysAge_years"],
          "DNAmPhysAge_years",
          cohort
        )
      }
    )
  })
}

# normalized horvath-online (Horvath1 bmiq). skip unless the cohort leaves zero scoring probes absent.
for (cohort_i in staged_cohorts) {
  for (id_i in Filter(is_normalized_horvath, names(mc_catalog))) {
    local({
      cohort <- cohort_i
      clock_id <- id_i
      test_that(
        paste0("parity (horvath normalized): ", clock_id, " @ ", cohort),
        {
          skip_if_not_installed("betanorm")
          norm_on <- stats::setNames(TRUE, clock_id)
          DNAm <- cohort_betas(
            cohort_cons[[cohort]],
            clock_cpgs(clock_id, normalize = norm_on)
          )

          # any absent scoring probe confounds the comparison with the oracle fill.
          absent <- setdiff(clock_scoring_cpgs(clock_id), colnames(DNAm))
          skip_if_not(
            length(absent) == 0L,
            paste0(
              length(absent),
              " scoring CpGs absent -- ",
              HORVATH_ONLINE_GAP
            )
          )

          res <- calc_clocks(
            DNAm,
            clock_id,
            pheno = cohort_pheno(cohort),
            normalize = norm_on,
            min_clocks_coverage = 0,
            min_samples_coverage = 0
          )
          exp <- expected_scores(clock_id, cohort)
          got <- as.numeric(res$scores[, clock_id][exp$sample_id])
          expect_false(anyNA(got))

          key <- paste0(clock_id, "@", cohort)
          tol <- HORVATH_NORM_TOL[[key]]
          if (is.null(tol)) {
            stop("No recorded residual snapshot for ", key, call. = FALSE)
          }
          expect_lt(
            max(abs(got - exp$value)),
            tol[["abs"]],
            label = sprintf("%s max_abs_diff", key)
          )
          expect_lt(
            rel_diff(got, exp$value),
            tol[["rel"]],
            label = sprintf("%s max_rel_diff", key)
          )
        }
      )
    })
  }
}

# degraded coverage against the DunedinPACE reference. needs the tier flag, not duckdb.
if (parity_on) {
  test_that("DunedinPACE matches danbelsky/DunedinPACE through a holed panel", {
    skip_if_not_installed("betanorm")

    norm_panel <- names(clock_norm_target("DunedinPACE"))
    score_panel <- clock_scoring_cpgs("DunedinPACE")
    norm_only <- setdiff(norm_panel, score_panel)

    # complete miss: 3 of 173 scoring CpGs and 1000 background-only ones
    absent <- c(score_panel[1:3], norm_only[1:1000])
    DNAm <- random_betas(setdiff(norm_panel, absent), n = 10L)

    # partial miss on both panels. reference vendor-fills rare probes, we never do
    holed_score <- score_panel[4:8]
    holed_norm <- norm_only[1001:1005]
    for (j in seq_along(holed_score)) {
      DNAm[j, holed_score[[j]]] <- NA_real_
      DNAm[j, holed_norm[[j]]] <- NA_real_
    }

    got <- calc_clocks(DNAm, "DunedinPACE")
    ref <- DunedinPACE::PACEProjector(t(DNAm))[["DunedinPACE"]]
    expect_false(anyNA(ref))
    expect_equal(got$scores[, "DunedinPACE"], ref[rownames(DNAm)])

    # absent cpgs filled from the target. present-but-holed ones cohort-mean filled.
    cov <- got$coverage$per_clock[[1]]$DunedinPACE
    expect_equal(cov$score_imputed_full, 3L)
    expect_equal(cov$norm_imputed_full, 1003L)
    expect_equal(cov$score_imputed_partial, 5L)
    expect_equal(cov$norm_imputed_partial, 10L)
    expect_equal(cov$score_dropped, 0L)
    expect_equal(cov$norm_dropped, 0L)
  })
}
