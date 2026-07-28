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

# the cohort's whole array (473k-866k probes, ~350-550 MB). Only for clocks whose
# recipe takes moments over the input matrix -- see needs_full_panel().
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

# Every fixture is graded on both axes and must clear both. The two are not
# redundant: the absolute bound is the only one with meaning near zero, and the
# relative bound is the only one with meaning at large magnitude.
#
# Only the absolute bound is scale-sensitive, so only it varies by block. The
# relative bound is scale-free and is 1e-10 everywhere, with no exceptions --
# relaxing a block is a statement about units, never about correctness.
PARITY_REL_TOL <- 1e-10

# External-pack biomarkers reach ~3.3e6 (PCB2M), where 1e-10 absolute sits below
# the float floor before any arithmetic. Measured worst absolute miss over all 56
# pack rows is 2.05e-8, at 6.3e-15 relative -- noise, not error. 1e-6 clears that
# by ~50x and is still ~3e-13 relative at PCB2M's scale.
PARITY_ABS_TOL <- c(core = 1e-10, fitage = 1e-10, packs = 1e-6, horvath = 1e-10)

# a clock whose goldens came from the Horvath online calculator. Read off the
# declared oracle, so re-generating a fixture upstream retires the block on its
# own -- the way the DNAmFitAge skip should have and did not.
is_horvath_online <- function(id) {
  any(vapply(
    clock_fixtures(id) %||% list(),
    function(fx) identical(as.character(fx[["oracle"]]), "horvath_online"),
    logical(1)
  ))
}

# A `sample_scale` op z-scores each sample over EVERY probe in the input matrix,
# so the scoring panel is not a sufficient input: subsetting first moves each
# sample's mean/sd and the score with it (~1.8e1 off, 82% relative, measured).
# The oracle ran on the whole array, so the fixture must too. The predicate is
# the package's own -- calc_clocks() reads the same declared op to decide which
# matrix to hand that clock's branch, so the two cannot drift.
needs_full_panel <- clock_needs_full_panel

# which block a clock is graded in. Derived from the catalog, never a clock list.
parity_block <- function(id) {
  if (is_horvath_online(id)) {
    "horvath"
  } else if (clock_is_external(id)) {
    "packs"
  } else if (identical(clock_group_id(id), "DNAmFitAge")) {
    "fitage"
  } else {
    "core"
  }
}

# max, not median: a median passes while half the samples are arbitrarily wrong,
# which is the same blindness that got correlation retired as a gate.
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

skip_if_no_cohort <- function(cohort) {
  testthat::skip_if_not(
    parity_on,
    "parity tier off (set MC_PARITY=1, e.g. via dev test_parity())"
  )
  testthat::skip_if(
    is.null(cohort_cons[[cohort]]),
    paste0(cohort, " fixture not staged")
  )
}

# known gaps: skip so the suite stays green. Empty today -- Zhang2019 left on
# 2026-07-25 once the tier started feeding it the whole array (needs_full_panel).
KNOWN_PARITY_GAPS <- character(0)

# The Horvath online calculator filled every completely-absent probe server-side
# with a per-probe constant it does not publish, and we hold only the raw betas.
# Measured: the residual tracks the absent-probe count and nothing else -- every
# (clock, cohort) pair with zero absent probes already agrees to ~1e-8 relative
# (Hannum/DNAmCystatinC @ 450K, DNAmTL/DNAmTIMP1/Horvath2 @ EPICv1), which no
# wrong coefficient or wrong matmul survives. So the gap is in the oracle's
# input, not our arithmetic, and no tolerance can express that honestly.
# Horvath1 additionally needs BMIQ, which we deliberately vendor no gold
# standard for. See DECISIONS 2026-07-25.
HORVATH_ONLINE_GAP <- paste0(
  "horvath_online oracle -- server-side fill of absent probes is unpublished. ",
  "Pairs with no absent probes already match to ~1e-8"
)

# whole-group gaps, keyed by group id. Empty today -- DNAmFitAge is no longer
# skipped, it runs as its own block below. Kept as a separate map because group
# ids and clock ids share a namespace (DNAmFitAge is both), so one flat map
# would be ambiguous about which a key meant.
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

# (clock, cohort) pairs upstream declares a fixture for, for one block.
parity_targets <- function(block) {
  out <- list()
  for (id in names(mc_catalog)) {
    if (!identical(parity_block(id), block)) {
      next
    }
    for (fx in clock_fixtures(id) %||% list()) {
      out[[length(out) + 1L]] <- list(
        id = id,
        cohort = as.character(fx[["cohort"]])
      )
    }
  }
  out
}

run_parity_target <- function(clock_id, cohort) {
  skip_if_no_cohort(cohort)
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
  DNAm <- if (any(vapply(seq_ids, needs_full_panel, logical(1)))) {
    cohort_betas_full(cohort_cons[[cohort]])
  } else {
    cpgs <- panels_union(clock_panels(seq_ids, packs))
    cohort_betas(cohort_cons[[cohort]], cpgs)
  }
  # parity gates numbers, not coverage policy
  res <- calc_clocks(
    DNAm,
    request,
    pheno = cohort_pheno(cohort),
    from = packs,
    min_clocks_coverage = 0,
    min_samples_coverage = 0
  )
  # routed member scores land on the alias column for that sex's samples
  expect_parity(res$scores[, request], clock_id, cohort)
}

# The blocks below are GENERATED from clock_fixtures(), so a fixture upstream
# drops produces no test at all -- not a skip, not a failure, just two fewer
# passes in a green run. This census is the guard on that generator: it turns a
# missing declaration into a failure. It needs no duckdb and no staged cohort,
# only the committed catalog, so it is gated on the tier flag alone.
test_that("every clock declares a fixture for every registry cohort", {
  skip_if_not(parity_on, "parity tier off (set MC_PARITY=1)")

  ids <- names(mc_catalog)
  declared <- lapply(stats::setNames(ids, ids), function(id) {
    unique(vapply(
      clock_fixtures(id) %||% list(),
      function(fx) as.character(fx[["cohort"]]),
      character(1)
    ))
  })

  # sex-routed aliases carry no fixture of their own -- their members do.
  # Derived from the routing registry, so a new family needs no edit here.
  aliases <- unique(unlist(sex_routed_members()$alias))
  expect_setequal(names(Filter(function(x) !length(x), declared)), aliases)

  # every other clock: both cohorts. Partial coverage is a gap, not a pass.
  rest <- setdiff(ids, aliases)
  incomplete <- rest[!vapply(
    rest,
    function(id) setequal(declared[[id]], PARITY_COHORTS),
    logical(1)
  )]
  expect_equal(incomplete, character(0))

  # a cohort outside the registry would generate targets that can never run
  expect_setequal(unique(unlist(declared)), PARITY_COHORTS)
})

# block 1 -- bundled clocks outside the DNAmFitAge family
for (target in parity_targets("core")) {
  local({
    clock_id <- target$id
    cohort <- target$cohort
    test_that(paste0("parity: ", clock_id, " @ ", cohort), {
      run_parity_target(clock_id, cohort)
    })
  })
}

# block 2 -- the DNAmFitAge family, at the core tolerances. Green as of
# 2026-07-25: upstream regenerated these from author code, so the old group skip
# was hiding 28 passing tests.
for (target in parity_targets("fitage")) {
  local({
    clock_id <- target$id
    cohort <- target$cohort
    test_that(paste0("parity (fitage): ", clock_id, " @ ", cohort), {
      run_parity_target(clock_id, cohort)
    })
  })
}

# block 3 -- external packs (PCClocks, SystemsAge, PCBrainAge). Relaxed absolute
# bound only; the relative bound is the same 1e-10 the other blocks get.
for (target in parity_targets("packs")) {
  local({
    clock_id <- target$id
    cohort <- target$cohort
    test_that(paste0("parity (packs): ", clock_id, " @ ", cohort), {
      run_parity_target(clock_id, cohort)
    })
  })
}

# block 4 -- Horvath online calculator oracles. Skipped wholesale, at the core
# tolerances: the goldens were produced on betas the server filled and (for
# DNAmAge) BMIQ'd, and neither is published, so the gap is not ours to close.
# Kept as its own block so the skip is one visible statement rather than 30
# scattered ones -- and so it retires itself if the oracle ever changes.
for (target in parity_targets("horvath")) {
  local({
    clock_id <- target$id
    cohort <- target$cohort
    test_that(paste0("parity (horvath online): ", clock_id, " @ ", cohort), {
      run_parity_target(clock_id, cohort)
    })
  })
}

# both PhysAge composites in one call, per cohort
for (cohort_i in PARITY_COHORTS) {
  local({
    cohort <- cohort_i
    test_that(
      paste0("PhysAge composites match the author fixtures @ ", cohort),
      {
        skip_if_no_cohort(cohort)
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
