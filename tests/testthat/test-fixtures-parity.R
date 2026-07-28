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

# packs: abs tol 1e-6 (scale ~3e6). see DECISIONS 2026-07-25
PARITY_ABS_TOL <- c(core = 1e-10, fitage = 1e-10, packs = 1e-6, horvath = 1e-10)

# horvath-online oracle clocks (declared, not listed)
is_horvath_online <- function(id) {
  any(vapply(
    clock_fixtures(id) %||% list(),
    function(fx) identical(as.character(fx[["oracle"]]), "horvath_online"),
    logical(1)
  ))
}

# sample_scale needs the full array (same predicate as calc_clocks)
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

# known gaps (clock- or clock@cohort-keyed). empty today.
KNOWN_PARITY_GAPS <- character(0)

# horvath block skipped: oracle filled absent probes server-side (DECISIONS 2026-07-25)
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
    ext_data = packs,
    min_clocks_coverage = 0,
    min_samples_coverage = 0
  )
  # routed member scores land on the alias column for that sex's samples
  expect_parity(res$scores[, request], clock_id, cohort)
}

# census: every catalog clock declares a fixture per cohort (guards the generator)
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

  # sex-routed aliases declare no fixture (members do)
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

# block 2 -- DNAmFitAge family at core tolerances
for (target in parity_targets("fitage")) {
  local({
    clock_id <- target$id
    cohort <- target$cohort
    test_that(paste0("parity (fitage): ", clock_id, " @ ", cohort), {
      run_parity_target(clock_id, cohort)
    })
  })
}

# block 3 -- external packs (relaxed abs tol only)
for (target in parity_targets("packs")) {
  local({
    clock_id <- target$id
    cohort <- target$cohort
    test_that(paste0("parity (packs): ", clock_id, " @ ", cohort), {
      run_parity_target(clock_id, cohort)
    })
  })
}

# block 4 -- horvath online oracles (skipped wholesale)
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
