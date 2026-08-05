# chunk-safety: row subsets scored off fixed cohort facts match whole-cohort rows

# a cohort carrying every missingness shape at once
chunk_cohort <- function(ids, n = 30L) {
  spec <- mc_spec(ids)
  panel <- spec$needed_union
  DNAm <- random_betas(panel, n = n)

  # 1. ordinary partial NA, spread across blocks
  DNAm[1:5, panel[1]] <- NA
  # 2. partial cohort-wide but all-NA inside the first block
  DNAm[seq_len(n %/% 3L), panel[2]] <- NA
  # 3. fully absent everywhere -> vendor ref or dropped by policy
  DNAm[, panel[3]] <- NA

  list(spec = spec, DNAm = DNAm, panel = panel)
}

# stitch per-block fragments by clock id, restore cohort row order
stitch_fragments <- function(fragments, sample_id) {
  fragments <- lapply(fragments, function(f) Filter(Negate(is.null), f))
  ids <- unique(unlist(lapply(fragments, names), use.names = FALSE))
  out <- lapply(ids, function(p) {
    m <- do.call(rbind, lapply(fragments, function(f) f[[p]]))
    m[match(sample_id, rownames(m)), , drop = FALSE]
  })
  names(out) <- ids
  out
}

# whole vs three blocks off the same facts -- both finalize after assembly
split_score <- function(spec, DNAm, blocks) {
  facts <- mc_cohort(DNAm, spec, pheno = NULL)
  whole <- score_cohort(DNAm, spec, facts)
  parts <- lapply(blocks, function(i) {
    score_cohort(DNAm[i, , drop = FALSE], spec, facts)
  })
  fragments <- function(field) {
    stitch_fragments(lapply(parts, function(p) p[[field]]), facts$sample_id)
  }
  list(
    facts = facts,
    whole = whole,
    parts = parts,
    whole_scores = finalize_cross_sample(whole$scores, whole$pending),
    chunked_scores = finalize_cross_sample(
      fragments("scores"),
      fragments("pending")
    )
  )
}

# one clock per chunk-sensitive shape: linear, cohort reduction, and row-moment banking.
MIXED <- c("Hannum", "DNAmPhysAge", "Zhang2019EN")

# one cohort, scored whole and in three blocks. built on first use.
chunk_cache <- new.env(parent = emptyenv())
chunk_run <- function() {
  if (is.null(chunk_cache[["fx"]])) {
    # mc_spec() notes the full-matrix clock in the request. not the point here.
    suppressMessages({
      cohort <- chunk_cohort(MIXED)
      chunk_cache[["fx"]] <- list(
        cohort = cohort,
        run = split_score(cohort$spec, cohort$DNAm, list(1:10, 11:20, 21:30))
      )
    })
  }
  chunk_cache[["fx"]]
}

test_that("scoring a row subset equals scoring the whole cohort", {
  skip_on_cran()
  fx <- chunk_run()
  cohort <- fx$cohort
  run <- fx$run

  # both chunk-sensitive shapes are actually in the request
  expect_true(length(cohort$spec$cross_sample) > 0)
  expect_true(length(cohort$spec$moment_domains) > 0)

  # every clock: block-scored cohort finalizes to the single-pass answer
  for (id in cohort$spec$sequence) {
    expect_equal(run$chunked_scores[[id]], run$whole_scores[[id]])
  }
})

test_that("the scoring loop defers exactly the declared cohort-reducing set", {
  skip_on_cran()
  fx <- chunk_run()
  cohort <- fx$cohort
  run <- fx$run

  # deferred clocks leave the loop with an intermediate and no score ...
  expect_equal(sort(names(run$whole$pending)), sort(cohort$spec$cross_sample))
  for (id in cohort$spec$cross_sample) {
    expect_null(run$whole$scores[[id]])
    # ... which is per-sample and wider than one column
    expect_equal(nrow(run$whole$pending[[id]]), nrow(cohort$DNAm))
    expect_true(ncol(run$whole$pending[[id]]) > 1L)
  }

  # and the finalize turns every one of them into an ordinary score column
  for (id in cohort$spec$sequence) {
    expect_equal(dim(run$whole_scores[[id]]), c(nrow(cohort$DNAm), 1L))
  }
})

test_that("the sample-axis split comes off the catalog, not a clock list", {
  skip_on_cran()
  # only a declared cohort-reducing recipe is cross_sample. split is total over the pool.
  all_ids <- resolve_clocks("all")
  split <- split_cross_sample(all_ids)
  expect_equal(sort(c(split$per_sample, split$cross_sample)), sort(all_ids))
  expect_false(any(is.na(vapply(split$cross_sample, clock_cross_sample_at, 1))))
  expect_true(all(is.na(vapply(split$per_sample, clock_cross_sample_at, 1))))

  # sex-routed alias inherits the axis from its members
  aliases <- unique(unname(sex_routed_members()$alias))
  inherited <- vapply(
    aliases,
    function(a) {
      members <- unlist(clock_routing(a), use.names = FALSE)
      isTRUE(clock_is_cross_sample(a)) ==
        any(vapply(members, clock_is_cross_sample, logical(1)))
    },
    logical(1)
  )
  expect_true(all(inherited))
})

test_that("cohort facts classify a column the first block cannot see", {
  skip_on_cran()
  fx <- chunk_run()
  cohort <- fx$cohort
  facts <- fx$run$facts

  # partial cohort-wide, so it is a fill column even though block 1 is all NA
  expect_true(cohort$panel[2] %in% names(facts$partial_fill))
  # fully absent, so it is not usable and takes no fill value
  expect_false(cohort$panel[3] %in% facts$usable_cols)
  expect_false(cohort$panel[3] %in% names(facts$partial_fill))
})

test_that("coverage assembles from blocks by concatenate, never by sum", {
  skip_on_cran()
  fx <- chunk_run()
  cohort <- fx$cohort
  run <- fx$run

  for (id in cohort$spec$sequence) {
    whole <- run$whole$coverage$per_clock[[id]]

    # panel-derived fields are CpG counts. every block agrees with the single pass.
    for (p in run$parts) {
      expect_equal(
        p$coverage$per_clock[[id]]$score_imputed_partial,
        whole$score_imputed_partial
      )
    }
    expect_equal(
      run$parts[[1]]$coverage$per_clock[[id]]$score_present,
      whole$score_present
    )

    # per-sample miss concatenates over disjoint rows
    miss_whole <- run$whole$coverage$sample_miss$score[[id]]
    miss_parts <- unlist(lapply(
      run$parts,
      function(p) p$coverage$sample_miss$score[[id]]
    ))
    expect_equal(miss_parts[names(miss_whole)], miss_whole)
  }
})

test_that("the clocks gate throws out of mc_cohort, before anything is scored", {
  skip_on_cran()
  spec <- mc_spec("Hannum")
  panel <- spec$needed_union
  DNAm <- random_betas(panel, n = 6L)
  # strip most of the panel -- coverage falls under the default floor
  DNAm <- DNAm[, seq_len(length(panel) %/% 4L), drop = FALSE]

  expect_error(mc_cohort(DNAm, spec, pheno = NULL))
})
