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

test_that("a resolved panel's positions index the usable set they came from", {
  skip_on_cran()
  # DunedinPACE brings a norm panel, so both roles are non-empty
  spec <- mc_spec(c("Hannum", "DunedinPACE"))
  DNAm <- random_betas(panels_union(spec$panels), n = 4L)
  facts <- mc_cohort(DNAm, spec)
  usable <- facts$usable_cols

  # one assertion per panel role over the whole request, not a loop per clock
  pull <- function(field) {
    unlist(
      lapply(facts$cpg_list$per_clock, function(x) x[[field]]),
      use.names = FALSE
    )
  }
  expect_equal(usable[pull("score_present_idx")], pull("score_present"))
  expect_equal(usable[pull("norm_present_idx")], pull("norm_present"))
  expect_true(length(pull("norm_present")) > 0L)
})

test_that("a position outside its axis is refused rather than indexed", {
  skip_on_cran()
  # 0 would silently drop an element and a high position yields NA
  block <- list(usable_idx = 1:5)
  expect_error(block_cols(0L, block))
  expect_error(block_cols(6L, block))

  fx <- chunk_run()
  # same CpGs, different order: every position stays in range
  desynced <- fx$run$facts
  desynced$usable_cols <- rev(desynced$usable_cols)
  expect_error(mc_block(fx$cohort$DNAm, fx$cohort$spec, desynced))
})

test_that("a panel's positions resolve to the panel's own CpGs in the block", {
  skip_on_cran()
  fx <- chunk_run()
  facts <- fx$run$facts
  block <- mc_block(fx$cohort$DNAm, fx$cohort$spec, facts)
  parts <- facts$cpg_list$panel_index$score$parts
  pull <- function(x, f) unlist(lapply(x, function(p) p[[f]]), use.names = FALSE)

  # the matmul reads coef by name off cols, so the columns must be those CpGs
  got <- lapply(parts, function(p) {
    colnames(observed_panel(p$present, p$present_idx, block)$values)
  })
  expect_equal(unlist(got, use.names = FALSE), pull(parts, "present"))

  # the cohort-mean overlay filters names and positions with one mask
  cached <- lapply(parts, function(p) {
    cached_cols(p$present, p$present_idx, block)
  })
  # the cohort carries partial NA, so this is not vacuous
  expect_true(length(pull(cached, "cols")) > 0L)
  expect_equal(facts$usable_cols[pull(cached, "idx")], pull(cached, "cols"))
})

test_that("an id join refuses a repeated key and an unmatched one", {
  skip_on_cran()
  # a repeated right key silently takes the first match
  expect_error(id_index("a", c("a", "a"), "T"))
  # an unmatched key silently yields an NA row
  expect_error(id_index(c("a", "b"), "a", "T"))
  # and the two policies that are not a defect
  expect_equal(id_index(c("b", "a", "c"), c("a", "b"), "T", "drop"), c(2L, 1L))
  expect_equal(
    id_index(c("b", "a", "c"), c("a", "b"), "T", "na"),
    c(2L, 1L, NA)
  )
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
