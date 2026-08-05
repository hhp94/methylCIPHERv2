# rbind.mc_result: gates, batch labels, assembly, opt-in re-finalize.

# disjoint sample blocks from one cohort, same CpG set
bind_blocks <- function(clocks = "Hannum", n = 12L, blocks = 3L) {
  DNAm <- random_betas(clock_cpgs(clocks), n = n)
  idx <- split(seq_len(n), rep(seq_len(blocks), length.out = n))
  list(
    DNAm = DNAm,
    parts = unname(lapply(idx, function(i) DNAm[i, , drop = FALSE]))
  )
}

bind_records <- function(clocks = "Hannum", n = 12L, blocks = 3L) {
  fx <- bind_blocks(clocks, n, blocks)
  fx$records <- lapply(fx$parts, function(m) calc_clocks(m, clocks))
  fx
}

# the two default shapes, built once each. tests that only read them share a draw.
BIND3 <- bind_records()
BIND2 <- bind_records(blocks = 2L)

# cross-sample clock is derived. this file is the main gate on `pending`.
XS_ID <- Filter(clock_is_cross_sample, resolve_clocks("all"))[[1L]]

# whole-cohort vs per-block records for the cross-sample clock
xs_split <- function(n = 12L, blocks = list(1:6, 7:12)) {
  DNAm <- random_betas(clock_cpgs(XS_ID), n = n)
  list(
    DNAm = DNAm,
    whole = calc_clocks(DNAm, XS_ID),
    parts = lapply(blocks, function(i) {
      calc_clocks(DNAm[i, , drop = FALSE], XS_ID)
    })
  )
}

test_that("disjoint records bind into one labelled union", {
  fx <- BIND3
  out <- do.call(rbind, fx$records)

  expect_s3_class(out, "mc_result")
  expect_setequal(rownames(out$scores), rownames(fx$DNAm))
  expect_equal(colnames(out$scores), colnames(fx$records[[1]]$scores))

  # a bind is not a re-run: every row is the score its own record carried
  first <- fx$records[[1]]$scores
  expect_equal(out$scores[rownames(first), , drop = FALSE], first)

  # one batch per record, per sample, and the two axes agree
  expect_equal(length(unique(out$provenance$mc_batch_id)), 3L)
  expect_equal(length(out$provenance$mc_batch_id), nrow(fx$DNAm))
  expect_setequal(
    names(out$coverage$per_clock),
    unique(out$provenance$mc_batch_id)
  )

  # label is derived from the ids. rbind mints nothing.
  expect_equal(
    names(out$coverage$per_clock),
    unlist(lapply(fx$records, function(r) names(r$coverage$per_clock)))
  )

  # pheno stays aligned to the sample axis, by id
  expect_equal(out$pheno$ID, out$provenance$sample_id)
})

test_that("rbind refuses what the caller chose differently", {
  fx <- BIND2

  # overlapping sample ids, whole or partial
  expect_error(rbind(fx$records[[1]], fx$records[[1]]))
  both <- calc_clocks(fx$DNAm, "Hannum")
  expect_error(rbind(both, fx$records[[1]]))

  # differing score columns
  cpgs <- union(clock_cpgs("Hannum"), clock_cpgs("PhenoAge"))
  DNAm <- random_betas(cpgs, n = 8L)
  expect_error(rbind(
    calc_clocks(DNAm[1:4, , drop = FALSE], "Hannum"),
    calc_clocks(DNAm[5:8, , drop = FALSE], "PhenoAge")
  ))

  # differing pheno_id
  bl <- bind_blocks(blocks = 2L)
  ph <- function(m, id) {
    out <- data.frame(x = rownames(m), stringsAsFactors = FALSE)
    names(out) <- id
    out
  }
  mk <- function(m, id) {
    calc_clocks(m, "Hannum", pheno = ph(m, id), pheno_id = id)
  }
  expect_error(rbind(mk(bl$parts[[1]], "ID"), mk(bl$parts[[2]], "sample")))

  # differing normalize -- state the decision directly (bmiq unfit on U(0,1))
  other <- fx$records[[2]]
  other$provenance$normalized <- "Hannum"
  expect_error(rbind(fx$records[[1]], other))
})

test_that("rbind is associative, and argument names are dropped not adopted", {
  fx <- bind_records(n = 16L, blocks = 4L)
  flat <- do.call(rbind, fx$records)

  expect_equal(Reduce(rbind, fx$records), flat)
  expect_equal(
    rbind(
      rbind(fx$records[[1]], fx$records[[2]]),
      rbind(fx$records[[3]], fx$records[[4]])
    ),
    flat
  )
  expect_equal(nrow(flat$scores), 16L)
  expect_equal(length(unique(flat$provenance$mc_batch_id)), 4L)

  # every record's samples still share the one label they arrived with
  for (rec in fx$records) {
    got <- flat$provenance$mc_batch_id[
      match(rownames(rec$scores), flat$provenance$sample_id)
    ]
    expect_equal(unique(got), unique(rec$provenance$mc_batch_id))
  }

  # the split() idiom the feature exists for: names ignored, never a label
  named <- rbind(early = fx$records[[1]], late = fx$records[[2]])
  expect_equal(named, rbind(fx$records[[1]], fx$records[[2]]))
})

test_that("a batch label is a function of the sample id set and nothing else", {
  skip_on_cran()
  DNAm <- random_betas(clock_cpgs("Hannum"), n = 6L)
  a <- calc_clocks(DNAm, "Hannum")

  # same ids scored twice -> same label, whatever the betas did in between
  expect_equal(
    calc_clocks(DNAm, "Hannum")$provenance$mc_batch_id,
    a$provenance$mc_batch_id
  )

  # the id *set*, not its order
  shuffled <- calc_clocks(DNAm[rev(seq_len(6L)), , drop = FALSE], "Hannum")
  expect_equal(
    unique(shuffled$provenance$mc_batch_id),
    unique(a$provenance$mc_batch_id)
  )

  # different ids -> different label. sim_DNAm(suffix =) is how a caller gets there.
  s1 <- sim_DNAm("Hannum", n = 4L, suffix = "T1")
  s2 <- sim_DNAm("Hannum", n = 4L, suffix = "T2")
  expect_equal(length(intersect(rownames(s1$DNAm), rownames(s2$DNAm))), 0L)
  out <- rbind(
    calc_clocks(s1$DNAm, "Hannum", pheno = s1$pheno),
    calc_clocks(s2$DNAm, "Hannum", pheno = s2$pheno)
  )
  expect_equal(length(names(out$coverage$per_clock)), 2L)
})

test_that("pheno never carries row names, whatever came in", {
  skip_on_cran()
  # negative .row_names_info means automatic/compact row names.
  is_auto <- function(df) .row_names_info(df) < 0

  n <- 6L
  DNAm <- random_betas(clock_cpgs("Hannum"), n = n)
  ph <- data.frame(ID = rownames(DNAm), stringsAsFactors = FALSE)
  rownames(ph) <- paste0("row_", seq_len(n))

  expect_true(is_auto(calc_clocks(DNAm, "Hannum")$pheno))
  # reordered, so the id-join subsets out of order
  rev_ph <- ph[rev(seq_len(n)), , drop = FALSE]
  expect_true(is_auto(calc_clocks(DNAm, "Hannum", pheno = rev_ph)$pheno))

  # and the bind keeps it, without a reset of its own
  half <- function(i) {
    calc_clocks(
      DNAm[i, , drop = FALSE],
      "Hannum",
      pheno = ph[i, , drop = FALSE]
    )
  }
  bound <- rbind(half(1:3), half(4:6))
  expect_true(is_auto(bound$pheno))
})

test_that("the coverage frames span the bind, keyed by batch", {
  fx <- BIND3
  out <- do.call(rbind, fx$records)

  cc <- clocks_coverage(out)
  one <- clocks_coverage(fx$records[[1]])
  expect_equal(names(cc)[[1]], "clock_id")
  expect_equal(nrow(cc), nrow(one) * 3L)
  expect_setequal(unique(cc$mc_batch_id), unique(out$provenance$mc_batch_id))

  # a bound row is the row its own record carried
  b1 <- unique(fx$records[[1]]$provenance$mc_batch_id)
  expect_equal(
    cc[cc$mc_batch_id == b1 & cc$clock_id == "Hannum", "score_needed"],
    one[one$clock_id == "Hannum", "score_needed"]
  )

  # every sample once, under its own batch
  sc <- samples_coverage(out)
  expect_setequal(sc$id, rownames(fx$DNAm))
  expect_equal(
    sc$mc_batch_id[match(out$provenance$sample_id, sc$id)],
    out$provenance$mc_batch_id
  )
})

test_that("the batch label reaches an exit frame only when there is more than one", {
  fx <- BIND3
  # all four exits must gain and lose the column together
  exits <- function(x) {
    ph <- mc_pheno(x$provenance$sample_id, Age = mc_ages(nrow(x$scores)))
    list(
      as.data.frame(x),
      as.data.frame(x, long = FALSE),
      calc_accel(x, data = ph, long = FALSE),
      clocks_coverage(x),
      samples_coverage(x)
    )
  }

  # one batch means one repeated hash, which tells the reader nothing
  for (df in exits(fx$records[[1]])) {
    expect_false("mc_batch_id" %in% names(df))
  }
  # and all of them carry it, last, once there is something to tell apart
  for (df in exits(do.call(rbind, fx$records))) {
    expect_equal(names(df)[[length(df)]], "mc_batch_id")
  }
})

test_that("a probe all-NA in one batch is recorded there, not merged away", {
  cpgs <- clock_cpgs("Hannum")
  DNAm <- random_betas(cpgs, n = 12L)
  hole <- cpgs[[1]]
  # all-NA within block 1, ordinary elsewhere -- the systematic per-batch offset
  DNAm[1:6, hole] <- NA_real_

  r1 <- calc_clocks(DNAm[1:6, , drop = FALSE], "Hannum")
  r2 <- calc_clocks(DNAm[7:12, , drop = FALSE], "Hannum")
  cov <- rbind(r1, r2)$coverage$per_clock
  b1 <- unique(r1$provenance$mc_batch_id)
  b2 <- unique(r2$provenance$mc_batch_id)

  expect_true(hole %in% cov[[b1]]$Hannum$missing_cpgs)
  expect_false(hole %in% cov[[b2]]$Hannum$missing_cpgs)
  expect_equal(
    cov[[b1]]$Hannum$score_present,
    cov[[b2]]$Hannum$score_present - 1L
  )
})

test_that("binding a cohort-reducing clock says so, and an ordinary one does not", {
  parts <- xs_split()$parts
  # the id is pinned because the note's content is which columns are unresolved
  expect_message(do.call(rbind, parts), XS_ID, fixed = TRUE)

  # one record in, one batch out -- its reduction still spans its whole cohort
  expect_no_message(rbind(parts[[1]]))

  # an ordinary clock retains no pending, so there is nothing to say
  fx <- BIND2
  expect_equal(fx$records[[1]]$provenance$pending, list())
  expect_no_message(do.call(rbind, fx$records))
})

test_that("refinalize composes: it reads pending and never consumes it", {
  id <- XS_ID
  fx <- xs_split(n = 18L, blocks = list(1:6, 7:12, 13:18))
  whole <- fx$whole
  r <- fx$parts
  by_row <- rownames(whole$scores)

  # bound, the column is three within-batch reductions -- not the cohort number
  bound <- suppressMessages(do.call(rbind, r))
  expect_setequal(names(bound$provenance$pending), id)
  expect_equal(nrow(bound$provenance$pending[[id]]), 18L)
  expect_false(isTRUE(all.equal(
    unname(bound$scores[by_row, id]),
    unname(whole$scores[, id])
  )))

  # re-finalize partway, bind more, re-finalize again
  suppressMessages({
    step <- refinalize_clocks(rbind(r[[1]], r[[2]]))
    out <- refinalize_clocks(rbind(r[[3]], step))
  })

  # pending accumulates across the interleaving rather than being spent
  expect_equal(nrow(step$provenance$pending[[id]]), 12L)
  expect_equal(nrow(out$provenance$pending[[id]]), 18L)
  expect_equal(out$scores[by_row, id], whole$scores[, id])

  # so a second call changes nothing, and the order of the folds does not either
  suppressMessages({
    expect_equal(refinalize_clocks(out)$scores, out$scores)
    flat <- refinalize_clocks(bound)
  })
  expect_equal(flat$scores[by_row, id], out$scores[by_row, id])
})

test_that("every finalizer resolves pending", {
  id <- XS_ID
  fx <- xs_split()
  whole <- fx$whole
  by_row <- rownames(whole$scores)

  suppressMessages({
    bound <- do.call(rbind, fx$parts)
    m <- as.matrix(bound)
    d <- as.data.frame(bound, long = FALSE)
    ph <- mc_pheno(bound$provenance$sample_id, Age = mc_ages(12L))
    acc <- calc_accel(bound, type = "diff", data = ph, long = FALSE)
  })

  # leaving the mc_result structure resolves pending, whichever exit is used
  expect_equal(unname(m[by_row, id]), unname(whole$scores[, id]))
  expect_equal(d[[id]][match(by_row, d[[1L]])], unname(whole$scores[, id]))
  # diff is score - Age, so it re-derives the same finalized column
  expect_equal(
    acc[[paste0(id, "_diff")]] + ph$Age[match(acc$ID, ph$ID)],
    unname(whole$scores[acc$ID, id])
  )
})
