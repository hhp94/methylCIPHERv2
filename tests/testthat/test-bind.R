# rbind.mc_result: gates, batch labels, assembly, opt-in re-finalize.

# one cohort cut into disjoint sample blocks, each scored on its own. the
# blocks share a CpG set and differ only in which samples they hold, which is
# what a user projecting with clock_cpgs() and blocking produces
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

test_that("disjoint records bind into one labelled union", {
  fx <- bind_records()
  out <- do.call(rbind, fx$records)

  expect_s3_class(out, "mc_result")
  expect_setequal(rownames(out$scores), rownames(fx$DNAm))
  expect_equal(colnames(out$scores), colnames(fx$records[[1]]$scores))

  # a bind is not a re-run: every row is the score its own record carried
  first <- fx$records[[1]]$scores
  expect_equal(out$scores[rownames(first), , drop = FALSE], first)

  # one batch per record, per sample, and the two axes agree
  expect_equal(length(unique(out$provenance$batch)), 3L)
  expect_equal(length(out$provenance$batch), nrow(fx$DNAm))
  expect_setequal(names(out$coverage$per_clock), unique(out$provenance$batch))

  # the label is derived from the ids, so each record kept the one it was born
  # with -- rbind mints nothing
  expect_equal(
    names(out$coverage$per_clock),
    unlist(lapply(fx$records, function(r) names(r$coverage$per_clock)))
  )
})

test_that("a batch label is a function of the sample ids and nothing else", {
  DNAm <- random_betas(clock_cpgs("Hannum"), n = 6L)
  a <- calc_clocks(DNAm, "Hannum")
  b <- calc_clocks(DNAm, "Hannum")
  # same ids scored twice -> same label, whatever the betas did in between
  expect_equal(a$provenance$batch, b$provenance$batch)

  # the id *set*, not its order: re-scoring a reordered block keeps the label
  shuffled <- calc_clocks(DNAm[rev(seq_len(6L)), , drop = FALSE], "Hannum")
  expect_equal(unique(shuffled$provenance$batch), unique(a$provenance$batch))

  # different ids -> different label
  other <- DNAm
  rownames(other) <- paste0(rownames(other), "_T2")
  expect_false(
    identical(
      calc_clocks(other, "Hannum")$provenance$batch,
      a$provenance$batch
    )
  )
})

test_that("colliding sample ids throw", {
  fx <- bind_records(blocks = 2L)
  expect_error(rbind(fx$records[[1]], fx$records[[1]]))
  # partial overlap is the same fault
  both <- calc_clocks(fx$DNAm, "Hannum")
  expect_error(rbind(both, fx$records[[1]]))
})

test_that("differing score columns throw", {
  n <- 8L
  cpgs <- union(clock_cpgs("Hannum"), clock_cpgs("PhenoAge"))
  DNAm <- random_betas(cpgs, n = n)
  r1 <- calc_clocks(DNAm[1:4, , drop = FALSE], "Hannum")
  r2 <- calc_clocks(DNAm[5:8, , drop = FALSE], "PhenoAge")
  expect_error(rbind(r1, r2))
})

test_that("differing pheno_id throws", {
  fx <- bind_blocks(blocks = 2L)
  ph <- function(m, id) {
    out <- data.frame(x = rownames(m), stringsAsFactors = FALSE)
    names(out) <- id
    out
  }
  r1 <- calc_clocks(
    fx$parts[[1]],
    "Hannum",
    pheno = ph(fx$parts[[1]], "ID"),
    pheno_id = "ID"
  )
  r2 <- calc_clocks(
    fx$parts[[2]],
    "Hannum",
    pheno = ph(fx$parts[[2]], "sample"),
    pheno_id = "sample"
  )
  expect_error(rbind(r1, r2))
})

test_that("pheno never carries row names, whatever came in", {
  # negative .row_names_info means automatic/compact -- no second identity
  # alongside the id column
  is_auto <- function(df) .row_names_info(df) < 0

  n <- 6L
  DNAm <- random_betas(clock_cpgs("Hannum"), n = n)
  expect_true(is_auto(calc_clocks(DNAm, "Hannum")$pheno))

  ph <- data.frame(ID = rownames(DNAm), stringsAsFactors = FALSE)
  rownames(ph) <- paste0("row_", seq_len(n))
  expect_true(is_auto(calc_clocks(DNAm, "Hannum", pheno = ph)$pheno))

  # reordered, so the id-join subsets out of order
  rev_ph <- ph[rev(seq_len(n)), , drop = FALSE]
  expect_true(is_auto(calc_clocks(DNAm, "Hannum", pheno = rev_ph)$pheno))

  # and the bind keeps it, without a reset of its own
  bound <- rbind(
    calc_clocks(
      DNAm[1:3, , drop = FALSE],
      "Hannum",
      pheno = ph[1:3, , drop = FALSE]
    ),
    calc_clocks(
      DNAm[4:6, , drop = FALSE],
      "Hannum",
      pheno = ph[4:6, , drop = FALSE]
    )
  )
  expect_true(is_auto(bound$pheno))
  expect_equal(bound$pheno$ID, bound$provenance$sample_id)
})

test_that("passing a pheno and omitting it give the same carried columns", {
  fx <- bind_blocks(blocks = 2L)
  r1 <- calc_clocks(
    fx$parts[[1]],
    "Hannum",
    pheno = data.frame(ID = rownames(fx$parts[[1]]), stringsAsFactors = FALSE)
  )
  r2 <- calc_clocks(fx$parts[[2]], "Hannum")

  # an omitted pheno is materialized to the id column, so the two records are
  # the same shape and bind -- there is nothing left to refuse
  expect_equal(names(r1$pheno), "ID")
  expect_equal(names(r2$pheno), "ID")
  out <- rbind(r1, r2)
  expect_equal(nrow(out$pheno), nrow(out$scores))
  expect_equal(out$pheno$ID, out$provenance$sample_id)
})

test_that("records normalized differently do not bind", {
  fx <- bind_records(blocks = 2L)
  # the gate reads the recorded decision, so state it directly rather than
  # standing up a BMIQ fit that U(0,1) betas cannot support
  fx$records[[2]]$provenance$normalized <- "Hannum"
  expect_error(rbind(fx$records[[1]], fx$records[[2]]))
})

test_that("flat rbind, Reduce and do.call agree", {
  fx <- bind_records()
  flat <- rbind(fx$records[[1]], fx$records[[2]], fx$records[[3]])
  expect_equal(Reduce(rbind, fx$records), flat)
  expect_equal(do.call(rbind, fx$records), flat)
})

test_that("a named list binds -- split() names are not a labelling attempt", {
  n <- 12L
  DNAm <- random_betas(clock_cpgs("Hannum"), n = n)
  # the canonical blocking idiom: split() names its result by factor level
  blocks <- split(seq_len(n), rep(1:3, each = 4))
  recs <- lapply(blocks, function(i) {
    calc_clocks(DNAm[i, , drop = FALSE], "Hannum")
  })
  expect_true(!is.null(names(recs)))

  out <- do.call(rbind, recs)
  expect_equal(nrow(out$scores), n)
  # the names are dropped, not adopted -- labels stay derived
  expect_equal(length(unique(out$provenance$batch)), 3L)
  expect_false(any(c("1", "2", "3") %in% out$provenance$batch))

  # naming by hand is the same: ignored, never a label
  fx <- bind_records(blocks = 2L)
  named <- rbind(early = fx$records[[1]], late = fx$records[[2]])
  expect_equal(named, do.call(rbind, fx$records))
})

test_that("sim_DNAm suffixes ids so two blocks are disjoint, and they bind", {
  a <- sim_DNAm("Hannum", n = 4L, suffix = "T1")
  b <- sim_DNAm("Hannum", n = 4L, suffix = "T2")

  expect_equal(a$suffix, "T1")
  expect_true(all(endsWith(rownames(a$DNAm), "_T1")))
  # the pheno id column is the same id source, not a parallel one
  expect_equal(a$pheno$ID, rownames(a$DNAm))
  expect_equal(length(intersect(rownames(a$DNAm), rownames(b$DNAm))), 0L)

  out <- rbind(
    calc_clocks(a$DNAm, "Hannum", pheno = a$pheno),
    calc_clocks(b$DNAm, "Hannum", pheno = b$pheno)
  )
  expect_equal(length(names(out$coverage$per_clock)), 2L)
  expect_equal(nrow(out$scores), 8L)

  # unsuffixed is the default, and leaves the ids alone
  plain <- sim_DNAm("Hannum", n = 4L)
  expect_null(plain$suffix)
  expect_equal(rownames(plain$DNAm), paste0("sample", 1:4))
})

test_that("re-association is exact: labels derive from ids, so nothing moves", {
  fx <- bind_records(n = 16L, blocks = 4L)
  flat <- do.call(rbind, fx$records)
  nested <- rbind(
    rbind(fx$records[[1]], fx$records[[2]]),
    rbind(fx$records[[3]], fx$records[[4]])
  )

  # the whole record agrees, not just the labels -- no renumbering to undo
  expect_equal(nested, flat)
  expect_equal(nrow(nested$scores), 16L)
  expect_equal(length(unique(nested$provenance$batch)), 4L)

  # every record's samples still share one batch
  for (rec in fx$records) {
    got <- nested$provenance$batch[
      match(rownames(rec$scores), nested$provenance$sample_id)
    ]
    expect_equal(length(unique(got)), 1L)
    expect_equal(unique(got), unique(rec$provenance$batch))
  }
})

test_that("clocks_coverage is one row per (clock, batch)", {
  fx <- bind_records()
  out <- do.call(rbind, fx$records)

  cc <- clocks_coverage(out)
  one <- clocks_coverage(fx$records[[1]])
  # the hash is the key this frame is on, but it reads as noise -- so it sits
  # at the end, not in front of the clock id
  expect_equal(names(cc)[[length(cc)]], "batch")
  expect_equal(names(cc)[[1]], "clock_id")
  expect_equal(nrow(cc), nrow(one) * 3L)
  expect_setequal(unique(cc$batch), unique(out$provenance$batch))

  b1 <- unique(fx$records[[1]]$provenance$batch)
  expect_equal(
    cc[cc$batch == b1 & cc$clock_id == "Hannum", "score_needed"],
    one[one$clock_id == "Hannum", "score_needed"]
  )
})

test_that("samples_coverage carries every sample once, under its own batch", {
  fx <- bind_records()
  out <- do.call(rbind, fx$records)
  sc <- samples_coverage(out)

  expect_setequal(sc$id, rownames(fx$DNAm))
  expect_equal(
    sc$batch[match(out$provenance$sample_id, sc$id)],
    out$provenance$batch
  )
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
  b1 <- unique(r1$provenance$batch)
  b2 <- unique(r2$provenance$batch)

  expect_true(hole %in% cov[[b1]]$Hannum$missing_cpgs)
  expect_false(hole %in% cov[[b2]]$Hannum$missing_cpgs)
  expect_equal(
    cov[[b1]]$Hannum$score_present,
    cov[[b2]]$Hannum$score_present - 1L
  )
})

test_that("an ordinary clock retains no pending, and re-finalize says so", {
  fx <- bind_records(blocks = 2L)
  expect_equal(fx$records[[1]]$provenance$pending, list())

  out <- do.call(rbind, fx$records)
  expect_equal(out$provenance$pending, list())
  expect_message(same <- refinalize_clocks(out))
  expect_equal(same$scores, out$scores)
})

test_that("re-finalizing a bound record reproduces the single-pass column", {
  id <- "DNAmPhysAge"
  DNAm <- random_betas(clock_cpgs(id), n = 12L)
  whole <- calc_clocks(DNAm, id)
  parts <- lapply(list(1:6, 7:12), function(i) {
    calc_clocks(DNAm[i, , drop = FALSE], id)
  })
  bound <- do.call(rbind, parts)

  # the per-sample intermediates survive the bind, stacked over every sample
  expect_setequal(names(bound$provenance$pending), id)
  expect_equal(nrow(bound$provenance$pending[[id]]), 12L)

  # bound, the column is two within-batch z-scores -- not the cohort number
  by_row <- rownames(whole$scores)
  expect_false(isTRUE(all.equal(
    unname(bound$scores[by_row, id]),
    unname(whole$scores[, id])
  )))

  suppressMessages(fixed <- refinalize_clocks(bound))
  expect_equal(fixed$scores[by_row, id], whole$scores[, id])
})

test_that("binding a cohort-reducing clock says so, and an ordinary one does not", {
  id <- "DNAmPhysAge"
  DNAm <- random_betas(clock_cpgs(id), n = 12L)
  parts <- lapply(list(1:6, 7:12), function(i) {
    calc_clocks(DNAm[i, , drop = FALSE], id)
  })
  # the clock is named from the record's own pending, not from a list here
  expect_message(do.call(rbind, parts), id, fixed = TRUE)

  # one record in, one batch out -- its reduction still spans its whole cohort
  expect_no_message(rbind(parts[[1]]))

  fx <- bind_records(blocks = 2L)
  expect_no_message(do.call(rbind, fx$records))
})

test_that("refinalize composes: it reads pending and never consumes it", {
  id <- "DNAmPhysAge"
  DNAm <- random_betas(clock_cpgs(id), n = 18L)
  whole <- calc_clocks(DNAm, id)
  r <- lapply(list(1:6, 7:12, 13:18), function(i) {
    calc_clocks(DNAm[i, , drop = FALSE], id)
  })

  # re-finalize partway, bind more, re-finalize again
  suppressMessages({
    step <- refinalize_clocks(rbind(r[[1]], r[[2]]))
    grown <- rbind(r[[3]], step)
    out <- refinalize_clocks(grown)
  })

  # pending accumulates across the interleaving rather than being spent
  expect_equal(nrow(step$provenance$pending[[id]]), 12L)
  expect_equal(nrow(out$provenance$pending[[id]]), 18L)

  by_row <- rownames(whole$scores)
  expect_equal(out$scores[by_row, id], whole$scores[, id])

  # so a second call changes nothing, and the order of the folds does not either
  suppressMessages({
    expect_equal(refinalize_clocks(out)$scores, out$scores)
    flat <- refinalize_clocks(do.call(rbind, r))
  })
  expect_equal(flat$scores[by_row, id], out$scores[by_row, id])
})
