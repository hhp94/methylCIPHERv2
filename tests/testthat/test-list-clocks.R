# rank is pinned on one typo only
test_that("did_you_mean puts the intended clock first, per namespace", {
  expect_equal(
    did_you_mean("systemage", suggestion_pools()$clock)[[1L]],
    "SystemsAge"
  )
  expect_equal(
    did_you_mean("grimage", suggestion_pools()$group)[[1L]],
    "GrimAge"
  )
})

test_that("clock and group suggestions stay in their own namespace", {
  pools <- suggestion_pools()
  expect_false("all" %in% c(pools$clock, pools$group, names(pools$clock)))
  expect_setequal(unname(pools$group), unique(mc_index[["group_id"]]))
  expect_true(all(did_you_mean("systemage", pools$group) %in% pools$group))
  expect_equal(did_you_mean("fitage", pools$group)[[1L]], "DNAmFitAge")
  expect_length(suggestion_bullets(c("zzz", "systemage")), 6L)
})

test_that("case is a suggestion, never a resolution", {
  expect_equal(
    did_you_mean("SYSTEMSAGE", suggestion_pools()$clock)[[1L]],
    "SystemsAge"
  )
  expect_error(resolve_clocks("SYSTEMSAGE"))
  # exact id resolves exactly, and a group token beats a same-named clock
  expect_equal(resolve_clocks("non_prcPhenoAge"), "non_prcPhenoAge")
  expect_length(resolve_clocks("prcPhenoAge"), 2L)
})

test_that("a mistyped routed member is pointed at its alias", {
  routed <- sex_routed_members()
  member <- names(routed$alias)[[1L]]
  expect_equal(
    did_you_mean(
      substr(member, 1L, nchar(member) - 1L),
      suggestion_pools()$clock
    )[[1L]],
    unname(routed$alias[[member]])
  )
})

test_that("suggestions never name a clock the user cannot request", {
  pools <- suggestion_pools()
  routed <- names(sex_routed_members()$alias)
  hits <- unlist(lapply(
    c("systemage", "fitage", "DNAmFitAge_Femal", "zzz"),
    did_you_mean,
    pool = pools$clock
  ))
  expect_false(any(hits %in% routed))
})

test_that("every listed callable clock actually resolves", {
  lc <- list_clocks(all_columns = TRUE)
  expect_equal(nrow(lc), length(mc_index[["clock_id"]]))
  expect_setequal(lc[["clock_id"]][lc[["callable"]]], resolve_clocks("all"))
  routed <- lc[!lc[["callable"]], ]
  for (i in seq_len(nrow(routed))) {
    seq_ids <- resolve_clocks_sequence(
      resolve_clocks(routed[["request_as"]][[i]])
    )
    expect_true(routed[["clock_id"]][[i]] %in% seq_ids)
  }
})

test_that("group_size is what the group token expands to", {
  lc <- list_clocks(all_columns = TRUE)
  for (g in unique(lc[["group_id"]])) {
    expect_equal(
      length(resolve_clocks(g)),
      lc[["group_size"]][match(g, lc[["group_id"]])]
    )
  }
})

test_that("every keyword resolves to a non-empty set of callable clocks", {
  callable <- resolve_clocks("all")
  for (tag in names(MC_TAGS)) {
    ids <- resolve_clocks(tag)
    expect_gt(length(ids), 0L)
    expect_true(all(ids %in% callable))
  }
})

test_that("keyword names cannot collide with a clock or group id", {
  expect_false(any(names(MC_TAGS) %in% mc_index[["clock_id"]]))
  expect_false(any(names(MC_TAGS) %in% mc_index[["group_id"]]))
  expect_false("all" %in% names(MC_TAGS))
  expect_equal(names(MC_TAGS), tolower(names(MC_TAGS)))
})

test_that("a keyword naming a dropped token is a hard stop, not a short set", {
  local_mocked_bindings(MC_TAGS = list(broken = c("Hannum", "NoSuchGroup")))
  expect_error(resolve_clocks("broken"))
})

test_that("keyword membership matches the tags column", {
  lc <- list_clocks()
  for (tag in names(MC_TAGS)) {
    tagged <- lc[["clock_id"]][grepl(tag, lc[["tags"]], fixed = TRUE)]
    expect_setequal(
      lc[["request_as"]][match(tagged, lc[["clock_id"]])],
      resolve_clocks(tag)
    )
    expect_setequal(list_clocks(tag = tag)[["clock_id"]], tagged)
  }
})

test_that("the default menu is narrow, and all_columns adds the rest", {
  narrow <- list_clocks()
  wide <- list_clocks(all_columns = TRUE)

  expect_equal(nrow(narrow), nrow(wide))
  expect_true(all(names(narrow) %in% names(wide)))
  expect_setequal(
    setdiff(names(wide), names(narrow)),
    c("callable", "group_size", "batch_dependent", "normalize")
  )
  # a filter does not change which columns come back
  expect_equal(names(list_clocks(group = "Dunedin")), names(narrow))
})

test_that("callable is dropped because request_as already carries it", {
  lc <- list_clocks(all_columns = TRUE)
  expect_equal(lc[["callable"]], lc[["request_as"]] == lc[["clock_id"]])
})

test_that("list_clock_tags returns the registry as a value", {
  out <- withVisible(list_clock_tags())
  expect_true(out$visible)
  expect_equal(out$value, MC_TAGS)
})

test_that("list_clocks filters, and rejects an unknown group", {
  expect_true(all(list_clocks(group = "GrimAge")[["group_id"]] == "GrimAge"))
  expect_gt(nrow(list_clocks(pattern = "horvath")), 0L)
  expect_equal(nrow(list_clocks(pattern = "nothing-matches-this")), 0L)
  # names a real group, so a broken suggestion pool is not mistaken for a reject
  expect_error(list_clocks(group = "Horvat"), "Horvath")
})
