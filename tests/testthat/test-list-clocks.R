test_that("list_clocks filters, rejects an unknown group, and widens on demand", {
  narrow <- list_clocks()
  wide <- list_clocks(all_columns = TRUE)

  expect_equal(nrow(narrow), nrow(wide))
  expect_setequal(
    setdiff(names(wide), names(narrow)),
    c("callable", "group_size", "batch_dependent", "normalize")
  )

  # a filter narrows rows, never columns
  expect_true(all(list_clocks(group = "GrimAge")[["group_id"]] == "GrimAge"))
  expect_equal(names(list_clocks(group = "Dunedin")), names(narrow))
  expect_gt(nrow(list_clocks(pattern = "horvath")), 0L)
  expect_equal(nrow(list_clocks(pattern = "nothing-matches-this")), 0L)

  # the tag filter selects what the same token resolves to
  tag <- names(MC_TAGS)[[1L]]
  expect_setequal(list_clocks(tag = tag)[["request_as"]], resolve_clocks(tag))

  # regex kept: names a real group so a broken suggestion is not a plain reject.
  expect_error(list_clocks(group = "Horvat"), "Horvath")
})

test_that("list_clock_tags returns the registry as a value", {
  out <- withVisible(list_clock_tags())
  expect_true(out$visible)
  expect_equal(out$value, MC_TAGS)
})

test_that("a token resolves exactly -- never by case or abbreviation", {
  skip_on_cran()
  expect_error(resolve_clocks("SYSTEMSAGE"))
  expect_error(resolve_clocks("Sys"))

  # an exact id resolves exactly, and a group token beats a same-named clock
  expect_equal(resolve_clocks("non_prcPhenoAge"), "non_prcPhenoAge")
  expect_setequal(
    resolve_clocks("prcPhenoAge"),
    mc_index[["clock_id"]][mc_index[["group_id"]] == "prcPhenoAge"]
  )

  # a routing target is refused by name
  expect_error(resolve_clocks(names(sex_routed_members()$alias)[[1L]]))
})

test_that("the callable pool holds exactly the clocks a user can request", {
  skip_on_cran()
  lc <- list_clocks(all_columns = TRUE)
  expect_equal(nrow(lc), length(mc_index[["clock_id"]]))
  expect_setequal(lc[["clock_id"]][lc[["callable"]]], resolve_clocks("all"))

  # a routed member is reachable only as its alias's dependency
  routed <- lc[!lc[["callable"]], ]
  reach <- vapply(
    seq_len(nrow(routed)),
    function(i) {
      seq_ids <- resolve_clocks_sequence(
        resolve_clocks(routed[["request_as"]][[i]])
      )
      routed[["clock_id"]][[i]] %in% seq_ids
    },
    logical(1L)
  )
  expect_true(all(reach))

  # group_size is what the group token expands to
  expect_equal(
    length(resolve_clocks("GrimAge")),
    lc[["group_size"]][match("GrimAge", lc[["group_id"]])]
  )
})

test_that("every keyword resolves to a non-empty set of callable clocks", {
  skip_on_cran()
  callable <- resolve_clocks("all")
  for (tag in names(MC_TAGS)) {
    ids <- resolve_clocks(tag)
    expect_gt(length(ids), 0L)
    expect_true(all(ids %in% callable))
  }

  # keyword names live in their own namespace
  expect_false(any(
    names(MC_TAGS) %in%
      c(mc_index[["clock_id"]], mc_index[["group_id"]], "all")
  ))

  # a keyword naming a dropped token is a hard stop, not a short set
  local_mocked_bindings(MC_TAGS = list(broken = c("Hannum", "NoSuchGroup")))
  expect_error(resolve_clocks("broken"))
})

test_that("a suggestion points at the intended clock, never at a routed member", {
  skip_on_cran()
  pools <- suggestion_pools()
  routed <- sex_routed_members()$alias

  expect_equal(did_you_mean("systemage", pools$clock)[[1L]], "SystemsAge")
  expect_equal(did_you_mean("fitage", pools$group)[[1L]], "DNAmFitAge")
  expect_false("all" %in% c(pools$clock, pools$group))

  # a mistyped routed member is pointed at its alias
  member <- names(routed)[[1L]]
  expect_equal(
    did_you_mean(
      substr(member, 1L, nchar(member) - 1L),
      pools$clock
    )[[1L]],
    unname(routed[[member]])
  )
  hits <- unlist(lapply(
    c("systemage", "fitage", "DNAmFitAge_Femal", "zzz"),
    did_you_mean,
    pool = pools$clock
  ))
  expect_false(any(hits %in% names(routed)))

  # both bad tokens are reported, and the good suggestion reaches the text
  bullets <- suggestion_bullets(c("zzz", "systemage"))
  expect_true(any(grepl("zzz", bullets, fixed = TRUE)))
  expect_true(any(grepl("SystemsAge", bullets, fixed = TRUE)))
})
