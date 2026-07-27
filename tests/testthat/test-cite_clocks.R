test_that("every catalog clock reaches at least one citation", {
  ids <- mc_index[["clock_id"]]
  aliases <- unique(unname(sex_routed_members()$alias))

  # aliases are package-minted and cite through their donor
  expect_true(all(setdiff(ids, aliases) %in% mc_citations[["clock_id"]]))
  donors <- vapply(
    mc_catalog[aliases],
    function(e) as.character(e[["donor_clock_id"]]),
    character(1L)
  )
  expect_true(all(donors %in% mc_citations[["clock_id"]]))

  expect_true(all(mc_index[["n_citations"]] >= 1L))
  expect_equal(
    mc_index[["n_citations"]][match("Hannum", ids)],
    sum(mc_citations[["clock_id"]] == "Hannum")
  )
  # exactly one primary per cited clock
  expect_true(all(
    tapply(
      mc_citations[["role"]] == "primary",
      mc_citations[["clock_id"]],
      sum
    ) ==
      1L
  ))
})

test_that("cite_clocks speaks the same clock tokens as calc_clocks", {
  one <- cite_clocks("Hannum")
  expect_s3_class(one, "mc_citation")
  expect_equal(unique(as.data.frame(one)$clock_id), "Hannum")
  expect_true(any(grepl("^@", one$bibtex)))

  grp <- cite_clocks("GrimAge")
  expect_setequal(
    unique(as.data.frame(grp)$clock_id),
    resolve_clocks("GrimAge")
  )
  expect_setequal(
    unique(as.data.frame(cite_clocks("all"))$clock_id),
    resolve_clocks("all")
  )

  # one entry per distinct paper, however many clocks share it
  many <- cite_clocks(c("Hannum", "Horvath1", "PhenoAge"))
  expect_equal(
    sum(grepl("^@", many$bibtex)),
    length(unique(as.data.frame(many)$bib_key))
  )

  expect_error(cite_clocks(42))
  expect_error(cite_clocks("DNAmFitAge_Female"))
})

test_that("a sex-routed alias cites through its donor", {
  alias <- unique(unname(sex_routed_members()$alias))[[1L]]
  donor <- mc_catalog[[alias]][["donor_clock_id"]]
  links <- as.data.frame(cite_clocks(alias))

  expect_equal(unique(links$clock_id), alias)
  expect_equal(
    links$bib_key,
    mc_citations$bib_key[mc_citations$clock_id == donor]
  )
})

test_that("a result cites the clocks it reported", {
  sim <- sim_DNAm(c("Hannum", "PhenoAge"), n = 3L, Age = TRUE, Female = TRUE)
  res <- calc_clocks(sim$DNAm, c("Hannum", "PhenoAge"), pheno = sim$pheno)
  cit <- cite_clocks(res)

  expect_setequal(unique(as.data.frame(cit)$clock_id), colnames(res$scores))

  # toBibtex is the export path: every key survives a round-trip to disk
  path <- withr::local_tempfile(fileext = ".bib")
  writeLines(toBibtex(cit), path, useBytes = TRUE)
  written <- grep("^@", readLines(path), value = TRUE)
  expect_true(all(
    unique(as.data.frame(cit)$bib_key) %in%
      sub("^@[^{]+\\{([^,]+),.*$", "\\1", written)
  ))
})

test_that("every shipped bib_key resolves in clocks.bib", {
  bib <- system.file("bibliography", "clocks.bib", package = "methylCIPHERv2")
  skip_if(!nzchar(bib) || !file.exists(bib))
  keys <- sub(
    "^@[^{]+\\{([^,]+),.*$",
    "\\1",
    grep(
      "^@",
      readLines(bib, warn = FALSE),
      value = TRUE
    )
  )
  expect_true(all(unique(mc_citations[["bib_key"]]) %in% trimws(keys)))
})
