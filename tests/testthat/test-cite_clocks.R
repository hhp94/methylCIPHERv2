# the whole-catalog citation frame, built once
ALL_CITES <- as.data.frame(cite_clocks("all"))

# the vendored .bib, read once (absent from a bare source checkout)
BIB_PATH <- system.file(
  "bibliography",
  "clocks.bib",
  package = "methylCIPHERv2"
)
BIB_LINES <- if (nzchar(BIB_PATH) && file.exists(BIB_PATH)) {
  readLines(BIB_PATH, warn = FALSE, encoding = "UTF-8")
}

test_that("cite_clocks speaks the same clock tokens as calc_clocks", {
  one <- cite_clocks("Hannum")
  expect_s3_class(one, "mc_citation")
  expect_equal(unique(as.data.frame(one)$clock_id), "Hannum")
  expect_true(any(grepl("^@", one$bibtex)))

  expect_setequal(
    unique(as.data.frame(cite_clocks("GrimAge"))$clock_id),
    resolve_clocks("GrimAge")
  )
  expect_setequal(unique(ALL_CITES$clock_id), resolve_clocks("all"))

  # one entry per distinct paper, however many clocks share it
  many <- cite_clocks(c("Hannum", "Horvath1", "PhenoAge"))
  expect_equal(
    sum(grepl("^@", many$bibtex)),
    length(unique(as.data.frame(many)$bib_key))
  )

  expect_error(cite_clocks(42)) # the default method refuses
  expect_error(cite_clocks("DNAmFitAge_Female")) # a routing target is not callable
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

test_that("every catalog clock reaches a citation, an alias through its donor", {
  skip_on_cran()
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

  # exactly one primary per cited clock
  expect_true(all(
    tapply(
      mc_citations[["role"]] == "primary",
      mc_citations[["clock_id"]],
      sum
    ) ==
      1L
  ))

  # and the alias's own links are the donor's
  alias <- aliases[[1L]]
  links <- as.data.frame(cite_clocks(alias))
  expect_equal(unique(links$clock_id), alias)
  expect_equal(
    links$bib_key,
    mc_citations$bib_key[
      mc_citations$clock_id == mc_catalog[[alias]][["donor_clock_id"]]
    ]
  )
})

test_that("the citation frame carries the paper's own fields", {
  skip_on_cran()
  paper_cols <- c("title", "author", "year", "journal", "doi", "url")
  expect_true(all(paper_cols %in% names(ALL_CITES)))
  # every .bib entry declares these, so an NA here means the join lost a key
  for (col in paper_cols) {
    expect_false(anyNA(ALL_CITES[[col]]))
  }

  skip_if(is.null(BIB_LINES))
  keys <- sub(
    "^@[^{]+\\{([^,]+),.*$",
    "\\1",
    grep("^@", BIB_LINES, value = TRUE)
  )
  expect_true(all(unique(mc_citations[["bib_key"]]) %in% trimws(keys)))
})
