# external clock-data asset tests (file://, live network opt-in)

# digest is a Suggests-only dep (content-addresses fake pack filenames)
skip_if_not_installed("digest")

# fake external pack on disk, returns its provenance row
fake_asset <- function(dir, group = "FakeGroup", payload = NULL) {
  if (is.null(payload)) {
    payload <- list(
      group_id = group,
      encoding = "canonical_matrices",
      encoding_version = 3L,
      cpgs = c("cg00000029", "cg00000108"),
      coefficient_matrix = matrix(1:4, 2)
    )
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  phash <- digest::digest(payload, algo = "sha256")
  file <- sprintf("%s-%s.qs2", tolower(group), phash)
  rtag <- sub("\\.qs2$", "", file) # tag = filename stem (bare hex tags rejected by github)
  src <- file.path(dir, file)
  qs2::qs_save(payload, src)
  list(
    group_id = group,
    release_tag = rtag,
    file = file,
    size_bytes = as.numeric(file.size(src)),
    encoding = payload$encoding,
    encoding_version = payload$encoding_version,
    n_clocks = 1L,
    n_cpgs = length(payload$cpgs),
    .payload = payload,
    .src = as.character(src)
  )
}

# mock provenance registry + file:// download URLs
local_fake_registry <- function(rows, .env = parent.frame()) {
  if (!is.null(rows$group_id)) {
    rows <- stats::setNames(list(rows), rows$group_id)
  }
  testthat::local_mocked_bindings(
    mc_external_groups = function() names(rows),
    mc_asset = function(group_id) {
      row <- rows[[group_id]]
      if (is.null(row)) {
        stop("Not an external clock group: ", group_id, call. = FALSE)
      }
      row
    },

    mc_asset_url = function(row) {
      paste0(
        "file:///",
        normalizePath(row$.src, winslash = "/", mustWork = FALSE)
      )
    },
    .env = .env
  )
  invisible(rows)
}

test_that("the assets dir resolves arg > option > env > default", {
  withr::local_options(mc.assets_dir = NULL)
  withr::local_envvar(MC_ASSETS_DIR = NA)
  expect_equal(get_mc_assets_dir(), mc_default_assets_dir())

  withr::local_envvar(MC_ASSETS_DIR = "from-env")
  expect_equal(get_mc_assets_dir(), path.expand("from-env"))

  withr::local_options(mc.assets_dir = "from-option")
  expect_equal(get_mc_assets_dir(), path.expand("from-option"))

  # an explicit source beats both layers
  expect_equal(mc_resolve_assets_dir("from-arg"), path.expand("from-arg"))
})

test_that("set_mc_assets_dir() sets, creates, restores, and rejects non-paths", {
  withr::local_options(mc.assets_dir = NULL)
  withr::local_envvar(MC_ASSETS_DIR = NA)

  dir <- withr::local_tempdir()
  old <- set_mc_assets_dir(dir)
  expect_null(old)
  expect_equal(get_mc_assets_dir(), as.character(fs::path_expand(dir)))

  # the returned value is what restores the previous state
  set_mc_assets_dir(old)
  expect_equal(get_mc_assets_dir(), mc_default_assets_dir())

  # create the dir if missing (bad path fails here, not mid-download)
  nested <- file.path(withr::local_tempdir(), "a", "b")
  set_mc_assets_dir(nested)
  expect_true(dir.exists(nested))

  set_mc_assets_dir(NULL)
  expect_equal(get_mc_assets_dir(), mc_default_assets_dir())
  for (bad in list(5, c("a", "b"), "", NA_character_)) {
    expect_error(set_mc_assets_dir(bad))
  }
})

test_that("list_mc_assets() answers what exists, what is staged, what is reclaimable", {
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)

  # readable before anything is on disk, and without prompting or fetching
  before <- list_mc_assets()
  expect_equal(before$group_id, "FakeGroup")
  expect_false(before$downloaded)
  expect_equal(before$superseded, 0L)
  expect_gt(as.numeric(before$size), 0)
  expect_length(list.files(assets), 0)

  suppressMessages(download_mc_assets(ask = FALSE))
  writeLines(
    "an older pack",
    file.path(assets, sprintf("fakegroup-%s.qs2", strrep("a", 64)))
  )

  after <- list_mc_assets()
  expect_true(after$downloaded)
  expect_equal(after$superseded, 1L)
  expect_gt(as.numeric(after$superseded_size), 0)

  # sizes stay numeric, so a caller can total them
  expect_true(is.numeric(after$size))
  expect_error(list_mc_assets("NotAClockGroup"))
})

test_that("the shipped registry covers the three external groups", {
  expect_setequal(
    mc_external_groups(),
    c("SystemsAge", "PCClocks", "PCBrainAge")
  )
  for (gid in mc_external_groups()) {
    row <- mc_asset(gid)
    expect_equal(row$group_id, gid)
    expect_gt(row$size_bytes, 0)
  }
})

test_that("registry lookups reject unknown ids and resolve group sets", {
  expect_error(mc_asset("NotAClockGroup"))
  expect_error(mc_resolve_groups(c("PCClocks", "Nope")))
  expect_setequal(mc_resolve_groups("all"), mc_external_groups())
  expect_setequal(mc_resolve_groups(NULL), mc_external_groups())
})

test_that("load_mc_assets() refuses to fetch unprompted in a non-interactive session", {
  skip_if(interactive())
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)

  expect_error(load_mc_assets("FakeGroup"))
  expect_length(list.files(assets), 0)
})

test_that("download_mc_assets() fetches, verifies, and leaves no scratch files", {
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)

  paths <- suppressMessages(download_mc_assets(ask = FALSE))
  expect_equal(unname(basename(paths)), row$file)
  expect_true(file.exists(file.path(assets, row$file)))
  expect_false(any(grepl(".part", list.files(assets), fixed = TRUE)))

  file.remove(row$.src)
  expect_silent(suppressMessages(download_mc_assets(ask = FALSE)))
})

test_that("a download failure is reported with the URL and leaves nothing behind", {
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  file.remove(row$.src) # the "server" 404s

  expect_error(download_mc_assets(ask = FALSE))
  expect_length(list.files(assets), 0)
})

test_that("load_mc_assets() downloads missing packs on consent and returns a named registry", {
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)

  packs <- suppressMessages(load_mc_assets("FakeGroup", ask = FALSE))
  expect_named(packs, "FakeGroup")
  expect_equal(packs[["FakeGroup"]], row$.payload)
  expect_true(file.exists(file.path(assets, row$file)))
})

test_that("load_mc_assets() rejects a corrupt staged file via the qs2 checksum", {
  assets <- withr::local_tempdir()
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  writeLines("not a qs2 file", file.path(assets, row$file))
  expect_error(load_mc_assets("FakeGroup", from = assets))
})

test_that("load_mc_assets() with an explicit path is a closed set: no download, missing is fatal", {
  assets <- withr::local_tempdir() # deliberately empty
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)

  expect_error(load_mc_assets("FakeGroup", from = assets, ask = FALSE))
  expect_length(list.files(assets), 0)
})

test_that("load_mc_assets() resolves in-memory pack(s) without touching disk", {
  dir <- withr::local_tempdir()
  a <- fake_asset(dir, group = "GroupA")
  b <- fake_asset(dir, group = "GroupB")
  local_fake_registry(stats::setNames(list(a, b), c("GroupA", "GroupB")))

  expect_equal(
    load_mc_assets("GroupA", from = a$.payload)[["GroupA"]],
    a$.payload
  )
  expect_error(load_mc_assets("GroupB", from = a$.payload))
  expect_warning(
    res <- load_mc_assets("GroupA", from = list(a$.payload, b$.payload))
  )
  expect_named(res, "GroupA")
  expect_equal(res[["GroupA"]], a$.payload)
})

test_that("clear_mc_assets() removes staged packs only on explicit consent", {
  skip_if(interactive())
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)

  # nothing staged: reports and is a no-op
  expect_message(clear_mc_assets())

  suppressMessages(download_mc_assets(ask = FALSE))
  # unprompted deletion is refused non-interactively (file survives)
  expect_error(clear_mc_assets())
  expect_true(file.exists(file.path(assets, row$file)))

  removed <- suppressMessages(clear_mc_assets(ask = FALSE))
  expect_equal(basename(removed), row$file)
  expect_false(file.exists(file.path(assets, row$file)))
})

test_that("clear_mc_assets() reclaims superseded packs and spares everything else", {
  skip_if(interactive())
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  suppressMessages(download_mc_assets(ask = FALSE))

  # what an earlier sync left behind: same declared stem, older content hash
  superseded <- file.path(
    assets,
    sprintf("fakegroup-%s.qs2", strrep("a", 64))
  )
  writeLines("an older pack", superseded)

  # neither of these is ours: a foreign stem, and a file with no content address
  bystanders <- file.path(
    assets,
    c(sprintf("othergroup-%s.qs2", strrep("b", 64)), "notes.txt")
  )
  for (f in bystanders) {
    writeLines("keep me", f)
  }

  expect_named(mc_stale_files(), "FakeGroup")
  expect_equal(unname(basename(mc_stale_files())), basename(superseded))

  # clear means clear: the current pack and the superseded one both go
  removed <- suppressMessages(clear_mc_assets(ask = FALSE))
  expect_setequal(basename(removed), basename(c(superseded, row$file)))
  expect_false(file.exists(file.path(assets, row$file)))
  expect_false(file.exists(superseded))

  # but only ours -- a foreign stem and an uncontent-addressed file survive
  expect_true(all(file.exists(bystanders)))
  expect_length(mc_stale_files(), 0)
})

test_that("the delete prompt counts downloaded and superseded packs apart", {
  skip_if(interactive())
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  suppressMessages(download_mc_assets(ask = FALSE))
  writeLines(
    "an older pack",
    file.path(assets, sprintf("fakegroup-%s.qs2", strrep("a", 64)))
  )

  # non-interactive refusal is where the clear summary is observable
  expect_error(clear_mc_assets(), "1 downloaded pack and 1 superseded pack")
})

test_that("download -> load -> clear round trips and leaves the assets dir empty", {
  skip_if(interactive())
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  a <- fake_asset(withr::local_tempdir(), group = "GroupA")
  b <- fake_asset(withr::local_tempdir(), group = "GroupB")
  local_fake_registry(stats::setNames(list(a, b), c("GroupA", "GroupB")))

  paths <- suppressMessages(download_mc_assets(ask = FALSE))
  expect_true(all(file.exists(paths)))
  expect_false(any(grepl(".part", list.files(assets), fixed = TRUE)))

  # staged: loads from disk, needs no consent even with ask = TRUE
  packs <- load_mc_assets("all")
  expect_named(packs, c("GroupA", "GroupB"))
  expect_equal(packs[["GroupA"]], a$.payload)
  expect_equal(packs[["GroupB"]], b$.payload)

  removed <- suppressMessages(clear_mc_assets(ask = FALSE))
  expect_setequal(basename(removed), c(a$file, b$file))
  expect_false(any(file.exists(paths)))
  expect_length(mc_staged_files("all"), 0)
  expect_length(list.files(assets), 0)

  # really gone: the next load would have to download, so it refuses
  expect_error(load_mc_assets("all"))

  # an empty request stays empty (it is not "all")
  expect_length(load_mc_assets(character(0)), 0)
})

test_that("a non-path `from` errors instead of silently hitting the assets dir", {
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  suppressMessages(download_mc_assets(ask = FALSE))
  staged <- file.path(assets, row$file)
  expect_true(file.exists(staged))

  # a loaded pack names no directory -- the dir resolver must reject it
  expect_error(mc_resolve_assets_dir(row$.payload))
  expect_error(mc_resolve_assets_dir(5))
  expect_error(mc_resolve_assets_dir(c("a", "b")))
  expect_error(mc_resolve_assets_dir(""))
  expect_error(load_mc_assets("FakeGroup", from = 5))
  expect_true(file.exists(staged))
})

test_that("`ask` is a strict flag -- only FALSE consents", {
  skip_if(interactive())
  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  bad_flags <- list(NA, NULL, "yes", 1, c(TRUE, TRUE))

  for (bad in bad_flags) {
    expect_error(download_mc_assets(ask = bad))
    expect_error(load_mc_assets("FakeGroup", ask = bad))
  }
  expect_length(list.files(assets), 0) # nothing fetched under a bad flag

  suppressMessages(download_mc_assets(ask = FALSE))
  staged <- file.path(assets, row$file)
  for (bad in bad_flags) {
    expect_error(clear_mc_assets(ask = bad))
  }
  expect_true(file.exists(staged)) # and nothing deleted under one either
})

test_that("the real PCBrainAge release asset downloads and verifies", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not(
    nzchar(Sys.getenv("MC_TEST_NETWORK")),
    "set MC_TEST_NETWORK=1 to run live download tests"
  )

  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  packs <- suppressMessages(load_mc_assets("PCBrainAge", ask = FALSE))
  pack <- packs[["PCBrainAge"]]
  row <- mc_asset("PCBrainAge")
  expect_length(pack$cpgs, row$n_cpgs)
  expect_equal(nrow(pack$coefficient_matrix), row$n_cpgs)
})
