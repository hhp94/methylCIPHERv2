# external clock-data asset tests (file://, live network opt-in)

# digest is a Suggests-only dep (content-addresses fake pack filenames)
skip_if_not_installed("digest")

# silence asset-verb chatter (cli + download.file stderr) when it is not under test.
quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(expr), type = "message")
  out
}

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

# standing setup: temp assets dir, one staged fake pack per group, mocked registry.
local_assets <- function(
  groups = "FakeGroup",
  set_option = TRUE,
  .env = parent.frame()
) {
  assets <- withr::local_tempdir(.local_envir = .env)
  if (set_option) {
    withr::local_options(mc.assets_dir = assets, .local_envir = .env)
  }
  rows <- lapply(groups, function(g) {
    fake_asset(withr::local_tempdir(.local_envir = .env), group = g)
  })
  names(rows) <- groups
  local_fake_registry(if (length(rows) == 1L) rows[[1L]] else rows, .env = .env)
  list(dir = assets, row = rows[[1L]], rows = rows)
}

# what an earlier sync left behind: same declared stem, older content hash
stage_superseded <- function(assets, group = "FakeGroup") {
  path <- file.path(
    assets,
    sprintf("%s-%s.qs2", tolower(group), strrep("a", 64))
  )
  writeLines("an older pack", path)
  path
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
  expect_equal(
    mc_resolve_assets_dir("ext-data-arg"),
    path.expand("ext-data-arg")
  )
})

test_that("set_mc_assets_dir() sets, creates, restores, and rejects non-paths", {
  withr::local_options(mc.assets_dir = NULL)
  withr::local_envvar(MC_ASSETS_DIR = NA)

  dir <- withr::local_tempdir()
  old <- set_mc_assets_dir(dir)
  # the previous override, not the resolved dir: there was none
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

  # set-then-restore must return the previous override, not the resolved dir.
  env_dir <- withr::local_tempdir()
  withr::local_envvar(MC_ASSETS_DIR = env_dir)
  was <- set_mc_assets_dir(withr::local_tempdir())
  expect_null(was)
  set_mc_assets_dir(was)
  expect_equal(get_mc_assets_dir(), as.character(fs::path_expand(env_dir)))
})

test_that("list_mc_assets() answers what is staged without fetching or prompting", {
  skip_on_cran()
  fx <- local_assets()
  assets <- fx$dir

  before <- list_mc_assets()
  expect_equal(before$group_id, "FakeGroup")
  expect_false(before$downloaded)
  expect_equal(before$superseded, 0L)
  expect_length(list.files(assets), 0) # read-only: nothing fetched

  quietly(download_mc_assets(ask = FALSE))
  stage_superseded(assets)

  after <- list_mc_assets()
  expect_true(after$downloaded)
  expect_equal(after$superseded, 1L)
  expect_true(is.numeric(after$size)) # so a caller can total them
  expect_error(list_mc_assets("NotAClockGroup"))
})

test_that("the shipped registry covers exactly the groups holding external clocks", {
  skip_on_cran()
  # group is external when any clock of it is.
  declared <- unique(unlist(lapply(mc_catalog, function(e) {
    if (isTRUE(e[["external_group"]])) e[["group_id"]] else NULL
  })))
  expect_setequal(mc_external_groups(), declared)
  expect_gt(length(declared), 0)
})

test_that("a group token is matched exactly, and an empty selection is never all", {
  skip_on_cran()
  expect_error(mc_asset("NotAClockGroup"))
  expect_error(mc_resolve_groups(c("PCClocks", "Nope")))
  expect_error(mc_resolve_groups("Sys")) # an abbreviation is not a group id

  # several allowed and deduplicated, and "all" absorbs the rest
  expect_equal(mc_resolve_groups(c("PCClocks", "PCClocks")), "PCClocks")
  expect_setequal(mc_resolve_groups(c("PCClocks", "all")), mc_external_groups())

  # empty selection selects nothing. never means every group.
  expect_equal(mc_resolve_groups(NULL), character(0))
  expect_equal(mc_resolve_groups(character(0)), character(0))
  expect_length(load_mc_assets(character(0)), 0L)
})

test_that("a consented fetch stages the pack, verifies it, and leaves no scratch", {
  skip_on_cran()
  fx <- local_assets()
  assets <- fx$dir
  row <- fx$row

  quietly(download_mc_assets(ask = FALSE))
  expect_true(file.exists(file.path(assets, row$file)))
  expect_false(any(grepl(".part", list.files(assets), fixed = TRUE)))

  # already staged: no second fetch, so the "server" can go away
  file.remove(row$.src)
  packs <- quietly(load_mc_assets("FakeGroup", ask = FALSE))
  expect_named(packs, "FakeGroup")
  expect_equal(packs[["FakeGroup"]], row$.payload)
})

test_that("a failed fetch leaves nothing behind", {
  skip_on_cran()
  fx <- local_assets()
  assets <- fx$dir
  file.remove(fx$row$.src) # the "server" 404s

  # download.file() warns then returns failed status. mc_fetch() turns that into the abort.
  expect_error(suppressWarnings(quietly(download_mc_assets(ask = FALSE))))
  expect_length(list.files(assets), 0)
})

test_that("load_mc_assets() rejects a corrupt staged file via the qs2 checksum", {
  skip_on_cran()
  assets <- withr::local_tempdir()
  row <- fake_asset(withr::local_tempdir())
  local_fake_registry(row)
  writeLines("not a qs2 file", file.path(assets, row$file))
  expect_error(load_mc_assets("FakeGroup", ext_data = assets))
})

test_that("an explicit ext_data is a closed set, and a non-path is an error", {
  fx <- local_assets()
  assets <- fx$dir
  row <- fx$row
  empty <- withr::local_tempdir()

  # a path names a closed set. missing pack is fatal, never downloaded.
  expect_error(load_mc_assets("FakeGroup", ext_data = empty, ask = FALSE))
  expect_length(list.files(empty), 0)
  expect_length(list.files(assets), 0)

  # a value that is not a path must not silently fall back to the assets dir
  quietly(download_mc_assets(ask = FALSE))
  staged <- file.path(assets, row$file)
  expect_error(load_mc_assets("FakeGroup", ext_data = 5))
  expect_true(file.exists(staged))
})

test_that("load_mc_assets() resolves in-memory pack(s) without touching disk", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  a <- fake_asset(dir, group = "GroupA")
  b <- fake_asset(dir, group = "GroupB")
  local_fake_registry(stats::setNames(list(a, b), c("GroupA", "GroupB")))

  expect_equal(
    load_mc_assets("GroupA", ext_data = a$.payload)[["GroupA"]],
    a$.payload
  )
  expect_error(load_mc_assets("GroupB", ext_data = a$.payload))
  expect_warning(
    res <- load_mc_assets("GroupA", ext_data = list(a$.payload, b$.payload))
  )
  expect_named(res, "GroupA")

  # a loaded pack names no directory, so the dir resolver must refuse one
  expect_error(mc_resolve_assets_dir(a$.payload))
})

test_that("clear_mc_assets() reclaims the current and superseded packs only", {
  skip_if(interactive())
  fx <- local_assets()
  assets <- fx$dir
  row <- fx$row

  # nothing staged: reports and is a no-op
  expect_message(clear_mc_assets())

  quietly(download_mc_assets(ask = FALSE))
  superseded <- stage_superseded(assets)

  # neither of these is ours: a foreign stem, and a file with no content address
  bystanders <- file.path(
    assets,
    c(sprintf("othergroup-%s.qs2", strrep("b", 64)), "notes.txt")
  )
  for (f in bystanders) {
    writeLines("keep me", f)
  }

  # clear means clear: the current pack and the superseded one both go
  removed <- suppressMessages(clear_mc_assets(ask = FALSE))
  expect_length(removed, 2L)
  expect_false(file.exists(file.path(assets, row$file)))
  expect_false(file.exists(superseded))

  # but only ours -- a foreign stem and an uncontent-addressed file survive
  expect_true(all(file.exists(bystanders)))
})

test_that("the consent gate fails closed -- only ask = FALSE moves bytes", {
  skip_if(interactive())
  fx <- local_assets()
  assets <- fx$dir
  row <- fx$row

  # unprompted, a non-interactive session cannot answer, so nothing is fetched
  expect_error(load_mc_assets("FakeGroup"))
  expect_length(list.files(assets), 0)

  bad_flags <- list(NA, NULL, "yes", 1, c(TRUE, TRUE))
  for (bad in bad_flags) {
    expect_error(download_mc_assets(ask = bad))
    expect_error(load_mc_assets("FakeGroup", ask = bad))
  }
  expect_length(list.files(assets), 0) # nothing fetched under a bad flag

  quietly(download_mc_assets(ask = FALSE))
  staged <- file.path(assets, row$file)
  expect_error(clear_mc_assets()) # deletion is refused unprompted too
  for (bad in bad_flags) {
    expect_error(clear_mc_assets(ask = bad))
  }
  expect_true(file.exists(staged)) # and nothing deleted under either
})

test_that("the real PCBrainAge release asset downloads and verifies", {
  skip_on_cran()
  # opt-in flag gates first. skip_if_offline() is a live DNS lookup.
  skip_if_not(
    nzchar(Sys.getenv("MC_TEST_NETWORK")),
    "set MC_TEST_NETWORK=1 to run live download tests"
  )
  skip_if_offline()

  assets <- withr::local_tempdir()
  withr::local_options(mc.assets_dir = assets)
  packs <- quietly(load_mc_assets("PCBrainAge", ask = FALSE))
  pack <- packs[["PCBrainAge"]]
  row <- mc_asset("PCBrainAge")
  expect_length(pack$cpgs, row$n_cpgs)
  expect_equal(nrow(pack$coefficient_matrix), row$n_cpgs)
})
