# External clock-data packs: content-addressed qs2 files fetched on demand.

MC_DEFAULT_RELEASE_REPO <- "hhp94/methylCIPHERv2"

# registry

mc_external_groups <- function() {
  assets <- mc_provenance[["external_assets"]]
  if (is.null(assets)) character(0) else names(assets)
}

# One provenance row for an external group.
mc_asset <- function(group_id) {
  row <- mc_provenance[["external_assets"]][[group_id]]
  if (is.null(row)) {
    stop(
      "Not an external clock group: ",
      group_id,
      "\nKnown groups: ",
      paste(mc_external_groups(), collapse = ", "),
      call. = FALSE
    )
  }
  row
}

# "all" or a vector of known group ids.
mc_resolve_groups <- function(groups) {
  if (is.null(groups) || identical(groups, "all") || !length(groups)) {
    return(mc_external_groups())
  }
  groups <- unique(as.character(groups))
  for (g in groups) {
    mc_asset(g)
  } # errors on any unknown id
  groups
}

# Public release-asset URL (option override for forks/testing).
mc_asset_url <- function(row) {
  repo <- getOption("mc.release_repo", MC_DEFAULT_RELEASE_REPO)
  sprintf(
    "https://github.com/%s/releases/download/%s/%s",
    repo,
    row[["release_tag"]],
    row[["file"]]
  )
}

# cache location

# CRAN-sanctioned per-user cache directory.
mc_default_cache_dir <- function() {
  path.expand(tools::R_user_dir("methylCIPHERv2", which = "cache"))
}

nz1 <- function(x) length(x) == 1L && !is.na(x) && nzchar(x)

# Active cache dir: assets arg, session option, env, then default.
mc_cache_dir <- function(assets = NULL) {
  if (nz1(assets)) {
    return(path.expand(assets))
  }
  opt <- getOption("mc.cache_dir")
  if (nz1(opt)) {
    return(path.expand(opt))
  }
  env <- Sys.getenv("MC_CACHE_DIR", unset = "")
  if (nz1(env)) {
    return(path.expand(env))
  }
  mc_default_cache_dir()
}

# helpers

mc_bytes <- function(x) {
  format(structure(as.numeric(x), class = "object_size"), units = "auto")
}

# Which packs are present in the cache (query only).
mc_cached_files <- function(groups = "all", assets = NULL) {
  groups <- mc_resolve_groups(groups)
  dir <- mc_cache_dir(assets)
  files <- vapply(
    groups,
    function(g) file.path(dir, mc_asset(g)[["file"]]),
    character(1)
  )
  files[file.exists(files)]
}

# download

# Fetch one pack: stage, validate qs2 checksum, atomic rename.
mc_fetch <- function(row, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  url <- mc_asset_url(row)
  dest <- file.path(dir, row[["file"]])
  tmp <- paste0(dest, ".part")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  # Large packs need a longer whole-transfer timeout.
  old_to <- options(timeout = max(getOption("timeout", 60L), 1800L))
  on.exit(options(old_to), add = TRUE)

  tryCatch(
    # quiet in tests / R CMD check; progress bar when interactive.
    utils::download.file(
      url,
      destfile = tmp,
      mode = "wb",
      quiet = !interactive()
    ),
    error = function(e) {
      stop(
        "Download failed for ",
        row[["group_id"]],
        ": ",
        conditionMessage(e),
        "\nURL: ",
        url,
        call. = FALSE
      )
    }
  )
  qs2::qs_read(tmp, validate_checksum = TRUE)
  file.rename(tmp, dest)
  dest
}

# Consent for downloading missing packs, or stop.
mc_consent <- function(rows, dir, ask) {
  if (!length(rows) || !isTRUE(ask)) {
    return(invisible(TRUE))
  }
  labels <- vapply(
    rows,
    function(r) {
      sprintf("%s (%s)", r[["group_id"]], mc_bytes(r[["size_bytes"]]))
    },
    character(1)
  )
  total <- mc_bytes(sum(vapply(
    rows,
    function(r) as.numeric(r[["size_bytes"]]),
    numeric(1)
  )))
  if (!interactive()) {
    stop(
      "Refusing to download ",
      paste(labels, collapse = ", "),
      " [",
      total,
      " total] without confirmation in a non-interactive session.\n",
      "Pass ask = FALSE to consent, or pre-stage the file(s) and point `assets` at them.",
      call. = FALSE
    )
  }
  ok <- utils::askYesNo(sprintf(
    "Download %d clock-data pack(s) [%s total] to\n  %s\n  %s\n?",
    length(rows),
    total,
    dir,
    paste(labels, collapse = "\n  ")
  ))
  if (!isTRUE(ok)) {
    stop(
      "Download declined for ",
      paste(
        vapply(rows, function(r) r[["group_id"]], character(1)),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Pre-fetch packs into the cache (skips already-present files).
mc_data_download <- function(groups = "all", assets = NULL, ask = TRUE) {
  groups <- mc_resolve_groups(groups)
  dir <- mc_cache_dir(assets)
  rows <- lapply(groups, mc_asset)
  files <- vapply(rows, function(r) file.path(dir, r[["file"]]), character(1))
  missing <- !file.exists(files)
  if (any(missing)) {
    mc_consent(rows[missing], dir, ask)
    for (i in which(missing)) {
      files[i] <- mc_fetch(rows[[i]], dir)
    }
  }
  invisible(files)
}

# load

# Canonicalize assets: NULL (open), cache-dir path, or loaded pack registry.
mc_canonicalize_assets <- function(assets) {
  if (is.null(assets)) {
    return(NULL)
  }
  if (is.character(assets)) {
    if (length(assets) != 1L || !nzchar(assets)) {
      stop("`assets` path must be a single non-empty string.", call. = FALSE)
    }
    return(assets)
  }
  is_pack <- function(x) is.list(x) && !is.null(x[["group_id"]])
  if (is_pack(assets)) {
    return(stats::setNames(list(assets), assets[["group_id"]]))
  }
  if (
    is.list(assets) &&
      length(assets) &&
      all(vapply(assets, is_pack, logical(1)))
  ) {
    return(stats::setNames(
      assets,
      vapply(assets, function(p) as.character(p[["group_id"]]), character(1))
    ))
  }
  stop(
    "`assets` must be NULL, a cache-dir path, a loaded pack, or a list of loaded packs.",
    call. = FALSE
  )
}

# Load packs for needed groups (open set may download; closed set never does).
load_mc_assets <- function(groups, assets = NULL, ask = TRUE) {
  groups <- unique(as.character(groups))
  groups <- groups[nzchar(groups)]
  if (!length(groups)) {
    return(stats::setNames(list(), character(0)))
  }
  for (g in groups) {
    mc_asset(g)
  } # errors on any unknown id

  canon <- mc_canonicalize_assets(assets)

  # Closed set from in-memory pack(s).
  if (is.list(canon)) {
    packs <- lapply(groups, function(g) {
      pack <- canon[[g]]
      if (is.null(pack)) {
        stop(
          "load_mc_assets(): external group '",
          g,
          "' is needed but not present in `assets` (closed set; no download).",
          call. = FALSE
        )
      }
      pack
    })
    extra <- setdiff(names(canon), groups)
    if (length(extra)) {
      warning(
        "load_mc_assets(): ignoring provided asset(s) not needed by the plan: ",
        paste(extra, collapse = ", "),
        call. = FALSE
      )
    }
    return(stats::setNames(packs, groups))
  }

  # Path (closed) or NULL (open): both read from a cache dir; only open may download.
  closed <- !is.null(canon)
  dir <- mc_cache_dir(canon)
  rows <- lapply(groups, mc_asset)
  files <- vapply(rows, function(r) file.path(dir, r[["file"]]), character(1))
  missing <- !file.exists(files)

  if (any(missing)) {
    if (closed) {
      stop(
        "load_mc_assets(): external group(s) ",
        paste(groups[missing], collapse = ", "),
        " not found in `assets` dir '",
        dir,
        "' (closed set; no download).",
        call. = FALSE
      )
    }
    mc_consent(rows[missing], dir, ask)
    for (i in which(missing)) {
      files[i] <- mc_fetch(rows[[i]], dir)
    }
  }

  # qs2's own checksum is the integrity guard; the filename is the content address.
  packs <- lapply(files, qs2::qs_read, validate_checksum = TRUE)
  stats::setNames(packs, groups)
}

# clear (stub; deletion is consent-gated)

# TODO: interactive delete flow. Never unlinks automatically.
clear_clock_cache <- function(groups = "all", assets = NULL) {
  files <- mc_cached_files(groups, assets)
  if (!length(files)) {
    message("No cached clock data to clear in ", mc_cache_dir(assets), ".")
    return(invisible(character(0)))
  }
  message(
    "Cached clock data (",
    mc_bytes(sum(file.size(files))),
    ") in ",
    mc_cache_dir(assets),
    ":\n  ",
    paste(basename(files), collapse = "\n  "),
    "\nDeletion is not yet wired up; remove these files manually for now."
  )
  invisible(files)
}
