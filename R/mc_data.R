# external clock-data assets (content-addressed qs2 packs)

# content-addressed tail of a pack filename: <stem>-<sha256>.qs2
MC_ASSET_SUFFIX <- "-[0-9a-f]{64}\\.qs2$"

mc_external_groups <- function() {
  assets <- mc_provenance[["external_assets"]]
  if (is.null(assets)) character(0) else names(assets)
}

# provenance row for one external group
mc_asset <- function(group_id) {
  row <- mc_provenance[["external_assets"]][[group_id]]
  if (is.null(row)) {
    cli::cli_abort(
      c(
        "{.val {group_id}} is not a known external clock group.",
        "i" = "Available groups: {.val {mc_external_groups()}}."
      ),
      call = NULL
    )
  }
  row
}

# "all" or a vector of known group ids, deduplicated.
mc_resolve_groups <- function(groups) {
  # empty selection selects nothing. only "all" means all.
  if (is.null(groups) || !length(groups)) {
    return(character(0))
  }
  # exact match only, like list_clocks(tag =). no partial matching.
  checkmate::assert_subset(groups, c("all", mc_external_groups()))
  groups <- unique(groups)
  if ("all" %in% groups) {
    return(mc_external_groups())
  }
  groups
}

# public release-asset URL
mc_asset_url <- function(row) {
  repo <- getOption("mc.release_repo", "hhp94/methylCIPHERv2")
  sprintf(
    "https://github.com/%s/releases/download/%s/%s",
    repo,
    row[["release_tag"]],
    row[["file"]]
  )
}

# default per-user assets dir (R_user_dir cache, never data)
mc_default_assets_dir <- function() {
  as.character(fs::path_expand(
    tools::R_user_dir("methylCIPHERv2", which = "cache")
  ))
}

# single non-empty path string
is_path_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

# active assets dir: `ext_data` arg, option, env, then default
mc_resolve_assets_dir <- function(ext_data = NULL) {
  if (!is.null(ext_data)) {
    if (!is_path_string(ext_data)) {
      cli::cli_abort(
        c(
          "{.arg ext_data} should be {.code NULL} or a single assets-dir path
           (got {.cls {class(ext_data)[[1L]]}} of length {length(ext_data)}).",
          "i" = "A loaded asset does not name a directory.",
          "i" = "Pass a loaded asset to {.fn load_mc_assets} instead."
        ),
        call = NULL
      )
    }
    return(as.character(fs::path_expand(ext_data)))
  }
  opt <- getOption("mc.assets_dir")
  if (is_path_string(opt)) {
    return(as.character(fs::path_expand(opt)))
  }
  env <- Sys.getenv("MC_ASSETS_DIR", unset = "")
  if (is_path_string(env)) {
    return(as.character(fs::path_expand(env)))
  }
  mc_default_assets_dir()
}

# the assets dir in effect for this session
#' Active Assets Directory
#'
#' Returns the directory that this session uses for clock assets.
#'
#' @details
#' The directory comes from the `mc.assets_dir` option, the `MC_ASSETS_DIR`
#' environment variable, or the default cache directory, in that order. Set
#' an override with [set_mc_assets_dir()].
#'
#' @returns A string. The absolute path to the assets directory in effect
#'   for this session.
#'
#' @seealso
#' - [set_mc_assets_dir()] to choose another directory.
#' - [list_mc_assets()] for the assets in the directory.
#' - [download_mc_assets()] to write the assets to disk.
#' - [load_mc_assets()] to read the assets into memory.
#' - [clear_mc_assets()] to delete the assets from disk.
#'
#' @examplesIf interactive()
#' get_mc_assets_dir()
#'
#' @export
get_mc_assets_dir <- function() {
  mc_resolve_assets_dir()
}

# set the session assets dir (NULL clears). returns previous value invisibly
#' Assets Directory Override
#'
#' Sets or clears the directory this session uses for clock assets.
#'
#' @param path A string. The new assets directory. Default is `NULL`, which
#'   clears the override.
#'
#' @details
#' `set_mc_assets_dir()` creates the directory if it does not exist yet. It
#' stops if the directory cannot be created or is not writable. Clearing the
#' override falls back to the `mc.assets_dir` option, the `MC_ASSETS_DIR`
#' environment variable, or the default cache directory, in that order.
#'
#' The return value restores the previous state exactly, in the manner of
#' [setwd()]. It is the override that was in place, and it is `NULL` when no
#' override was set. Passing `NULL` back clears the override, rather than
#' pinning it to the directory that happened to be in effect.
#'
#' For the directory in effect, which is always a path, call
#' [get_mc_assets_dir()].
#'
#' @returns A string. The override that was in place, returned invisibly, or
#'   `NULL` when no override was set.
#'
#' @seealso
#' - [get_mc_assets_dir()] for the directory in effect.
#' - [list_mc_assets()] for the assets in the directory.
#' - [download_mc_assets()] to write the assets to disk.
#' - [load_mc_assets()] to read the assets into memory.
#' - [clear_mc_assets()] to delete the assets from disk.
#'
#' @examplesIf interactive()
#' old <- set_mc_assets_dir(tempdir())
#' get_mc_assets_dir()
#' set_mc_assets_dir(old)
#'
#' @export
set_mc_assets_dir <- function(path = NULL) {
  # returns the previous override (NULL means no override was set).
  old <- getOption("mc.assets_dir")
  if (is.null(path)) {
    options(mc.assets_dir = NULL)
    return(invisible(old))
  }
  if (!is_path_string(path)) {
    cli::cli_abort(
      "{.arg path} should be {.code NULL} or a single non-empty string
       (got {.cls {class(path)[[1L]]}} of length {length(path)}).",
      call = NULL
    )
  }
  dir <- as.character(fs::path_expand(path))
  tryCatch(
    fs::dir_create(dir),
    error = function(e) {
      cli::cli_abort(
        c(
          "The assets directory {.path {dir}} cannot be created.
           {conditionMessage(e)}.",
          "i" = "Pass a writable path to {.arg path}."
        ),
        call = NULL
      )
    }
  )
  if (!isTRUE(unname(fs::file_access(dir, "write")))) {
    cli::cli_abort(
      c(
        "The assets directory {.path {dir}} is not writable.",
        "i" = "Check the permissions on that directory.",
        "i" = "To use another path, pass it to {.arg path}."
      ),
      call = NULL
    )
  }
  options(mc.assets_dir = dir)
  invisible(old)
}

# staged paths for provenance rows, in order
mc_pack_paths <- function(dir, rows) {
  if (!length(rows)) {
    return(character(0))
  }
  files <- vapply(rows, function(r) r[["file"]], character(1))
  as.character(fs::path(dir, files))
}

# aligned label/size lines for cli_verbatim. consent prompt shows a capped head.
mc_manifest_lines <- function(labels, sizes) {
  if (!length(labels)) {
    return(character(0))
  }
  sizes <- as.numeric(sizes)
  cells <- format(fs::fs_bytes(c(sizes, sum(sizes))))
  w <- max(nchar(labels), nchar("total"))
  keep <- seq_len(min(length(labels), MC_MSG_CAP))
  out <- sprintf("  %-*s  %s", w, labels[keep], cells[keep])
  n_more <- length(labels) - length(keep)
  if (n_more) {
    out <- c(out, cli::format_inline("  and {n_more} more asset{?s}"))
  }
  if (length(labels) > 1L) {
    out <- c(out, sprintf("  %-*s  %s", w, "total", cells[[length(cells)]]))
  }
  out
}

# one bullet per pack. brace in a file name cannot become a cli template.
mc_manifest_bullets <- function(labels, sizes) {
  sizes <- as.numeric(sizes)
  human <- trimws(format(fs::fs_bytes(sizes)))
  items <- capped_bullets(seq_along(labels), function(i) {
    vapply(
      i,
      function(k) cli::format_inline("{.val {labels[[k]]}} ({human[[k]]})"),
      character(1L)
    )
  })
  if (length(labels) > 1L) {
    total <- trimws(format(fs::fs_bytes(sum(sizes))))
    items <- c(items, "*" = cli::format_inline("total ({total})"))
  }
  items
}

# cli renders the context, askYesNo asks one short line
mc_ask_yes_no <- function(header, labels, sizes, dir, question) {
  cli::cli_inform(c("i" = "{header}", "i" = "Assets dir: {.path {dir}}"))
  cli::cli_verbatim(mc_manifest_lines(labels, sizes))
  utils::flush.console()
  isTRUE(utils::askYesNo(question))
}

# stage, validate, atomic rename
mc_fetch <- function(row, dir) {
  fs::dir_create(dir)
  url <- mc_asset_url(row)
  dest <- fs::path(dir, row[["file"]])
  tmp <- paste0(as.character(dest), ".part")
  on.exit(if (fs::file_exists(tmp)) fs::file_delete(tmp), add = TRUE)

  old_to <- options(timeout = max(getOption("timeout", 60L), 1800L))
  on.exit(options(old_to), add = TRUE)

  status <- tryCatch(
    utils::download.file(
      url,
      destfile = tmp,
      mode = "wb",
      quiet = !interactive()
    ),
    error = function(e) {
      cli::cli_abort(
        c(
          "The asset for {.val {row[['group_id']]}} did not download.
           {conditionMessage(e)}.",
          "i" = "URL: {.url {url}}",
          "i" = "Run {.fn download_mc_assets} again to retry.",
          "i" = "Or fetch that URL by hand and point {.arg ext_data} at the
                 directory."
        ),
        call = NULL
      )
    }
  )
  if (!identical(as.integer(status), 0L)) {
    cli::cli_abort(
      c(
        "The asset for {.val {row[['group_id']]}} did not download
         (status {status}).",
        "i" = "URL: {.url {url}}",
        "i" = "Run {.fn download_mc_assets} again to retry.",
        "i" = "Or fetch that URL by hand and point {.arg ext_data} at the
               directory."
      ),
      call = NULL
    )
  }
  payload <- qs2::qs_read(tmp, validate_checksum = TRUE)
  fs::file_move(tmp, dest)
  list(path = as.character(dest), payload = payload)
}

# consent for downloading missing packs
mc_consent <- function(rows, dir, ask) {
  if (!length(rows) || !ask) {
    return(invisible(TRUE))
  }
  ids <- vapply(rows, function(r) r[["group_id"]], character(1))
  sizes <- vapply(rows, function(r) as.numeric(r[["size_bytes"]]), numeric(1))

  if (!interactive()) {
    cli::cli_abort(
      c(
        "{length(rows)} clock-data asset{?s} cannot be downloaded in a
         non-interactive session.",
        "i" = "Assets dir: {.path {dir}}",
        mc_manifest_bullets(ids, sizes),
        "i" = "Pass {.code ask = FALSE} to allow the download.",
        "i" = "Or put the file{cli::qty(length(rows))}{?s} in that directory
               yourself, then point {.arg ext_data} at it."
      ),
      call = NULL
    )
  }
  ok <- mc_ask_yes_no(
    header = cli::format_inline(
      "Download {length(rows)} clock-data asset{?s}:"
    ),
    labels = ids,
    sizes = sizes,
    dir = dir,
    question = sprintf("Download %d asset(s)?", length(rows))
  )
  if (!ok) {
    cli::cli_abort(
      c(
        "You cancelled the download of {.val {ids}}.",
        "i" = "Run {.fn download_mc_assets} again to answer the prompt a
               second time.",
        "i" = "Or point {.arg ext_data} at a directory that holds the
               asset{cli::qty(ids)}{?s}."
      ),
      call = NULL
    )
  }
  invisible(TRUE)
}

# pre-fetch packs into the assets dir
#' Clock Asset Downloads
#'
#' Downloads the clock assets for one or more groups into the assets
#' directory.
#'
#' @inheritParams mc-params
#'
#' @details
#' Only an asset that is missing from the assets directory is fetched.
#' `download_mc_assets()` asks for consent before a download, refuses in a
#' non-interactive session, and treats `ask = FALSE` as consent to proceed
#' without asking.
#'
#' @returns A character vector. The path to each requested asset file, named
#'   by group id. Returned invisibly.
#'
#' @seealso
#' - [list_mc_assets()] for the assets already on disk.
#' - [load_mc_assets()] to read the assets into memory.
#' - [clear_mc_assets()] to delete the assets from disk.
#' - [get_mc_assets_dir()] for the directory in effect.
#' - [set_mc_assets_dir()] to choose another directory.
#'
#' @examplesIf interactive()
#' download_mc_assets("SystemsAge")
#' download_mc_assets("SystemsAge", ask = FALSE)
#'
#' @export
download_mc_assets <- function(groups = "all", ask = TRUE) {
  checkmate::assert_flag(ask)
  groups <- mc_resolve_groups(groups)
  dir <- mc_resolve_assets_dir()
  rows <- lapply(groups, mc_asset)
  files <- mc_pack_paths(dir, rows)
  missing <- unname(!fs::file_exists(files))
  if (any(missing)) {
    mc_consent(rows[missing], dir, ask)
    for (i in which(missing)) {
      files[[i]] <- mc_fetch(rows[[i]], dir)[["path"]]
    }
  }
  invisible(stats::setNames(files, groups))
}

# currently declared packs present in the assets dir, named by group id
mc_staged_files <- function(groups = "all") {
  groups <- mc_resolve_groups(groups)
  dir <- mc_resolve_assets_dir()
  files <- mc_pack_paths(dir, lapply(groups, mc_asset))
  files <- stats::setNames(files, groups)
  files[fs::file_exists(files)]
}

# content-address stem from the declared filename, or NA
mc_asset_stem <- function(row) {
  file <- as.character(row[["file"]])
  stem <- sub(MC_ASSET_SUFFIX, "", file)
  if (identical(stem, file)) NA_character_ else stem
}

# superseded packs left by a moved payload_hash (reclaim only)
mc_stale_files <- function(groups = "all") {
  groups <- mc_resolve_groups(groups)
  dir <- mc_resolve_assets_dir()
  empty <- stats::setNames(character(0), character(0))
  if (!fs::dir_exists(dir)) {
    return(empty)
  }
  on_disk <- as.character(fs::path_file(
    fs::dir_ls(dir, type = "file", recurse = FALSE)
  ))
  if (!length(on_disk)) {
    return(empty)
  }
  addressed <- on_disk[grepl(MC_ASSET_SUFFIX, on_disk)]
  out <- character(0)
  labels <- character(0)
  for (g in groups) {
    row <- mc_asset(g)
    stem <- mc_asset_stem(row)
    if (is.na(stem)) {
      next
    }
    hit <- addressed[startsWith(addressed, paste0(stem, "-"))]
    hit <- setdiff(hit, as.character(row[["file"]]))
    if (!length(hit)) {
      next
    }
    out <- c(out, as.character(fs::path(dir, hit)))
    labels <- c(labels, rep(g, length(hit)))
  }
  stats::setNames(out, labels)
}

# display label for a superseded pack
mc_stale_labels <- function(stale) {
  if (!length(stale)) {
    return(character(0))
  }
  paste0(names(stale), " (superseded)")
}

# browsable table: what exists, how big, what is downloaded, what is reclaimable
#' Clock Asset Inventory
#'
#' Lists the declared clock assets and their status in the assets
#' directory.
#'
#' @inheritParams mc-params
#'
#' @returns A data.frame. One row for each requested group, with its clock
#'   count, CpG count, asset size, download status, and the count and total
#'   size of its superseded assets.
#'
#' @seealso
#' - [list_clocks()] for the clocks in the catalog.
#' - [list_clock_tags()] for the tags a `clocks` value accepts.
#' - [clock_cpgs()] for the CpGs a set of clocks needs.
#' - [download_mc_assets()] to write the assets to disk.
#' - [load_mc_assets()] to read the assets into memory.
#' - [clear_mc_assets()] to delete the assets from disk.
#' - [get_mc_assets_dir()] for the directory in effect.
#' - [set_mc_assets_dir()] to choose another directory.
#'
#' @examplesIf interactive()
#' list_mc_assets()
#' list_mc_assets("SystemsAge")
#'
#' @export
list_mc_assets <- function(groups = "all") {
  groups <- mc_resolve_groups(groups)
  dir <- mc_resolve_assets_dir()
  rows <- lapply(groups, mc_asset)
  stale <- mc_stale_files(groups)
  stale_size <- if (length(stale)) {
    as.numeric(fs::file_size(stale))
  } else {
    numeric(0)
  }
  from_row <- function(field, mode) {
    vapply(rows, function(r) as.vector(r[[field]], mode), vector(mode, 1L))
  }

  out <- data.frame(
    group_id = groups,
    n_clocks = from_row("n_clocks", "integer"),
    n_cpgs = from_row("n_cpgs", "integer"),
    downloaded = unname(fs::file_exists(mc_pack_paths(dir, rows))),
    superseded = vapply(
      groups,
      function(g) sum(names(stale) == g),
      integer(1L),
      USE.NAMES = FALSE
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # assign after data.frame() so fs_bytes class survives
  out[["size"]] <- fs::fs_bytes(from_row("size_bytes", "double"))
  out[["superseded_size"]] <- fs::fs_bytes(vapply(
    groups,
    function(g) sum(stale_size[names(stale) == g]),
    numeric(1L),
    USE.NAMES = FALSE
  ))
  out[c(
    "group_id",
    "n_clocks",
    "n_cpgs",
    "size",
    "downloaded",
    "superseded",
    "superseded_size"
  )]
}

# ext_data: NULL (open), assets-dir path, or loaded pack registry
mc_canonicalize_ext_data <- function(ext_data) {
  if (is.null(ext_data)) {
    return(NULL)
  }
  # a path is validated once, by mc_resolve_assets_dir()
  if (is.character(ext_data)) {
    return(ext_data)
  }
  is_pack <- function(x) is.list(x) && !is.null(x[["group_id"]])
  if (is_pack(ext_data)) {
    return(stats::setNames(list(ext_data), ext_data[["group_id"]]))
  }
  if (
    is.list(ext_data) &&
      length(ext_data) &&
      all(vapply(ext_data, is_pack, logical(1)))
  ) {
    return(stats::setNames(
      ext_data,
      vapply(ext_data, function(p) as.character(p[["group_id"]]), character(1))
    ))
  }
  cli::cli_abort(
    "{.arg ext_data} should be {.code NULL}, an assets-dir path, a loaded
     asset, or a list of loaded assets.",
    call = NULL
  )
}

# load packs for needed groups into memory
#' Loaded Clock Assets
#'
#' Reads the assets for one or more clock groups into memory.
#'
#' @inheritParams mc-params
#' @param groups A character vector. The asset groups to load. One or more of
#'   `"PCBrainAge"`, `"PCClocks"`, `"SystemsAge"` and `"Zhang2019"`, or
#'   `"all"` for every group. Repeated values are ignored, and an empty vector
#'   loads nothing.
#'
#' @inheritSection mc-params The assets directory
#'
#' @details
#' `load_mc_assets()` asks for consent before a download, refuses in a
#' non-interactive session, and treats `ask = FALSE` as consent to proceed
#' without asking.
#'
#' An asset in `ext_data` for a group that was not requested is not used, and
#' `load_mc_assets()` warns about it.
#'
#' @returns A named list. It holds the loaded asset for each requested group,
#'   in the order of `groups`.
#'
#' @seealso
#' - [list_mc_assets()] for the assets already on disk.
#' - [download_mc_assets()] to write the assets to disk.
#' - [clear_mc_assets()] to delete the assets from disk.
#' - [get_mc_assets_dir()] for the directory in effect.
#' - [set_mc_assets_dir()] to choose another directory.
#'
#' @examplesIf interactive()
#' assets <- load_mc_assets("SystemsAge")
#' names(assets)
#'
#' @export
load_mc_assets <- function(groups, ext_data = NULL, ask = TRUE) {
  checkmate::assert_flag(ask)
  # one reading of `groups` for every verb. empty is the common case.
  groups <- mc_resolve_groups(groups)
  if (!length(groups)) {
    return(stats::setNames(list(), character(0)))
  }
  rows <- lapply(groups, mc_asset)

  canon <- mc_canonicalize_ext_data(ext_data)

  if (is.list(canon)) {
    packs <- lapply(groups, function(g) {
      pack <- canon[[g]]
      if (is.null(pack)) {
        cli::cli_abort(
          c(
            "The {.val {g}} asset is not in {.arg ext_data}.",
            "i" = "A list of loaded assets is fixed, so no asset is downloaded.",
            "i" = "Add the {.val {g}} asset to {.arg ext_data}.",
            "i" = "Or pass an assets directory path, or {.code NULL}, to allow a
                   download."
          ),
          call = NULL
        )
      }
      pack
    })
    extra <- setdiff(names(canon), groups)
    if (length(extra)) {
      cli::cli_warn(
        c(
          "{length(extra)} asset{?s} in {.arg ext_data}
           {cli::qty(extra)}{?is/are} not used:
           {.val {capped_vals(extra)}}.",
          "i" = "{.fn load_mc_assets} reads only the
                 asset{cli::qty(extra)}{?s} for the groups you asked for.",
          "i" = "Run {.fn list_mc_assets} to see the declared groups."
        ),
        call = NULL
      )
    }
    return(stats::setNames(packs, groups))
  }

  # path = closed set (no download), NULL = open set
  closed <- !is.null(canon)
  dir <- mc_resolve_assets_dir(canon)
  files <- mc_pack_paths(dir, rows)
  missing <- unname(!fs::file_exists(files))
  packs <- vector("list", length(groups))

  if (any(missing)) {
    if (closed) {
      gone <- groups[missing]
      cli::cli_abort(
        c(
          "The {.val {gone}} asset{cli::qty(gone)}{?s}
           {cli::qty(gone)}{?is/are} not in {.path {dir}}.",
          "i" = "A path in {.arg ext_data} is fixed, so no asset is
                 downloaded.",
          "i" = "Put the missing file{cli::qty(gone)}{?s} in that directory.",
          "i" = "Or pass {.code ext_data = NULL} to allow a download."
        ),
        call = NULL
      )
    }
    mc_consent(rows[missing], dir, ask)
    for (i in which(missing)) {
      got <- mc_fetch(rows[[i]], dir)
      files[[i]] <- got[["path"]]
      packs[[i]] <- got[["payload"]]
    }
  }

  unread <- vapply(packs, is.null, logical(1))
  packs[unread] <- lapply(files[unread], qs2::qs_read, validate_checksum = TRUE)
  stats::setNames(packs, groups)
}

# count phrase for the clear prompt (each part binds its own plural)
mc_delete_summary <- function(n_downloaded, n_stale) {
  parts <- character(0)
  if (n_downloaded) {
    parts <- c(parts, cli::format_inline("{n_downloaded} downloaded asset{?s}"))
  }
  if (n_stale) {
    parts <- c(parts, cli::format_inline("{n_stale} superseded asset{?s}"))
  }
  paste(parts, collapse = " and ")
}

# consent for deleting packs
mc_consent_delete <- function(files, dir, ask, n_stale = 0L) {
  if (!ask) {
    return(TRUE)
  }
  sizes <- as.numeric(fs::file_size(files))
  labels <- names(files)
  if (is.null(labels)) {
    labels <- as.character(fs::path_file(files))
  }
  what <- mc_delete_summary(length(files) - n_stale, n_stale)

  if (!interactive()) {
    cli::cli_abort(
      c(
        "{what} cannot be deleted in a non-interactive session.",
        "i" = "Assets dir: {.path {dir}}",
        mc_manifest_bullets(labels, sizes),
        "i" = "Pass {.code ask = FALSE} to allow the deletion."
      ),
      call = NULL
    )
  }
  mc_ask_yes_no(
    header = paste0("Delete ", what, ":"),
    labels = labels,
    sizes = sizes,
    dir = dir,
    question = sprintf("Delete %d asset(s)?", length(files))
  )
}

# remove every pack this package put in the assets dir
#' Clock Asset Removal
#'
#' Deletes the clock assets for one or more groups from the assets
#' directory.
#'
#' @inheritParams mc-params
#'
#' @details
#' `clear_mc_assets()` removes both the currently declared assets and every
#' superseded asset left behind by an earlier sync, for the requested groups.
#' It asks for consent before deleting, refuses in a non-interactive
#' session, and treats `ask = FALSE` as consent to proceed without asking.
#' If the assets directory holds no asset for the requested groups, no file
#' is removed.
#'
#' @returns A character vector. The path to each deleted asset file, returned
#'   invisibly. Empty when no file was removed.
#'
#' @seealso
#' - [list_mc_assets()] for the assets on disk and their size.
#' - [download_mc_assets()] to write the assets to disk.
#' - [load_mc_assets()] to read the assets into memory.
#' - [get_mc_assets_dir()] for the directory in effect.
#' - [set_mc_assets_dir()] to choose another directory.
#'
#' @examplesIf interactive()
#' clear_mc_assets("SystemsAge")
#'
#' @export
clear_mc_assets <- function(groups = "all", ask = TRUE) {
  checkmate::assert_flag(ask)
  dir <- mc_resolve_assets_dir()
  downloaded <- mc_staged_files(groups)
  stale <- mc_stale_files(groups)
  files <- c(downloaded, stats::setNames(stale, mc_stale_labels(stale)))
  if (!length(files)) {
    cli::cli_inform(c(
      "The assets directory {.path {dir}} holds no clock assets.",
      "i" = "To clear a different directory, set it with
             {.fn set_mc_assets_dir} first."
    ))
    return(invisible(character(0)))
  }
  if (!mc_consent_delete(files, dir, ask, n_stale = length(stale))) {
    cli::cli_inform("No file was removed from {.path {dir}}.")
    return(invisible(character(0)))
  }
  freed <- sum(as.numeric(fs::file_size(files)))
  # fs::file_delete() throws on failure
  fs::file_delete(files)
  cli::cli_inform(
    "{mc_delete_summary(length(downloaded), length(stale))}
     ({format(fs::fs_bytes(freed))}) {cli::qty(length(files))}{?was/were}
     deleted from {.path {dir}}."
  )
  invisible(files)
}
