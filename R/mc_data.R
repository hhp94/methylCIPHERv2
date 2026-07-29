# external clock-data assets (content-addressed qs2 packs)
MC_DEFAULT_RELEASE_REPO <- "hhp94/methylCIPHERv2"

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
        "{.val {group_id}} is not an external clock group.",
        "i" = "Known groups: {.val {mc_external_groups()}}."
      ),
      call = NULL
    )
  }
  row
}

# "all" or a vector of known group ids
mc_resolve_groups <- function(groups) {
  if (is.null(groups) || identical(groups, "all") || !length(groups)) {
    return(mc_external_groups())
  }
  groups <- unique(as.character(groups))
  for (g in groups) {
    mc_asset(g)
  }
  groups
}

# public release-asset URL
mc_asset_url <- function(row) {
  repo <- getOption("mc.release_repo", MC_DEFAULT_RELEASE_REPO)
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
          "{.arg ext_data} must be {.code NULL} or a single assets-dir path
           (got {class(ext_data)[[1L]]} of length {length(ext_data)}).",
          "i" = "A loaded pack names no directory -- pass packs only to
                 {.fn load_mc_assets}."
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
#' @export
get_mc_assets_dir <- function() {
  mc_resolve_assets_dir()
}

# set the session assets dir (NULL clears). returns previous value invisibly
#' @export
set_mc_assets_dir <- function(path = NULL) {
  old <- getOption("mc.assets_dir")
  if (is.null(path)) {
    options(mc.assets_dir = NULL)
    return(invisible(old))
  }
  if (!is_path_string(path)) {
    cli::cli_abort(
      "{.arg path} must be {.code NULL} or a single non-empty string
       (got {class(path)[[1L]]} of length {length(path)}).",
      call = NULL
    )
  }
  dir <- as.character(fs::path_expand(path))
  tryCatch(
    fs::dir_create(dir),
    error = function(e) {
      cli::cli_abort(
        "Can't create assets dir {.path {dir}}: {conditionMessage(e)}.",
        call = NULL
      )
    }
  )
  if (!isTRUE(unname(fs::file_access(dir, "write")))) {
    cli::cli_abort(
      "Assets dir {.path {dir}} is not writable.",
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

# aligned label/size lines (cli_verbatim only)
mc_manifest_lines <- function(labels, sizes) {
  if (!length(labels)) {
    return(character(0))
  }
  sizes <- as.numeric(sizes)
  cells <- format(fs::fs_bytes(c(sizes, sum(sizes))))
  w <- max(nchar(labels), nchar("total"))
  out <- sprintf("  %-*s  %s", w, labels, cells[seq_along(labels)])
  if (length(labels) > 1L) {
    out <- c(out, sprintf("  %-*s  %s", w, "total", cells[[length(cells)]]))
  }
  out
}

# one bullet per pack (no alignment)
mc_manifest_bullets <- function(labels, sizes) {
  sizes <- as.numeric(sizes)
  items <- sprintf("%s (%s)", labels, trimws(format(fs::fs_bytes(sizes))))
  if (length(labels) > 1L) {
    items <- c(
      items,
      sprintf("total (%s)", trimws(format(fs::fs_bytes(sum(sizes)))))
    )
  }
  bullets(items)
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
          "Download failed for {.val {row[['group_id']]}}:
           {conditionMessage(e)}.",
          "i" = "URL: {.url {url}}"
        ),
        call = NULL
      )
    }
  )
  if (!identical(as.integer(status), 0L)) {
    cli::cli_abort(
      c(
        "Download failed for {.val {row[['group_id']]}}
         (status {status}).",
        "i" = "URL: {.url {url}}"
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
        "Can't download {length(rows)} clock-data pack{?s} without
         confirmation in a non-interactive session.",
        "i" = "Assets dir: {.path {dir}}",
        mc_manifest_bullets(ids, sizes),
        "i" = "Pass {.code ask = FALSE} to consent, or pre-stage the files
               and point {.arg ext_data} at them."
      ),
      call = NULL
    )
  }
  ok <- mc_ask_yes_no(
    header = cli::format_inline(
      "Download {length(rows)} clock-data pack{?s}:"
    ),
    labels = ids,
    sizes = sizes,
    dir = dir,
    question = sprintf("Download %d pack(s)?", length(rows))
  )
  if (!ok) {
    cli::cli_abort(
      "Download declined for {.val {ids}}.",
      call = NULL
    )
  }
  invisible(TRUE)
}

# pre-fetch packs into the assets dir
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
    "{.arg ext_data} must be {.code NULL}, an assets-dir path, a loaded pack,
     or a list of loaded packs.",
    call = NULL
  )
}

# load packs for needed groups into memory
#' @export
load_mc_assets <- function(groups, ext_data = NULL, ask = TRUE) {
  checkmate::assert_flag(ask)
  groups <- unique(as.character(groups))
  groups <- groups[nzchar(groups)]
  if (identical(groups, "all")) {
    groups <- mc_external_groups()
  }
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
            "Need pack {.val {g}}, but it is not in {.arg ext_data}.",
            "i" = "Closed set -- no download. Include it in {.arg ext_data}
                   or pass an assets dir / {.code NULL}."
          ),
          call = NULL
        )
      }
      pack
    })
    extra <- setdiff(names(canon), groups)
    if (length(extra)) {
      cli::cli_warn(
        "Ignoring unused pack{?s} in {.arg ext_data}: {.val {extra}}.",
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
          "Pack{cli::qty(gone)}{?s} {.val {gone}} not found in
           {.path {dir}}.",
          "i" = "Closed set -- no download. Stage the file{cli::qty(gone)}{?s}
                 there, or pass {.code ext_data = NULL} to allow download."
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
    parts <- c(parts, cli::format_inline("{n_downloaded} downloaded pack{?s}"))
  }
  if (n_stale) {
    parts <- c(parts, cli::format_inline("{n_stale} superseded pack{?s}"))
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
        "Can't delete {what} without confirmation in a non-interactive
         session.",
        "i" = "Assets dir: {.path {dir}}",
        mc_manifest_bullets(labels, sizes),
        "i" = "Pass {.code ask = FALSE} to consent."
      ),
      call = NULL
    )
  }
  mc_ask_yes_no(
    header = paste0("Delete ", what, ":"),
    labels = labels,
    sizes = sizes,
    dir = dir,
    question = sprintf("Delete %d pack(s)?", length(files))
  )
}

# remove every pack this package put in the assets dir
#' @export
clear_mc_assets <- function(groups = "all", ask = TRUE) {
  checkmate::assert_flag(ask)
  dir <- mc_resolve_assets_dir()
  downloaded <- mc_staged_files(groups)
  stale <- mc_stale_files(groups)
  files <- c(downloaded, stats::setNames(stale, mc_stale_labels(stale)))
  if (!length(files)) {
    cli::cli_inform("No clock assets to clear in {.path {dir}}.")
    return(invisible(character(0)))
  }
  if (!mc_consent_delete(files, dir, ask, n_stale = length(stale))) {
    cli::cli_inform("Deletion declined -- nothing removed from {.path {dir}}.")
    return(invisible(character(0)))
  }
  freed <- sum(as.numeric(fs::file_size(files)))
  # fs::file_delete() throws on failure
  fs::file_delete(files)
  cli::cli_inform(
    "Removed {mc_delete_summary(length(downloaded), length(stale))}
     ({format(fs::fs_bytes(freed))}) from {.path {dir}}."
  )
  invisible(files)
}
