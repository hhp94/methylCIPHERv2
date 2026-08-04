# vendor methylCIPHER-meta into the package. source("data-raw/sync.R") then sync()

# setup

for (pkg in c(
  "jsonlite",
  "qs2",
  "usethis",
  "digest",
  "processx",
  "fs",
  "rlang"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package '", pkg, "'. Install it first.", call. = FALSE)
  }
}

`%||%` <- rlang::`%||%`

SYNC_SCRIPT <- file.path("data-raw", "sync.R")
ACCESSORS_FILE <- file.path("R", "accessors.R")

# shared constants and stack-operand rules from R/accessors.R.
# sourced into its own env (definitions only).
if (!file.exists(ACCESSORS_FILE)) {
  stop(
    "sync.R runs from the package root; ",
    ACCESSORS_FILE,
    " not found",
    call. = FALSE
  )
}
mc_runtime <- new.env(parent = environment())
source(ACCESSORS_FILE, local = mc_runtime)

asset_dir <- file.path("data-raw", "assets")
meta_dir <- file.path("data-raw", "methylCIPHER-meta")

# gitignored lockfile for external pack rebuild skip.
LOCKFILE <- file.path(asset_dir, "lockfile.rds")

META_REMOTE <- "https://github.com/hhp94/methylCIPHER-meta.git"

# external families as release assets, rest in sysdata
EXTERNAL_GROUPS <- c("SystemsAge", "PCClocks", "PCBrainAge")

# single external clocks inside an otherwise-bundled group (group in both buckets)
EXTERNAL_CLOCKS <- c("Zhang2019BLUP")

# bump when pack layout changes (new payload_hash).
EXTERNAL_ENCODING_VERSION <- 3L

# pack staleness keys on encoder code plus upstream bundle_hashes. re-source after editing.
sync_code_fingerprint <- function() {
  files <- c(SYNC_SCRIPT, ACCESSORS_FILE)
  if (!all(file.exists(files))) {
    # code we cannot read is never a cache hit
    return(NA_character_)
  }
  digest::digest(c(
    list(EXTERNAL_ENCODING_VERSION),
    lapply(files, parse, keep.source = FALSE)
  ))
}

# pin-only fields, excluded from the content-addressed pack
EXTERNAL_PIN_FIELDS <- c("source_git_sha", "manifest_generated_at_sha")

# field registries (meta JSON -> shipped catalog)

FIELD_REGISTRY <- c(
  "clock_id",
  "group_id",
  "weights_format",
  "computation_type",
  "pmid",
  "output_transform",
  "normalization",
  "imputation",
  "intercept",
  "covariates",
  "recipe",
  "components",
  "shared",
  "probe_sets",
  "code_ref",
  "code_deps",
  "definition",
  "n_cpgs",
  "n_cpgs_normalization",
  "license"
)

# build-time only, stripped after resolution. `n_cpgs` cross-checks the derived panel.
# `shared` stays: recipe operands resolve name -> file through it.
CATALOG_BUILD_ONLY_FIELDS <- c("file_refs", "n_cpgs")

trim_build_only_fields <- function(clocks) {
  lapply(clocks, function(e) {
    e[CATALOG_BUILD_ONLY_FIELDS] <- NULL
    e
  })
}

BIB_INST_PATH <- file.path("inst", "bibliography", "clocks.bib")

GROUP_FIELD_REGISTRY <- c("group_id", "members", "shared_tensors", "routing")

IMPUTATION_FIELDS <- c("policy", "ref")
COMPONENT_FIELDS <- c(
  "name",
  "file",
  "row_key",
  "col_key",
  "intercept",
  "covariates"
)
SHARED_FIELDS <- c("name", "file")
PROBE_SET_FIELDS <- c("name", "role", "file", "n")
RECIPE_STEP_DROP <- c("note")

# one fixtures[] entry per cohort (weights_extraction.md sec 7)
FIXTURE_FIELDS <- c(
  "cohort",
  "expected",
  "missing",
  "oracle",
  "server_col",
  "server_raw",
  "server_normalization"
)

# keep named fields, preserve explicit JSON nulls
keep_fields <- function(x, fields) {
  if (is.null(x) || !is.list(x)) {
    return(NULL)
  }
  out <- list()
  for (f in fields) {
    if (f %in% names(x)) {
      out[f] <- list(x[[f]])
    }
  }
  if (!length(out)) NULL else out
}

keep_fields_each <- function(xs, fields) {
  if (is.null(xs) || !is.list(xs)) {
    return(NULL)
  }
  lapply(xs, function(item) keep_fields(item, fields))
}

prune_recipe <- function(recipe) {
  if (is.null(recipe)) {
    return(NULL)
  }
  lapply(recipe, function(step) {
    if (!is.list(step)) {
      return(step)
    }
    drop <- intersect(RECIPE_STEP_DROP, names(step))
    if (length(drop)) {
      step[drop] <- NULL
    }
    step
  })
}

# upstream `fixtures` is a list, one entry per registry cohort
prune_fixtures <- function(fxs) {
  if (is.null(fxs) || !is.list(fxs) || !length(fxs)) {
    return(NULL)
  }
  out <- lapply(fxs, function(fx) keep_fields(fx, FIXTURE_FIELDS))
  out <- Filter(Negate(is.null), out)
  if (!length(out)) NULL else out
}

# declared probe_sets kept for every weights_format (DunedinPACE QN background)
prune_clock_meta <- function(meta) {
  out <- keep_fields(meta, FIELD_REGISTRY)
  if (is.null(out)) {
    return(list())
  }
  if ("imputation" %in% names(out)) {
    out[["imputation"]] <- keep_fields(out[["imputation"]], IMPUTATION_FIELDS)
  }
  if ("components" %in% names(out)) {
    out[["components"]] <- keep_fields_each(
      out[["components"]],
      COMPONENT_FIELDS
    )
  }
  if ("shared" %in% names(out)) {
    out[["shared"]] <- keep_fields_each(out[["shared"]], SHARED_FIELDS)
  }
  if ("probe_sets" %in% names(out)) {
    out[["probe_sets"]] <- keep_fields_each(
      out[["probe_sets"]],
      PROBE_SET_FIELDS
    )
  }
  if ("recipe" %in% names(out)) {
    out[["recipe"]] <- prune_recipe(out[["recipe"]])
  }
  out
}

prune_group_meta <- function(gmeta) {
  keep_fields(gmeta, GROUP_FIELD_REGISTRY) %||% list()
}

# bibliography

CITATION_FIELDS <- c("clock_id", "pmid", "role", "bib_key")
CITATION_ROLES <- c("primary", "cite_also")

# paper fields from clocks.bib onto the citation join. pmid stays on the csv.
BIB_FIELDS <- c(
  "title",
  "author",
  "year",
  "journal",
  "volume",
  "number",
  "pages",
  "doi",
  "url"
)

# clock -> paper join (1:N). meta pmid is the primary only
read_clock_citations <- function(repo_path) {
  path <- file.path(repo_path, "bibliography", "clock_citations.csv")
  df <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  for (col in CITATION_FIELDS) {
    if (!col %in% names(df)) {
      stop("clock_citations.csv is missing column '", col, "'", call. = FALSE)
    }
    df[[col]] <- trimws(df[[col]])
  }
  bad <- df[!nzchar(df[["bib_key"]]), , drop = FALSE]
  if (nrow(bad)) {
    stop(
      "clock_citations.csv has ",
      nrow(bad),
      " row(s) with an empty bib_key (upstream gap): ",
      paste(unique(bad[["clock_id"]]), collapse = ", "),
      call. = FALSE
    )
  }
  bad_role <- setdiff(unique(df[["role"]]), CITATION_ROLES)
  if (length(bad_role)) {
    stop(
      "clock_citations.csv has unknown role(s): ",
      paste(bad_role, collapse = ", "),
      call. = FALSE
    )
  }
  df[, CITATION_FIELDS, drop = FALSE]
}

# citation rows for released clocks, primary first (both roles kept)
build_citations_table <- function(citations, clock_ids) {
  df <- citations[citations[["clock_id"]] %in% clock_ids, , drop = FALSE]
  missing_cites <- setdiff(clock_ids, unique(df[["clock_id"]]))
  if (length(missing_cites)) {
    stop(
      "clock(s) in manifest.json with no clock_citations.csv row: ",
      paste(missing_cites, collapse = ", "),
      call. = FALSE
    )
  }
  n_primary <- tapply(
    df[["role"]] == "primary",
    df[["clock_id"]],
    sum
  )
  if (any(n_primary != 1L)) {
    stop(
      "clock(s) without exactly one primary citation: ",
      paste(names(n_primary)[n_primary != 1L], collapse = ", "),
      call. = FALSE
    )
  }
  df <- df[
    order(df[["clock_id"]], df[["role"]] != "primary", df[["bib_key"]]),
    ,
    drop = FALSE
  ]
  row.names(df) <- NULL
  df
}

# parse clocks.bib as upstream emits it (one field per line, fixed order).
BIB_INDENT <- "  "
BIB_ASSIGN <- " = {"

# present on every entry, whatever its type (lib_bib.py REQUIRED_FIELDS)
BIB_REQUIRED <- c("title", "author", "year", "pmid")

# brace-protected casing ({DNA}, {eLife}) is a BibTeX rendering hint, not data
bib_unbrace <- function(x) {
  gsub("}", "", gsub("{", "", x, fixed = TRUE), fixed = TRUE)
}

# "@article{Key_2022_35029144," -> "Key_2022_35029144". entry type is not checked
bib_entry_key <- function(line, where) {
  parts <- strsplit(line, "{", fixed = TRUE)[[1L]]
  if (
    length(parts) != 2L ||
      !startsWith(parts[[1L]], "@") ||
      !endsWith(parts[[2L]], ",")
  ) {
    stop(where, ": unreadable entry header: ", line, call. = FALSE)
  }
  key <- trimws(substr(parts[[2L]], 1L, nchar(parts[[2L]]) - 1L))
  if (!nzchar(key)) {
    stop(where, ": entry header declares no key: ", line, call. = FALSE)
  }
  key
}

# the field lines of one entry -> named character vector, unbraced
bib_entry_fields <- function(body, where) {
  out <- character()
  for (line in body) {
    if (!startsWith(line, BIB_INDENT)) {
      stop(where, ": field line is not indented: ", line, call. = FALSE)
    }
    txt <- substring(line, nchar(BIB_INDENT) + 1L)
    if (startsWith(txt, " ")) {
      stop(where, ": field line is over-indented: ", line, call. = FALSE)
    }
    # trailing comma on every field but the last
    if (endsWith(txt, "},")) {
      txt <- substr(txt, 1L, nchar(txt) - 1L)
    }
    if (!endsWith(txt, "}")) {
      stop(
        where,
        ": field line does not close its value: ",
        line,
        call. = FALSE
      )
    }
    at <- regexpr(BIB_ASSIGN, txt, fixed = TRUE)
    if (at < 1L) {
      stop(where, ": field line has no ' = {': ", line, call. = FALSE)
    }
    name <- substr(txt, 1L, at - 1L)
    value <- substr(txt, at + nchar(BIB_ASSIGN), nchar(txt) - 1L)
    if (name %in% names(out)) {
      stop(where, ": field '", name, "' is declared twice", call. = FALSE)
    }
    out[[name]] <- trimws(bib_unbrace(value))
  }
  missing <- setdiff(BIB_REQUIRED, names(out))
  if (length(missing)) {
    stop(
      where,
      ": entry is missing required field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out
}

# clocks.bib -> one row per entry, columns BIB_FIELDS + pmid + bib_key
read_bib_fields <- function(repo_path) {
  path <- file.path(repo_path, "bibliography", "clocks.bib")
  con <- file(path, encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  starts <- which(startsWith(lines, "@"))
  if (!length(starts)) {
    stop("clocks.bib has no entries", call. = FALSE)
  }
  ends <- c(starts[-1L] - 1L, length(lines))

  wanted <- c(BIB_FIELDS, "pmid")
  cells <- matrix(
    NA_character_,
    nrow = length(starts),
    ncol = length(wanted),
    dimnames = list(NULL, wanted)
  )
  keys <- character(length(starts))
  seen <- character()

  for (i in seq_along(starts)) {
    block <- lines[seq.int(starts[i], ends[i])]
    # the blank separator line belongs to no entry
    block <- block[seq_len(max(which(nzchar(trimws(block)))))]
    where <- paste0("clocks.bib entry at line ", starts[i])
    if (length(block) < 3L || !identical(block[[length(block)]], "}")) {
      stop(where, ": entry does not close with a lone '}'", call. = FALSE)
    }
    keys[[i]] <- bib_entry_key(block[[1L]], where)
    fields <- bib_entry_fields(block[-c(1L, length(block))], where)
    seen <- c(seen, names(fields))
    have <- intersect(wanted, names(fields))
    cells[i, have] <- fields[have]
  }

  # every field we read must appear on some entry.
  never <- setdiff(wanted, unique(seen))
  if (length(never)) {
    stop(
      "clocks.bib declares no ",
      paste(never, collapse = ", "),
      " field on any of its ",
      length(starts),
      " entries -- the emitter's field names moved",
      call. = FALSE
    )
  }

  dup <- unique(keys[duplicated(keys)])
  if (length(dup)) {
    stop(
      "clocks.bib has duplicate entry key(s): ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }

  out <- as.data.frame(cells, stringsAsFactors = FALSE)
  out[["bib_key"]] <- keys
  row.names(out) <- NULL
  out
}

# join paper fields onto citations. every key and both pmid copies must agree.
attach_bib_fields <- function(citations, bib) {
  absent <- setdiff(unique(citations[["bib_key"]]), bib[["bib_key"]])
  if (length(absent)) {
    stop(
      "clock_citations.csv names bib_key(s) missing from clocks.bib: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  idx <- match(citations[["bib_key"]], bib[["bib_key"]])

  csv_pmid <- trimws(citations[["pmid"]])
  bib_pmid <- trimws(bib[["pmid"]][idx])
  clash <- nzchar(csv_pmid) & !is.na(bib_pmid) & csv_pmid != bib_pmid
  if (any(clash)) {
    stop(
      "pmid disagrees between clock_citations.csv and clocks.bib for: ",
      paste(unique(citations[["bib_key"]][clash]), collapse = ", "),
      call. = FALSE
    )
  }

  for (f in BIB_FIELDS) {
    citations[[f]] <- bib[[f]][idx]
  }
  citations
}

vendor_bibliography <- function(repo_path) {
  src <- file.path(repo_path, "bibliography", "clocks.bib")
  fs::dir_create(dirname(BIB_INST_PATH))
  fs::file_copy(src, BIB_INST_PATH, overwrite = TRUE)
  sz <- file.info(BIB_INST_PATH)$size
  message(sprintf("sync: wrote %s (%.1f KB)", BIB_INST_PATH, sz / 1024))
  invisible(list(path = BIB_INST_PATH, size_bytes = sz))
}

# resolve SoT

git_exec <- function(..., dir = NULL) {
  args <- c(...)
  if (!is.null(dir)) {
    args <- c("-C", dir, args)
  }
  res <- processx::run("git", args, error_on_status = FALSE)
  if (res$status != 0L) {
    stop(
      "git ",
      paste(args, collapse = " "),
      " failed:\n",
      res$stderr,
      call. = FALSE
    )
  }
  if (nzchar(res$stdout)) strsplit(res$stdout, "\r?\n")[[1]] else character()
}

# git_exec when stdout is a single value.
git_value <- function(..., dir = NULL) {
  out <- trimws(git_exec(..., dir = dir))
  out <- out[nzchar(out)]
  if (length(out) != 1L) {
    stop(
      "git returned ",
      length(out),
      " lines where one value was expected:\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  out
}

# clone/fetch meta_dir and checkout a commit (NULL = origin tip).
resolve_source <- function(source_git_sha = NULL) {
  fs::dir_create(dirname(meta_dir))

  # discard a non-repo meta_dir (interrupted clone) and re-clone.
  is_repo <- dir.exists(file.path(meta_dir, ".git"))
  if (dir.exists(meta_dir) && !is_repo) {
    unlink(meta_dir, recursive = TRUE, force = TRUE)
  }

  if (!is_repo) {
    # clone to a temp sibling, then rename so a partial clone cannot poison meta_dir.
    message("sync: cloning ", META_REMOTE, " -> ", meta_dir)
    tmp <- paste0(meta_dir, ".tmp-", Sys.getpid())
    unlink(tmp, recursive = TRUE, force = TRUE)
    git_exec("clone", "--filter=blob:none", META_REMOTE, tmp)
    if (!file.rename(tmp, meta_dir)) {
      stop(
        "could not move the fresh clone ",
        tmp,
        " -> ",
        meta_dir,
        call. = FALSE
      )
    }
  } else {
    message("sync: fetching ", META_REMOTE, " into ", meta_dir)
    git_exec("remote", "set-url", "origin", META_REMOTE, dir = meta_dir)
    git_exec("fetch", "origin", "--tags", "--prune", dir = meta_dir)
  }

  if (!is.null(source_git_sha) && nzchar(source_git_sha)) {
    ref <- source_git_sha
  } else {
    # resolve origin/HEAD, or stop -- no branch-name fallback
    ref <- git_value(
      "symbolic-ref",
      "--short",
      "refs/remotes/origin/HEAD",
      dir = meta_dir
    )
  }
  message("sync: checkout ", ref)
  # force detach checkout, mirror is disposable
  git_exec("checkout", "-f", "--detach", ref, dir = meta_dir)

  path <- as.character(fs::path_real(meta_dir))
  sha <- git_value("rev-parse", "HEAD", dir = path)
  list(path = path, source_git_sha = sha)
}

# manifest

read_manifest <- function(repo_path) {
  jsonlite::fromJSON(
    file.path(repo_path, "manifest.json"),
    simplifyVector = FALSE
  )
}

# catalog

# list meta files under weights/, clock vs group by basename
list_meta_files <- function(repo_path) {
  metas <- list.files(
    file.path(repo_path, "weights"),
    pattern = "\\.meta\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  basenames <- basename(metas)
  list(
    clock = metas[basenames != "_group.meta.json"],
    group = metas[basenames == "_group.meta.json"]
  )
}

# manifest.json: one row per weights_status=done clock, with bundle_hash
manifest_clocks <- function(manifest) {
  rows <- manifest[["clocks"]] %||% list()
  if (!length(rows)) {
    stop("manifest.json declares no clocks", call. = FALSE)
  }
  out <- list()
  for (r in rows) {
    cid <- as.character(r[["clock_id"]] %||% NA_character_)
    if (is.na(cid) || !nzchar(cid)) {
      stop("manifest.json has a clock row without a clock_id", call. = FALSE)
    }
    # bundle_hash is the pack staleness key (missing -> rebuild)
    bh <- as.character(r[["bundle_hash"]] %||% NA_character_)
    if (is.na(bh) || !nzchar(bh)) {
      stop(
        "manifest.json row for '",
        cid,
        "' has no bundle_hash",
        call. = FALSE
      )
    }
    out[[cid]] <- list(
      bundle_hash = bh,
      out_sha256 = as.character(r[["out_sha256"]] %||% NA_character_),
      verification_status = as.character(
        r[["verification_status"]] %||% NA_character_
      )
    )
  }
  out
}

# covariates from recipe steps plus top-level mirror, plus Female for sex routing
extract_covariates <- function(meta) {
  cid <- as.character(meta[["clock_id"]] %||% NA_character_)
  covs <- covariate_names(meta[["covariates"]], paste0("clock '", cid, "'"))
  for (step in meta[["recipe"]] %||% list()) {
    covs <- c(
      covs,
      covariate_names(
        step[["covariates"]],
        paste0("clock '", cid, "' recipe op '", step[["op"]] %||% "?", "'")
      )
    )
  }
  unique(covs[nzchar(covs) & !is.na(covs)])
}

# one declared covariate list -> registry names (object = coefs, array = names)
covariate_names <- function(x, where) {
  if (is.null(x)) {
    return(character())
  }
  if (is.character(x)) {
    return(x)
  }
  if (!is.list(x)) {
    stop(
      where,
      ": `covariates` is neither an object nor an array",
      call. = FALSE
    )
  }
  if (!length(x)) {
    return(character())
  }
  nms <- names(x)
  if (is.null(nms)) {
    return(as.character(unlist(x, use.names = FALSE)))
  }
  if (!all(nzchar(nms))) {
    stop(where, ": `covariates` is a partially named object", call. = FALSE)
  }
  nms
}

# sample-axis classification (cross_sample vs per_sample ops)
CROSS_SAMPLE_OPS <- c("cohort_zscore")

# closed op vocabulary -- unknown ops stop the sync
KNOWN_OPS <- c(
  # per-sample: each row is transformed against frozen parameters only
  "linear",
  "linear_mean",
  "impute",
  "stack",
  "transform",
  "poly",
  "row_sum",
  "fitage_kdm",
  "epitoc2",
  "project",
  "center_scale",
  "sample_scale",
  # cross-sample (CROSS_SAMPLE_OPS above)
  "cohort_zscore"
)

# 1-based position of the first cross_sample recipe step, or NA
first_cross_sample_step <- function(meta) {
  cid <- as.character(meta[["clock_id"]] %||% NA_character_)
  ops <- vapply(
    meta[["recipe"]] %||% list(),
    function(s) as.character(s[["op"]] %||% NA_character_),
    character(1L)
  )
  unknown <- unique(ops[is.na(ops) | !(ops %in% KNOWN_OPS)])
  if (length(unknown)) {
    stop(
      "clock '",
      cid,
      "': recipe declares op(s) outside the known vocabulary: ",
      paste(unknown, collapse = ", "),
      ". Classify each as per-sample or cross-sample and add it to KNOWN_OPS ",
      "(and CROSS_SAMPLE_OPS if it is batch-derived).",
      call. = FALSE
    )
  }
  hits <- which(ops %in% CROSS_SAMPLE_OPS)
  if (!length(hits)) NA_integer_ else as.integer(hits[[1L]])
}

# other clocks this clock consumes, by clock_id (recipe `inputs`)
extract_clock_inputs <- function(meta) {
  ins <- character()
  for (step in meta[["recipe"]] %||% list()) {
    ins <- c(ins, as.character(unlist(step[["inputs"]] %||% character())))
  }
  unique(ins[nzchar(ins) & !is.na(ins)])
}

# relative path from repo root (forward slashes).
rel_from_repo <- function(abs_path, repo_path) {
  as.character(fs::path_rel(abs_path, start = repo_path))
}

# repo-relative weights/ path -> absolute.
resolve_repo_rel <- function(repo_path, rel) {
  if (is.null(rel) || !nzchar(rel)) {
    return(NA_character_)
  }
  file.path(repo_path, rel)
}

build_catalog <- function(repo_path, manifest) {
  files <- list_meta_files(repo_path)

  groups <- list()
  for (gp in files$group) {
    gmeta <- jsonlite::fromJSON(gp, simplifyVector = FALSE)
    gid <- as.character(gmeta[["group_id"]] %||% NA_character_)
    if (is.na(gid) || !nzchar(gid)) {
      stop("group meta declares no group_id: ", gp, call. = FALSE)
    }
    if (!is.null(groups[[gid]])) {
      stop("two group metas declare group_id '", gid, "'", call. = FALSE)
    }
    entry <- prune_group_meta(gmeta)
    if (!is.null(entry[["members"]])) {
      entry[["members"]] <- unlist(entry[["members"]], use.names = FALSE)
    }
    if (!is.null(entry[["shared_tensors"]])) {
      entry[["shared_tensors"]] <- unlist(
        entry[["shared_tensors"]],
        use.names = FALSE
      )
    }
    entry[["path"]] <- rel_from_repo(gp, repo_path)
    groups[[gid]] <- entry
  }

  # clock_id -> meta path, selected by the manifest
  by_id <- list()
  for (mp in files$clock) {
    cid <- as.character(
      jsonlite::fromJSON(mp, simplifyVector = FALSE)[["clock_id"]] %||%
        NA_character_
    )
    if (is.na(cid) || !nzchar(cid)) {
      stop("clock meta declares no clock_id: ", mp, call. = FALSE)
    }
    if (!is.null(by_id[[cid]])) {
      stop(
        "two clock metas declare clock_id '",
        cid,
        "': ",
        by_id[[cid]],
        " and ",
        mp,
        call. = FALSE
      )
    }
    by_id[[cid]] <- mp
  }
  released <- manifest_clocks(manifest)
  absent <- setdiff(names(released), names(by_id))
  if (length(absent)) {
    stop(
      "manifest.json names clock(s) with no meta on disk: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }

  citations <- attach_bib_fields(
    build_citations_table(
      read_clock_citations(repo_path),
      names(released)
    ),
    read_bib_fields(repo_path)
  )

  clocks <- list()
  for (cid in names(released)) {
    mp <- by_id[[cid]]
    meta <- jsonlite::fromJSON(mp, simplifyVector = FALSE)
    gid <- as.character(meta[["group_id"]] %||% NA_character_)
    wf <- as.character(meta[["weights_format"]] %||% NA_character_)
    if (is.na(gid) || !nzchar(gid)) {
      stop("clock '", cid, "' declares no group_id (", mp, ")", call. = FALSE)
    }
    if (is.na(wf) || !nzchar(wf)) {
      stop(
        "clock '",
        cid,
        "' declares no weights_format (",
        mp,
        ")",
        call. = FALSE
      )
    }

    entry <- prune_clock_meta(meta)

    if (!is.null(entry[["normalization"]])) {
      entry[["normalization"]] <- unlist(
        entry[["normalization"]],
        use.names = FALSE
      )
    }
    entry[["output_transform"]] <- as.character(
      entry[["output_transform"]] %||% "identity"
    )
    if (!is.null(entry[["computation_type"]])) {
      entry[["computation_type"]] <- as.character(entry[["computation_type"]])
    }
    # provenance pointer to the primary paper, not the citation set (mc_citations)
    entry[["pmid"]] <- as.character(entry[["pmid"]] %||% NA_character_)

    entry[["imputation_policy"]] <- as.character(
      entry[["imputation"]][["policy"]] %||% NA_character_
    )
    entry[["covariates_required"]] <- extract_covariates(meta)
    entry[["clock_inputs"]] <- extract_clock_inputs(meta)
    entry[["cross_sample_at"]] <- first_cross_sample_step(meta)
    # the field is per-clock. a group may be partly bundled and partly external
    entry[["external_group"]] <- cid %in%
      EXTERNAL_CLOCKS ||
      gid %in% EXTERNAL_GROUPS
    entry[["fixtures"]] <- prune_fixtures(meta[["fixtures"]])
    entry[["meta_path"]] <- rel_from_repo(mp, repo_path)
    entry[["bundle_hash"]] <- released[[cid]][["bundle_hash"]]
    entry[["verification_status"]] <- released[[cid]][["verification_status"]]
    entry[["file_refs"]] <- declared_tensors(meta)
    entry[["coef_path"]] <- if (identical(wf, "cpg_coefficient")) {
      coef_path(cid, gid)
    } else {
      NULL
    }

    clocks[[cid]] <- entry
  }

  list(
    clocks = clocks,
    groups = groups,
    citations = citations,
    source_git_sha = NA_character_,
    schema_version = manifest[["schema_version"]],
    n_clocks = length(clocks)
  )
}

# declared paths + declared tensor shapes

# declared paths: weights/ is vendored, papers/ is recognised and skipped
PATH_VENDOR_PREFIX <- "weights/"
PATH_POLICY_SKIP_PREFIX <- "papers/"

# normalize + validate one declared pointer (path and whether we vendor it)
declared_path <- function(x, field, cid) {
  p <- as.character(x %||% NA_character_)
  if (length(p) != 1L || is.na(p) || !nzchar(trimws(p))) {
    stop(
      "clock '",
      cid,
      "': declared `",
      field,
      "` is empty or not a single path",
      call. = FALSE
    )
  }
  p <- as.character(fs::path_norm(gsub("\\", "/", trimws(p), fixed = TRUE)))
  if (fs::is_absolute_path(p) || startsWith(p, "..")) {
    stop(
      "clock '",
      cid,
      "': declared `",
      field,
      "` is not repo-relative: ",
      p,
      call. = FALSE
    )
  }
  if (startsWith(p, PATH_VENDOR_PREFIX)) {
    return(list(path = p, vendor = TRUE))
  }
  if (startsWith(p, PATH_POLICY_SKIP_PREFIX)) {
    return(list(path = p, vendor = FALSE))
  }
  stop(
    "clock '",
    cid,
    "': declared `",
    field,
    "` points under neither ",
    PATH_VENDOR_PREFIX,
    " nor ",
    PATH_POLICY_SKIP_PREFIX,
    ": ",
    p,
    call. = FALSE
  )
}

# the vendored path, or stop. For pointers whose prefix is fixed by contract.
vendored_path <- function(x, field, cid) {
  d <- declared_path(x, field, cid)
  if (!d$vendor) {
    stop(
      "clock '",
      cid,
      "': declared `",
      field,
      "` must be under ",
      PATH_VENDOR_PREFIX,
      ": ",
      d$path,
      call. = FALSE
    )
  }
  d$path
}

# weights/{group_id}/{clock_id}.csv.gz -- derived coef path, never globbed
coef_path <- function(clock_id, group_id) {
  cid <- trimws(as.character(clock_id %||% NA_character_))
  gid <- trimws(as.character(group_id %||% NA_character_))
  if (is.na(cid) || !nzchar(cid) || is.na(gid) || !nzchar(gid)) {
    stop(
      "coef_path needs both clock_id and group_id; got '",
      cid,
      "' / '",
      gid,
      "'",
      call. = FALSE
    )
  }
  paste0(PATH_VENDOR_PREFIX, gid, "/", cid, ".csv.gz")
}

# path kind from the declaring field: tensor vs R source
TENSOR_KIND <- "tensor"
CODE_KIND <- "r_source"

# declared tensor shape: row_key is col 1, col_key the rest
tensor_spec <- function(
  kind = TENSOR_KIND,
  row_key = NULL,
  col_key = NULL,
  field = NA_character_
) {
  cols <- character()
  if (!is.null(col_key) && nzchar(as.character(col_key))) {
    cols <- trimws(strsplit(as.character(col_key), ",", fixed = TRUE)[[1L]])
  }
  list(
    kind = kind,
    field = as.character(field),
    row_key = if (is.null(row_key)) NA_character_ else as.character(row_key),
    col_key = cols
  )
}

# every vendored path a meta declares, as path -> tensor_spec
declared_tensors <- function(meta) {
  cid <- as.character(meta[["clock_id"]] %||% NA_character_)
  out <- list()
  add <- function(
    x,
    field,
    kind = TENSOR_KIND,
    row_key = NULL,
    col_key = NULL
  ) {
    d <- declared_path(x, field, cid)
    if (!d$vendor) {
      return(invisible(NULL))
    }
    out[[d$path]] <<- tensor_spec(kind, row_key, col_key, field)
    invisible(NULL)
  }

  if (identical(as.character(meta[["weights_format"]]), "cpg_coefficient")) {
    add(
      coef_path(meta[["clock_id"]], meta[["group_id"]]),
      "coef_path",
      row_key = "cpg",
      col_key = "coefficient"
    )
  }
  # vendor_mean ref is a scalar path or absent (empty is an error)
  if (!is.null(meta[["imputation"]][["ref"]])) {
    add(
      meta[["imputation"]][["ref"]],
      "imputation.ref",
      row_key = "cpg",
      col_key = "value"
    )
  }
  for (comp in meta[["components"]] %||% list()) {
    add(
      comp[["file"]],
      "components[].file",
      row_key = comp[["row_key"]],
      col_key = comp[["col_key"]]
    )
  }
  # probe_sets and shared: cpg list or cpg,value (header says which)
  for (ps in meta[["probe_sets"]] %||% list()) {
    add(ps[["file"]], "probe_sets[].file")
  }
  for (sh in meta[["shared"]] %||% list()) {
    add(sh[["file"]], "shared[].file")
  }
  if (!is.null(meta[["code_ref"]])) {
    add(meta[["code_ref"]], "code_ref", kind = CODE_KIND)
  }
  for (dep in meta[["code_deps"]] %||% list()) {
    add(dep, "code_deps[]", kind = CODE_KIND)
  }
  out
}

# tensor IO

# load a declared weights/ tensor and assert row_key / col_key vs header
read_tensor_csv <- function(path, spec = NULL) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  header <- names(df)

  # row_key names column 1, col_key asserted only when declared
  if (!is.null(spec)) {
    mismatch <- function(want) {
      stop(
        "Tensor ",
        path,
        " has header (",
        paste(header, collapse = ", "),
        ") but declares (",
        paste(want, collapse = ", "),
        ")",
        call. = FALSE
      )
    }
    row_key <- spec[["row_key"]]
    if (!is.na(row_key) && !identical(header[[1L]], row_key)) {
      mismatch(row_key)
    }
    col_key <- spec[["col_key"]]
    if (length(col_key) && !identical(header[-1L], col_key)) {
      mismatch(c(row_key, col_key))
    }
  }

  # empty declared tensor is an upstream/sync gap
  if (!nrow(df)) {
    stop("Tensor ", path, " declares no rows", call. = FALSE)
  }
  if (ncol(df) == 1L) {
    return(as.character(df[[1L]]))
  }
  key <- as.character(df[[1L]])
  if (anyDuplicated(key)) {
    stop(
      "Duplicate keys in ",
      path,
      " -- a named vector would silently collapse them",
      call. = FALSE
    )
  }
  if (ncol(df) == 2L) {
    if (!is.numeric(df[[2L]])) {
      stop(
        "Tensor ",
        path,
        " is two-column but its value column is not numeric",
        call. = FALSE
      )
    }
    return(stats::setNames(df[[2L]], key))
  }
  df
}

# sex-routed aliases

# mint sex-routed alias stems from `_group.meta.json` routing.sex
SEX_SUFFIX <- c(female = "_Female", male = "_Male")

# aliases use kind=sex_routed_alias, weights_format / computation_type stay NA
KIND_SEX_ROUTED_ALIAS <- "sex_routed_alias"

# {stem: {female, male}} for one group, empty when it does not route
group_sex_routes <- function(gside) {
  routing <- gside[["routing"]][["sex"]]
  if (is.null(routing)) {
    return(list())
  }
  stems <- list()
  for (sx in names(SEX_SUFFIX)) {
    for (id in as.character(unlist(routing[[sx]] %||% character()))) {
      stem <- sub(paste0(SEX_SUFFIX[[sx]], "$"), "", id)
      stems[[stem]][[sx]] <- id
    }
  }
  stems
}

# earliest cross_sample_at over an alias's members (min, or NA)
alias_cross_sample_at <- function(members) {
  at <- vapply(
    members,
    function(m) {
      v <- m[["cross_sample_at"]]
      if (is.null(v) || !length(v)) NA_integer_ else as.integer(v)[[1L]]
    },
    integer(1L)
  )
  if (all(is.na(at))) NA_integer_ else min(at, na.rm = TRUE)
}

# mint one alias clock per routed stem (no panel of its own)
attach_sex_routed_aliases <- function(catalog) {
  for (gid in names(catalog[["groups"]])) {
    routes <- group_sex_routes(catalog[["groups"]][[gid]])
    for (stem in names(routes)) {
      if (!is.null(catalog[["clocks"]][[stem]])) {
        stop(
          "sex-routed alias '",
          stem,
          "' collides with an existing clock_id.",
          call. = FALSE
        )
      }
      # both routing.sex members must be in the manifest, or skip / stop
      declared <- vapply(
        names(SEX_SUFFIX),
        function(sx) as.character(routes[[stem]][[sx]] %||% NA_character_),
        character(1L)
      )
      have <- !is.na(declared) &
        vapply(
          declared,
          function(id) !is.na(id) && !is.null(catalog[["clocks"]][[id]]),
          logical(1L)
        )
      if (!any(have)) {
        message(
          "sync: sex-routed stem '",
          stem,
          "' has no vendored members; no alias minted"
        )
        next
      }
      if (!all(have)) {
        stop(
          "sex-routed stem '",
          stem,
          "' is half-vendored: ",
          paste(names(SEX_SUFFIX)[!have], collapse = ", "),
          " missing from the catalog",
          call. = FALSE
        )
      }
      donor <- catalog[["clocks"]][[routes[[stem]]$female]]
      members <- catalog[["clocks"]][c(
        routes[[stem]]$female,
        routes[[stem]]$male
      )]
      catalog[["clocks"]][[stem]] <- list(
        clock_id = stem,
        group_id = gid,
        kind = KIND_SEX_ROUTED_ALIAS,
        weights_format = NA_character_,
        computation_type = NA_character_,
        output_transform = "identity",
        normalization = donor[["normalization"]],
        routing = routes[[stem]],
        # routing.sex is the alias dependency edge (both members scored first)
        clock_inputs = c(routes[[stem]]$female, routes[[stem]]$male),
        covariates_required = "Female",
        imputation_policy = donor[["imputation_policy"]],
        pmid = donor[["pmid"]],
        # alias cites through its donor
        donor_clock_id = donor[["clock_id"]],
        license = donor[["license"]],
        # cross_sample_at derived from members (never assume NA)
        cross_sample_at = alias_cross_sample_at(members),
        external_group = FALSE
      )
    }
  }
  catalog[["n_clocks"]] <- length(catalog[["clocks"]])
  catalog
}

# materialize

# group ids on one side of the bundled/external split (a group may be on both)
split_group_ids <- function(catalog, external) {
  gids <- vapply(
    catalog[["clocks"]],
    function(e) {
      if (isTRUE(e[["external_group"]]) == external) {
        as.character(e[["group_id"]])
      } else {
        NA_character_
      }
    },
    character(1L)
  )
  sort(unique(gids[!is.na(gids)]))
}

# per-group tensor payload for group_ids on one side of the external split.
build_group_bundles <- function(
  repo_path,
  catalog,
  group_ids,
  external = FALSE
) {
  bundles <- list()
  for (gid in group_ids) {
    member_ids <- names(catalog[["clocks"]])[
      vapply(
        catalog[["clocks"]],
        function(c) {
          identical(c[["group_id"]], gid) &&
            isTRUE(c[["external_group"]]) == external
        },
        logical(1L)
      )
    ]
    members <- catalog[["clocks"]][member_ids]

    # path -> declared shape, unioned over members plus group shared_tensors
    specs <- list()
    for (entry in members) {
      for (rel in names(entry[["file_refs"]])) {
        specs[[rel]] <- entry[["file_refs"]][[rel]]
      }
    }
    gside <- catalog[["groups"]][[gid]]
    for (rel in as.character(gside[["shared_tensors"]] %||% character())) {
      rel <- vendored_path(rel, "shared_tensors[]", gid)
      if (is.null(specs[[rel]])) {
        specs[[rel]] <- tensor_spec(field = "shared_tensors[]")
      }
    }

    tensors <- list()
    for (rel in names(specs)) {
      spec <- specs[[rel]]
      abs <- resolve_repo_rel(repo_path, rel)
      # missing declared path is a hard stop
      if (!file.exists(abs)) {
        stop(
          "Declared tensor missing from snapshot: ",
          rel,
          " (group ",
          gid,
          ", declared by ",
          spec[["field"]],
          ")",
          call. = FALSE
        )
      }
      # dispatch on the DECLARING FIELD, not the file extension.
      tensors[[rel]] <- switch(
        spec[["kind"]],
        # code_ref entry points plus their code_deps closure
        r_source = list(type = "r_source", text = readLines(abs, warn = FALSE)),
        tensor = read_tensor_csv(abs, spec),
        stop(
          "Declared path ",
          rel,
          " has unknown kind '",
          spec[["kind"]],
          "'",
          call. = FALSE
        )
      )
    }

    bundles[[gid]] <- list(
      group_id = gid,
      clocks = member_ids,
      tensors = tensors
    )
  }
  bundles
}

# materialize each clock's scoring CpGs as probe_sets {name, role, cpgs}

# row labels of one loaded tensor (column 1 is the key)
tensor_row_keys <- function(t, where) {
  if (is.null(t)) {
    stop(where, ": tensor is absent from the group bundle", call. = FALSE)
  }
  keys <- if (is.data.frame(t)) {
    as.character(t[[1L]])
  } else if (is.numeric(t) && !is.null(names(t))) {
    names(t)
  } else if (is.character(t) && is.null(dim(t))) {
    as.character(t)
  } else {
    stop(where, ": tensor shape carries no row keys", call. = FALSE)
  }
  if (!length(keys)) {
    stop(where, ": tensor has no row keys", call. = FALSE)
  }
  keys
}

# materialize one probe_set from a file pointer.
materialize_probe_set <- function(ps, tensors, cid = NA_character_) {
  where <- paste0(
    "probe_set '",
    ps[["name"]] %||% "?",
    "' (role ",
    ps[["role"]] %||% "?",
    ") of clock ",
    cid
  )
  rel <- vendored_path(ps[["file"]], "probe_sets[].file", cid)
  cpgs <- tensor_row_keys(tensors[[rel]], where)
  if (anyDuplicated(cpgs)) {
    stop(where, " contains duplicate CpG id(s): ", rel, call. = FALSE)
  }
  # keep file so accessors can fetch values, not just row keys
  list(
    name = ps[["name"]],
    role = ps[["role"]],
    file = rel,
    cpgs = cpgs
  )
}

# union of a clock's own cpg-keyed components.
own_component_cpgs <- function(entry, tensors, cid) {
  cpgs <- character()
  for (comp in entry[["components"]] %||% list()) {
    if (!identical(comp[["row_key"]], "cpg")) {
      next
    }
    rel <- vendored_path(comp[["file"]], "components[].file", cid)
    cpgs <- c(
      cpgs,
      tensor_row_keys(tensors[[rel]], paste0("component ", rel, " of ", cid))
    )
  }
  unique(cpgs[nzchar(cpgs) & !is.na(cpgs)])
}

add_scoring_probe_set <- function(entry, cpgs) {
  c(
    entry[["probe_sets"]] %||% list(),
    list(list(name = "scoring_derived", role = "scoring", cpgs = unique(cpgs)))
  )
}

# the resolved scoring panel of one clock (post resolve_group_scoring_probe_sets)
resolved_scoring_cpgs <- function(entry, cid) {
  hits <- Filter(
    function(ps) identical(as.character(ps[["role"]]), "scoring"),
    entry[["probe_sets"]] %||% list()
  )
  if (length(hits) != 1L) {
    stop(
      cid,
      ": expected 1 scoring probe_set, found ",
      length(hits),
      call. = FALSE
    )
  }
  as.character(hits[[1L]][["cpgs"]])
}

# the clock's own cpg-keyed tensors (coef file or row_key==cpg components)
own_scoring_cpgs <- function(entry, tensors, cid) {
  cp <- entry[["coef_path"]]
  if (!is.null(cp)) {
    return(tensor_row_keys(
      tensors[[cp]],
      paste0("coef tensor ", cp, " of ", cid)
    ))
  }
  own_component_cpgs(entry, tensors, cid)
}

# derived scoring panel must equal declared n_cpgs
assert_declared_n_cpgs <- function(entry, cpgs, cid) {
  n <- entry[["n_cpgs"]]
  if (is.null(n) || !length(n)) {
    stop("clock '", cid, "' declares no n_cpgs", call. = FALSE)
  }
  if (length(cpgs) != as.integer(n)) {
    stop(
      "clock '",
      cid,
      "': derived scoring panel has ",
      length(cpgs),
      " CpGs but the meta declares n_cpgs=",
      as.integer(n),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# materialize probe_sets, then resolve each clock's scoring panel via own tensors / inputs DAG
resolve_group_scoring_probe_sets <- function(catalog, bundles) {
  for (gid in names(bundles)) {
    tensors <- bundles[[gid]][["tensors"]]
    ids <- bundles[[gid]][["clocks"]]

    for (cid in ids) {
      declared <- catalog[["clocks"]][[cid]][["probe_sets"]]
      if (length(declared)) {
        catalog[["clocks"]][[cid]][["probe_sets"]] <- lapply(
          declared,
          materialize_probe_set,
          tensors = tensors,
          cid = cid
        )
      }
    }

    memo <- new.env(parent = emptyenv())
    panel_of <- function(cid, stack) {
      hit <- memo[[cid]]
      if (!is.null(hit)) {
        return(hit)
      }
      if (cid %in% stack) {
        stop(
          "Dependency cycle among clocks: ",
          paste(c(stack, cid), collapse = " -> "),
          call. = FALSE
        )
      }
      entry <- catalog[["clocks"]][[cid]]
      cpgs <- own_scoring_cpgs(entry, tensors, cid)
      if (!length(cpgs)) {
        deps <- setdiff(
          as.character(entry[["clock_inputs"]] %||% character()),
          cid
        )
        deps <- Filter(
          function(d) {
            identical(
              catalog[["clocks"]][[d]][["group_id"]],
              entry[["group_id"]]
            )
          },
          deps
        )
        cpgs <- unique(unlist(
          lapply(deps, panel_of, stack = c(stack, cid)),
          use.names = FALSE
        ))
      }
      cpgs <- unique(cpgs[nzchar(cpgs) & !is.na(cpgs)])
      memo[[cid]] <- cpgs
      cpgs
    }

    for (cid in ids) {
      cpgs <- panel_of(cid, character())
      assert_declared_n_cpgs(catalog[["clocks"]][[cid]], cpgs, cid)
      catalog[["clocks"]][[cid]][["probe_sets"]] <- add_scoring_probe_set(
        catalog[["clocks"]][[cid]],
        cpgs
      )
    }
  }
  catalog
}

# external asset encoding

# named numeric -> double[n] in cpgs order.
align_double <- function(x, cpgs, label) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(label, ": expected a named numeric vector", call. = FALSE)
  }
  # named numerics only -- unnamed would mis-align a coefficient column
  if (is.null(names(x))) {
    stop(label, ": unnamed numeric cannot be aligned to cpgs", call. = FALSE)
  }
  if (!setequal(names(x), cpgs)) {
    stop(label, ": probe set differs from canonical cpgs", call. = FALSE)
  }
  if (identical(names(x), cpgs)) {
    as.double(unname(x))
  } else {
    as.double(unname(x[cpgs]))
  }
}

# cbind named numeric tensors into a double matrix.
cbind_aligned <- function(tensors, rels, cpgs, col_names) {
  if (length(rels) != length(col_names)) {
    stop("cbind_aligned: rels and col_names length mismatch", call. = FALSE)
  }
  cols <- vector("list", length(rels))
  for (i in seq_along(rels)) {
    rel <- rels[[i]]
    if (is.null(tensors[[rel]])) {
      stop(
        "Missing tensor for matrix column ",
        col_names[[i]],
        ": ",
        rel,
        call. = FALSE
      )
    }
    cols[[i]] <- align_double(tensors[[rel]], cpgs, rel)
  }
  mat <- do.call(cbind, cols)
  storage.mode(mat) <- "double"
  colnames(mat) <- col_names
  mat
}

# every cpg-keyed tensor in the group (coefs, components, shared impute ref)
group_cpg_keyed_rels <- function(catalog, gid) {
  rels <- character()
  for (cid in names(catalog[["clocks"]])) {
    e <- catalog[["clocks"]][[cid]]
    if (!identical(as.character(e[["group_id"]] %||% ""), gid)) {
      next
    }
    if (!is.null(e[["coef_path"]])) {
      rels <- c(rels, as.character(e[["coef_path"]]))
    }
    for (comp in e[["components"]] %||% list()) {
      if (identical(comp[["row_key"]], "cpg")) {
        rels <- c(rels, vendored_path(comp[["file"]], "components[].file", cid))
      }
    }
    ref <- e[["imputation"]][["ref"]]
    if (!is.null(ref)) {
      rels <- c(rels, vendored_path(ref, "imputation.ref", cid))
    }
  }
  unique(rels)
}

# canonical probe order for an external pack (declared shared row order)
pack_canonical_cpgs <- function(tensors, catalog, gid) {
  rels <- group_cpg_keyed_rels(catalog, gid)
  if (!length(rels)) {
    stop(gid, ": group declares no cpg-keyed tensors", call. = FALSE)
  }
  orders <- lapply(rels, function(r) {
    tensor_row_keys(tensors[[r]], paste0(gid, " cpg-keyed tensor ", r))
  })
  cpgs <- orders[[1L]]
  if (anyDuplicated(cpgs)) {
    stop(
      gid,
      ": cpg-keyed tensor ",
      rels[[1L]],
      " repeats a CpG",
      call. = FALSE
    )
  }
  for (i in seq_along(orders)) {
    if (identical(orders[[i]], cpgs)) {
      next
    }
    stop(
      gid,
      ": cpg-keyed tensor ",
      rels[[i]],
      if (setequal(orders[[i]], cpgs)) {
        " lists the same CpGs in a different order than "
      } else {
        " covers a different probe set than "
      },
      rels[[1L]],
      call. = FALSE
    )
  }

  # a declared bare probe list over exactly this order is redundant with `cpgs`.
  drop_lists <- names(tensors)[vapply(
    tensors,
    function(x) {
      is.character(x) && is.null(dim(x)) && identical(as.character(x), cpgs)
    },
    logical(1L)
  )]
  list(cpgs = cpgs, drop_lists = drop_lists)
}

# leftover tensors that are not cpg-aligned bulk.
residual_tensors <- function(tensors, used_rels) {
  keep <- setdiff(names(tensors), used_rels)
  if (!length(keep)) {
    return(list())
  }
  tensors[keep]
}

# declared component file by name (absent or ambiguous is an error)
component_file <- function(entry, name, cid) {
  hits <- Filter(
    function(c) identical(as.character(c[["name"]]), as.character(name)),
    entry[["components"]] %||% list()
  )
  if (length(hits) != 1L) {
    stop(
      cid,
      ": expected exactly 1 component named '",
      name,
      "', found ",
      length(hits),
      call. = FALSE
    )
  }
  vendored_path(hits[[1L]][["file"]], "components[].file", cid)
}

# every member clock of a group, in catalog order
group_member_ids <- function(catalog, gid) {
  names(catalog[["clocks"]])[vapply(
    catalog[["clocks"]],
    function(e) identical(as.character(e[["group_id"]] %||% ""), gid),
    logical(1L)
  )]
}

# member clocks that own a coefficient file, clock_id -> coef_path
member_coef_files <- function(catalog, gid) {
  ids <- names(catalog[["clocks"]])[vapply(
    catalog[["clocks"]],
    function(e) {
      identical(as.character(e[["group_id"]] %||% ""), gid) &&
        !is.null(e[["coef_path"]])
    },
    logical(1L)
  )]
  # radix, not the locale collation -- this order feeds payload_hash.
  ids <- sort(ids, method = "radix")
  stats::setNames(
    vapply(
      ids,
      function(i) as.character(catalog[["clocks"]][[i]][["coef_path"]]),
      character(1L)
    ),
    ids
  )
}

# the one vendor-mean fill table the group's members declare
group_impute_ref <- function(catalog, gid) {
  refs <- character()
  for (cid in names(catalog[["clocks"]])) {
    e <- catalog[["clocks"]][[cid]]
    if (!identical(as.character(e[["group_id"]] %||% ""), gid)) {
      next
    }
    ref <- e[["imputation"]][["ref"]]
    if (!is.null(ref)) {
      refs <- c(refs, vendored_path(ref, "imputation.ref", cid))
    }
  }
  refs <- unique(refs)
  if (length(refs) != 1L) {
    stop(
      gid,
      ": expected exactly 1 declared imputation ref across members, found ",
      length(refs),
      call. = FALSE
    )
  }
  refs
}

encode_pcclocks <- function(bundle, catalog) {
  tensors <- bundle[["tensors"]]
  gid <- "PCClocks"
  resolved <- pack_canonical_cpgs(tensors, catalog, gid)
  cpgs <- resolved$cpgs

  # every PC* vector is a sibling clock's own coefficient file.
  coefs <- member_coef_files(catalog, gid)
  if (!length(coefs)) {
    stop(gid, ": no member coefficient tensors declared", call. = FALSE)
  }
  coef_rels <- unname(coefs)
  col_names <- names(coefs)

  impute_rel <- group_impute_ref(catalog, gid)

  used <- c(coef_rels, impute_rel, resolved$drop_lists)
  bundle[["cpgs"]] <- cpgs
  bundle[["coefficient_matrix"]] <- cbind_aligned(
    tensors,
    coef_rels,
    cpgs,
    col_names
  )
  bundle[["impute"]] <- align_double(tensors[[impute_rel]], cpgs, impute_rel)
  bundle[["tensors"]] <- residual_tensors(tensors, used)
  bundle[["encoding"]] <- "canonical_matrices"
  bundle
}

# stack operand -> column label, same rule the accessors use when reading.
systemsage_stack_labels <- function(entry) {
  stack <- Filter(
    function(s) identical(s[["op"]], "stack"),
    entry[["recipe"]] %||% list()
  )
  if (length(stack) != 1L) {
    stop(
      "SystemsAge: expected exactly 1 stack step, found ",
      length(stack),
      call. = FALSE
    )
  }
  mc_runtime[["stack_label_map"]](stack[[1L]], "SystemsAge")
}

# component a linear step multiplies to produce each labelled stack column
systemsage_system_components <- function(entry) {
  map <- systemsage_stack_labels(entry)
  steps <- Filter(
    function(s) {
      identical(s[["op"]], "linear") &&
        !is.null(s[["out"]]) &&
        as.character(s[["out"]]) %in% names(map)
    },
    entry[["recipe"]] %||% list()
  )
  stats::setNames(
    vapply(steps, function(s) as.character(s[["coef"]]), character(1L)),
    vapply(
      steps,
      function(s) unname(map[[as.character(s[["out"]])]]),
      character(1L)
    )
  )
}

encode_systemsage <- function(bundle, catalog) {
  tensors <- bundle[["tensors"]]
  gid <- "SystemsAge"
  resolved <- pack_canonical_cpgs(tensors, catalog, gid)
  cpgs <- resolved$cpgs

  # organs = sibling coef files, systems = composite's own components
  organs <- member_coef_files(catalog, gid)
  labels <- names(organs)
  composite <- catalog[["clocks"]][[gid]]
  sys_comp <- systemsage_system_components(composite)

  # every stack label must be a member clock id (catches Age_prediction)
  declared <- setdiff(unname(systemsage_stack_labels(composite)), gid)
  members <- setdiff(group_member_ids(catalog, gid), gid)
  if (!setequal(declared, members)) {
    stop(
      gid,
      ": stack column labels and group members disagree -- only in labels: ",
      paste(setdiff(declared, members), collapse = ", "),
      "; only in members: ",
      paste(setdiff(members, declared), collapse = ", "),
      call. = FALSE
    )
  }

  # the labels with a CpG front are exactly the coefficient-owning members
  if (!setequal(labels, names(sys_comp))) {
    stop(
      gid,
      ": organ members and raw_{system} recipe steps disagree -- only in ",
      "members: ",
      paste(setdiff(labels, names(sys_comp)), collapse = ", "),
      "; only in recipe: ",
      paste(setdiff(names(sys_comp), labels), collapse = ", "),
      call. = FALSE
    )
  }
  system_rels <- vapply(
    labels,
    function(l) component_file(composite, sys_comp[[l]], gid),
    character(1L),
    USE.NAMES = FALSE
  )
  organ_rels <- unname(organs)

  age_rel <- component_file(composite, "age_pc_coef", gid)
  impute_rel <- group_impute_ref(catalog, gid)

  for (r in c(organ_rels, system_rels, age_rel, impute_rel)) {
    if (is.null(tensors[[r]])) {
      stop(gid, ": missing tensor ", r, call. = FALSE)
    }
  }

  used <- c(
    organ_rels,
    system_rels,
    age_rel,
    impute_rel,
    resolved$drop_lists
  )
  bundle[["cpgs"]] <- cpgs
  bundle[["organs"]] <- cbind_aligned(tensors, organ_rels, cpgs, labels)
  bundle[["systems"]] <- cbind_aligned(tensors, system_rels, cpgs, labels)
  bundle[["age"]] <- align_double(tensors[[age_rel]], cpgs, age_rel)
  bundle[["impute"]] <- align_double(tensors[[impute_rel]], cpgs, impute_rel)
  bundle[["tensors"]] <- residual_tensors(tensors, used)
  bundle[["encoding"]] <- "canonical_matrices"
  bundle
}

encode_pcbrainage <- function(bundle, catalog) {
  tensors <- bundle[["tensors"]]
  gid <- "PCBrainAge"
  resolved <- pack_canonical_cpgs(tensors, catalog, gid)
  cpgs <- resolved$cpgs

  coefs <- member_coef_files(catalog, gid)
  impute_rel <- group_impute_ref(catalog, gid)
  for (r in c(unname(coefs), impute_rel)) {
    if (is.null(tensors[[r]])) {
      stop(gid, ": missing ", r, call. = FALSE)
    }
  }

  used <- c(unname(coefs), impute_rel, resolved$drop_lists)
  bundle[["cpgs"]] <- cpgs
  bundle[["coefficient_matrix"]] <- cbind_aligned(
    tensors,
    unname(coefs),
    cpgs,
    names(coefs)
  )
  bundle[["impute"]] <- align_double(tensors[[impute_rel]], cpgs, impute_rel)
  bundle[["tensors"]] <- residual_tensors(tensors, used)
  bundle[["encoding"]] <- "canonical_matrices"
  bundle
}

# zhang2019 BLUP arm: one dense coef vector. omit policy, no vendored ref.
encode_zhang2019 <- function(bundle, catalog) {
  gid <- "Zhang2019"
  ids <- as.character(bundle[["clocks"]] %||% character())
  if (length(ids) != 1L) {
    stop(
      gid,
      ": expected 1 external member, found ",
      length(ids),
      call. = FALSE
    )
  }
  entry <- catalog[["clocks"]][[ids[[1L]]]]
  rel <- entry[["coef_path"]]
  if (is.null(rel) || is.null(bundle[["tensors"]][[rel]])) {
    stop(gid, ": missing coefficient tensor for ", ids[[1L]], call. = FALSE)
  }

  bundle[["cpgs"]] <- resolved_scoring_cpgs(entry, ids[[1L]])
  # the coefficients stay a raw tensor, keyed by the coef_path clock_coefs() reads
  bundle[["encoding"]] <- "raw_tensors"
  bundle
}

encode_external_asset <- function(bundle, catalog) {
  gid <- bundle[["group_id"]] %||% NA_character_
  if (identical(gid, "PCClocks")) {
    encode_pcclocks(bundle, catalog)
  } else if (identical(gid, "SystemsAge")) {
    encode_systemsage(bundle, catalog)
  } else if (identical(gid, "PCBrainAge")) {
    encode_pcbrainage(bundle, catalog)
  } else if (identical(gid, "Zhang2019")) {
    encode_zhang2019(bundle, catalog)
  } else {
    stop("No external encoding for group_id=", gid, call. = FALSE)
  }
}

# runtime registry row for mc_provenance
external_asset_registry_row <- function(a) {
  # release_tag is <group>-<hash> (set by build_external_assets)
  tag <- as.character(a[["release_tag"]] %||% NA_character_)
  if (is.na(tag) || !nzchar(tag)) {
    stop(
      "external asset for group '",
      a[["group_id"]],
      "' has no release_tag",
      call. = FALSE
    )
  }
  list(
    group_id = a[["group_id"]],
    release_tag = tag,
    file = a[["file"]],
    size_bytes = a[["size_bytes"]],
    encoding = a[["encoding"]],
    encoding_version = a[["encoding_version"]],
    n_clocks = a[["n_clocks"]],
    n_cpgs = a[["n_cpgs"]]
  )
}

# sysdata

# flat per-clock index for list_clocks() filters.
build_index <- function(catalog) {
  clocks <- catalog[["clocks"]]
  ids <- names(clocks)
  citations <- catalog[["citations"]]

  scal <- function(field, default = NA_character_) {
    unname(vapply(
      clocks,
      function(e) {
        v <- e[[field]]
        if (is.null(v) || !length(v)) default else as.character(v)[[1L]]
      },
      character(1L)
    ))
  }
  lgl <- function(field) {
    unname(vapply(clocks, function(e) isTRUE(e[[field]]), logical(1L)))
  }

  # sample-axis boundary position (list_clocks filter is derived from it)
  cross_at <- unname(vapply(
    clocks,
    function(e) {
      v <- e[["cross_sample_at"]]
      if (is.null(v) || !length(v)) NA_integer_ else as.integer(v)[[1L]]
    },
    integer(1L)
  ))

  # citation count (alias resolves through donor)
  n_citations <- unname(vapply(
    clocks,
    function(e) {
      id <- as.character(e[["donor_clock_id"]] %||% e[["clock_id"]])
      sum(citations[["clock_id"]] == id)
    },
    integer(1L)
  ))

  idx <- data.frame(
    clock_id = ids,
    group_id = scal("group_id"),
    # `kind` holds package-minted classifications (sex-routed aliases)
    kind = scal("kind", "clock"),
    weights_format = scal("weights_format"),
    computation_type = scal("computation_type"),
    output_transform = scal("output_transform", "identity"),
    imputation_policy = scal("imputation_policy"),
    cross_sample_at = cross_at,
    batch_dependent = !is.na(cross_at),
    external_group = lgl("external_group"),
    # `pmid` is the primary paper only -- a scalar cannot hold the citation set
    pmid = scal("pmid"),
    n_citations = n_citations,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  idx$covariates_required <- unname(lapply(clocks, function(e) {
    as.character(e[["covariates_required"]] %||% character())
  }))
  idx$n_covariates <- lengths(idx$covariates_required)

  idx
}

# drop minted scoring probe_sets for external groups (panel lives in the pack)
drop_external_probe_cpgs <- function(clocks) {
  for (cid in names(clocks)) {
    if (!isTRUE(clocks[[cid]][["external_group"]])) {
      next
    }
    ps <- clocks[[cid]][["probe_sets"]]
    for (k in seq_along(ps)) {
      ps[[k]][["cpgs"]] <- NULL
    }
    clocks[[cid]][["probe_sets"]] <- ps
  }
  clocks
}

# name items by a declared key, stop on collisions (or missing keys when total).
# total=FALSE leaves unkeyed elements unnamed.
key_declarations <- function(items, field, cid, what, total = TRUE) {
  if (!length(items)) {
    return(items)
  }
  keys <- vapply(
    items,
    function(x) {
      v <- if (is.list(x)) x[[field]] else NULL
      if (is.null(v) || !length(v)) "" else as.character(v)[[1L]]
    },
    character(1L)
  )
  if (total && !all(nzchar(keys))) {
    stop(
      cid,
      ": ",
      sum(!nzchar(keys)),
      " of ",
      length(keys),
      " ",
      what,
      " declare no '",
      field,
      "'",
      call. = FALSE
    )
  }
  dup <- unique(keys[nzchar(keys) & duplicated(keys)])
  if (length(dup)) {
    stop(
      cid,
      ": ",
      what,
      " share a '",
      field,
      "': ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }
  stats::setNames(items, keys)
}

# recipe ops declared at most once (linear / linear_mean may repeat per stack col).
SINGLETON_OPS <- c("stack", "poly")

# probe_set norm roles and stack operand namespaces from the accessors.
NORM_ROLES <- mc_runtime[["NORM_ROLES"]]
STACK_NAMESPACES <- mc_runtime[["STACK_NAMESPACES"]]

# shape invariants asserted once at sync (named by clock).
assert_catalog_shape <- function(clocks) {
  for (cid in names(clocks)) {
    entry <- clocks[[cid]]
    recipe <- entry[["recipe"]] %||% list()
    ops <- vapply(
      recipe,
      function(s) as.character(s[["op"]] %||% "")[[1L]],
      character(1L)
    )

    for (op in SINGLETON_OPS) {
      n <- sum(ops == op)
      if (n > 1L) {
        stop(
          cid,
          " declares ",
          n,
          " '",
          op,
          "' ops; at most one is supported",
          call. = FALSE
        )
      }
    }

    n_bg <- length(intersect(
      NORM_ROLES,
      names(entry[["probe_sets"]])
    ))
    if (n_bg > 1L) {
      stop(
        cid,
        " declares ",
        n_bg,
        " normalization background probe_sets; at most one is supported",
        call. = FALSE
      )
    }

    # declared column labels must cover every operand, or the two zip wrong
    for (step in recipe[ops == "stack"]) {
      declared <- step[["columns"]]
      if (is.null(declared)) {
        next
      }
      n_op <- sum(lengths(lapply(
        STACK_NAMESPACES,
        function(k) as.character(unlist(step[[k]] %||% character()))
      )))
      n_lab <- length(as.character(unlist(declared)))
      if (n_lab != n_op) {
        stop(
          cid,
          ": stack declares ",
          n_lab,
          " column label(s) for ",
          n_op,
          " operand(s)",
          call. = FALSE
        )
      }
    }
  }
  clocks
}

# key the three per-clock declaration lists the accessors look up by name
key_catalog_lists <- function(clocks) {
  for (cid in names(clocks)) {
    entry <- clocks[[cid]]
    entry[["components"]] <- key_declarations(
      entry[["components"]],
      "name",
      cid,
      "components"
    )
    entry[["probe_sets"]] <- key_declarations(
      entry[["probe_sets"]],
      "role",
      cid,
      "probe_sets"
    )
    entry[["shared"]] <- key_declarations(
      entry[["shared"]],
      "name",
      cid,
      "shared tensors"
    )
    entry[["recipe"]] <- key_declarations(
      entry[["recipe"]],
      "out",
      cid,
      "recipe steps",
      total = FALSE
    )
    clocks[[cid]] <- entry
  }
  clocks
}

build_sysdata <- function(
  repo_path,
  catalog,
  ship_groups,
  external_assets = NULL
) {
  message(
    "sync: building shipped bundles for ",
    length(ship_groups),
    " groups..."
  )
  bundles <- build_group_bundles(repo_path, catalog, ship_groups)
  catalog <- resolve_group_scoring_probe_sets(catalog, bundles)
  catalog <- attach_sex_routed_aliases(catalog)
  catalog[["clocks"]] <- trim_build_only_fields(catalog[["clocks"]])
  catalog[["clocks"]] <- key_catalog_lists(catalog[["clocks"]])
  catalog[["clocks"]] <- assert_catalog_shape(catalog[["clocks"]])

  mc_catalog <- drop_external_probe_cpgs(catalog[["clocks"]])
  mc_groups <- catalog[["groups"]]
  mc_bundles <- bundles
  mc_index <- build_index(catalog)
  mc_citations <- catalog[["citations"]]
  ext_reg <- NULL
  if (!is.null(external_assets) && length(external_assets)) {
    ext_reg <- lapply(external_assets, external_asset_registry_row)
  }
  mc_provenance <- list(
    source_git_sha = catalog[["source_git_sha"]],
    schema_version = catalog[["schema_version"]],
    n_clocks = catalog[["n_clocks"]],
    n_ship_groups = length(ship_groups),
    ship_groups = ship_groups,
    external_groups = split_group_ids(catalog, TRUE),
    external_assets = ext_reg
  )

  usethis::use_data(
    mc_catalog,
    mc_groups,
    mc_bundles,
    mc_index,
    mc_citations,
    mc_provenance,
    internal = TRUE,
    overwrite = TRUE,
    compress = "xz"
  )

  path <- file.path("R", "sysdata.rda")
  sz <- file.info(path)$size
  message(sprintf(
    "sync: wrote %s (%.1f KB; %d clocks indexed)",
    path,
    sz / 1024,
    nrow(mc_index)
  ))
  invisible(list(
    path = path,
    size_bytes = sz,
    objects = c(
      "mc_catalog",
      "mc_groups",
      "mc_bundles",
      "mc_index",
      "mc_citations",
      "mc_provenance"
    )
  ))
}

# content-addressed external packs

# low ZSTD, no shuffle.
QS2_COMPRESS_LEVEL <- 1L
QS2_SHUFFLE <- FALSE

# canonical pack for hash + qs_save (weights only)
stable_external_payload <- function(bundle) {
  for (f in EXTERNAL_PIN_FIELDS) {
    bundle[[f]] <- NULL
  }

  # radix sort -- these orders feed payload_hash
  tensors <- bundle[["tensors"]] %||% list()
  if (length(tensors) && !is.null(names(tensors))) {
    tensors <- tensors[sort(names(tensors), method = "radix")]
  }

  clocks <- as.character(bundle[["clocks"]] %||% character())
  if (length(clocks)) {
    clocks <- sort(unique(clocks), method = "radix")
  }

  out <- list(
    encoding_version = as.integer(
      bundle[["encoding_version"]] %||% EXTERNAL_ENCODING_VERSION
    ),
    encoding = as.character(bundle[["encoding"]] %||% "canonical_matrices"),
    group_id = as.character(bundle[["group_id"]] %||% NA_character_),
    clocks = clocks,
    schema_version = bundle[["schema_version"]],
    cpgs = bundle[["cpgs"]],
    coefficient_matrix = bundle[["coefficient_matrix"]],
    organs = bundle[["organs"]],
    systems = bundle[["systems"]],
    age = bundle[["age"]],
    impute = bundle[["impute"]],
    tensors = if (length(tensors)) tensors else NULL
  )
  out[!vapply(out, is.null, logical(1L))]
}

# Stable across R versions. serialize() writes a 14 byte header whose bytes
# 7 to 10 carry the writer's R version, so hashing the raw stream re-addressed
# every pack on an R upgrade. digest's skip = "auto" drops that header.
payload_hash_of <- function(payload) {
  digest::digest(payload, algo = "sha256", serializeVersion = 2L)
}

# gitHub release target

# `owner/repo` from a plain slug, https URL, or scp-style remote
parse_github_owner_repo <- function(url) {
  x <- trimws(as.character(url %||% ""))
  if (!nzchar(x)) {
    return(NULL)
  }
  x <- sub("/+$", "", x) # trailing slash
  x <- sub("\\.git$", "", x)
  x <- sub("^[A-Za-z][A-Za-z0-9+.-]*://", "", x) # scheme
  x <- sub("^[^/@]+@", "", x) # scp-style user@
  x <- sub("^github\\.com[:/]", "", x) # host (only github.com)
  parts <- strsplit(x, "/", fixed = TRUE)[[1L]]
  # exactly owner/repo, each segment a plain GitHub name
  if (length(parts) != 2L || !all(grepl("^[A-Za-z0-9._-]+$", parts))) {
    return(NULL)
  }
  list(
    owner = parts[[1L]],
    repo = parts[[2L]],
    slug = paste0(parts[[1L]], "/", parts[[2L]])
  )
}

# release target: env MC_RELEASE_REPO, else the package origin remote
package_release_repo <- function() {
  env <- Sys.getenv("MC_RELEASE_REPO", unset = "")
  if (nzchar(env)) {
    parsed <- parse_github_owner_repo(env)
    if (is.null(parsed)) {
      stop(
        "MC_RELEASE_REPO is set to '",
        env,
        "' which is not owner/repo or a GitHub URL. ",
        "The publish target is never guessed and the owner is never inferred -- ",
        "fix the value, or unset it to use git remote origin.",
        call. = FALSE
      )
    }
    return(parsed)
  }
  url <- tryCatch(
    git_value("remote", "get-url", "origin"),
    error = function(e) NA_character_
  )
  parsed <- parse_github_owner_repo(url)
  if (is.null(parsed)) {
    stop(
      "Cannot resolve package GitHub repo for releases. Set MC_RELEASE_REPO=owner/repo ",
      "or configure git remote origin.",
      call. = FALSE
    )
  }
  parsed
}

package_release_target_commitish <- function() {
  env <- Sys.getenv("MC_RELEASE_TARGET", unset = "")
  if (nzchar(env)) {
    return(env)
  }
  br <- tryCatch(
    git_value("rev-parse", "--abbrev-ref", "HEAD"),
    error = function(e) NA_character_
  )
  if (is.na(br) || !nzchar(br) || identical(br, "HEAD")) {
    return(git_value("rev-parse", "HEAD"))
  }
  br
}

# upload via `uv run python data-raw/gh_upload.py` (idempotent by release_tag).
GH_UPLOAD_PY <- file.path("data-raw", "gh_upload.py")

uv_bin <- function() {
  w <- Sys.which("uv")
  if (!nzchar(w)) {
    stop(
      "`uv` not found on PATH (needed to run ",
      GH_UPLOAD_PY,
      " for uploads).",
      call. = FALSE
    )
  }
  w
}

# upload PAT from MC_UPLOAD_PAT only (no GITHUB_PAT fallback)
upload_pat <- function() {
  Sys.getenv("MC_UPLOAD_PAT", unset = "")
}

upload_external_assets <- function(assets) {
  if (!length(assets)) {
    message("sync: no external assets to upload")
    return(invisible(list()))
  }
  pat <- upload_pat()
  if (!nzchar(pat)) {
    stop(
      "upload=TRUE requires a GitHub token in MC_UPLOAD_PAT ",
      "(a fine-grained PAT with Contents:read/write on the release repo). ",
      "Set it in ~/.Renviron. MC_UPLOAD_PAT is the only source: GITHUB_PAT / ",
      "GITHUB_TOKEN / GH_TOKEN are broad tokens other tooling sets, and are ",
      "never used to publish releases.",
      call. = FALSE
    )
  }

  repo <- package_release_repo()
  target <- package_release_target_commitish()

  items <- lapply(assets, function(a) {
    fpath <- a[["path"]] %||% file.path(asset_dir, a[["file"]])
    list(
      group_id = as.character(a[["group_id"]] %||% ""),
      tag = as.character(a[["release_tag"]] %||% ""),
      path = as.character(fs::path_real(fpath)),
      name = as.character(a[["file"]] %||% basename(fpath))
    )
  })

  req <- list(
    slug = repo$slug,
    target_commitish = target,
    assets = unname(items)
  )
  json <- jsonlite::toJSON(req, auto_unbox = TRUE, null = "null")

  message(
    "sync: uploading ",
    length(items),
    " external asset(s) to ",
    repo$slug,
    " (target_commitish=",
    target,
    ") via PyGithub"
  )
  # feed PAT only to the child env, stage JSON for processx stdin
  req_file <- tempfile("uv-gh-req-", fileext = ".json")
  on.exit(unlink(req_file), add = TRUE)
  writeLines(json, req_file)

  proc <- processx::run(
    uv_bin(),
    c("run", "python", GH_UPLOAD_PY),
    stdin = req_file,
    env = c("current", GITHUB_TOKEN = pat, GH_TOKEN = pat),
    error_on_status = FALSE
  )
  if (proc$status != 0L) {
    stop(GH_UPLOAD_PY, " failed:\n", proc$stderr, call. = FALSE)
  }

  # gh_upload.py stdout must be exactly {"results": [...]}
  res <- tryCatch(
    jsonlite::fromJSON(proc$stdout, simplifyVector = FALSE),
    error = function(e) NULL
  )
  results <- if (is.list(res)) res[["results"]] else NULL
  if (!is.list(results)) {
    stop(
      GH_UPLOAD_PY,
      " exited 0 but did not return the documented {\"results\": [...]} object ",
      "on stdout, so what was published cannot be confirmed. stdout was:\n",
      proc$stdout,
      call. = FALSE
    )
  }

  reported <- vapply(
    results,
    function(r) as.character(r[["name"]] %||% NA_character_),
    character(1L)
  )
  unreported <- setdiff(
    vapply(items, function(i) i[["name"]], character(1L)),
    reported
  )
  if (length(unreported)) {
    stop(
      GH_UPLOAD_PY,
      " reported no outcome for: ",
      paste(unreported, collapse = ", "),
      call. = FALSE
    )
  }

  for (r in results) {
    message(
      "sync: ",
      r[["action"]],
      " ",
      r[["name"]],
      " -> ",
      repo$slug,
      " @ tag ",
      r[["tag"]]
    )
  }
  invisible(assets)
}

# build content-addressed external packs: <group>-<payload_hash>.qs2.
build_external_assets <- function(repo_path, catalog, external_groups) {
  fs::dir_create(asset_dir)
  assets <- list()

  for (gid in external_groups) {
    message("sync: building external asset for ", gid, "...")
    raw_bundle <- build_group_bundles(
      repo_path,
      catalog,
      gid,
      external = TRUE
    )[[gid]]
    # resolve probe sets before encoding, adopt catalog so pack and sysdata match
    catalog <- resolve_group_scoring_probe_sets(
      catalog,
      stats::setNames(list(raw_bundle), gid)
    )
    bundle <- encode_external_asset(raw_bundle, catalog)
    bundle[["schema_version"]] <- catalog[["schema_version"]]
    bundle[["encoding_version"]] <- EXTERNAL_ENCODING_VERSION

    payload <- stable_external_payload(bundle)
    phash <- payload_hash_of(payload)
    fname <- sprintf("%s-%s.qs2", tolower(gid), phash)
    fpath <- file.path(asset_dir, fname)
    # tag = filename stem (<group>-<hash>)
    rtag <- sub("\\.qs2$", "", fname)

    qs2::qs_save(
      payload,
      fpath,
      compress_level = QS2_COMPRESS_LEVEL,
      shuffle = QS2_SHUFFLE
    )

    sz <- file.info(fpath)$size
    n_cpgs <- length(payload[["cpgs"]] %||% character())
    message(sprintf(
      "sync: wrote %s (%.2f MB; payload_hash=%s; n_cpgs=%s)",
      fpath,
      sz / 1e6,
      phash,
      n_cpgs
    ))
    assets[[gid]] <- list(
      group_id = gid,
      path = fpath,
      file = fname,
      payload_hash = phash,
      release_tag = rtag,
      size_bytes = as.integer(sz),
      n_clocks = length(payload[["clocks"]] %||% character()),
      n_cpgs = n_cpgs,
      encoding = payload[["encoding"]] %||% "canonical_matrices",
      encoding_version = as.integer(
        payload[["encoding_version"]] %||% EXTERNAL_ENCODING_VERSION
      )
    )
  }

  # prune superseded local staging packs (published releases untouched)
  for (gid in external_groups) {
    keep <- as.character(assets[[gid]][["file"]] %||% "")
    siblings <- list.files(
      asset_dir,
      pattern = paste0("^", tolower(gid), "-[0-9a-f]+\\.qs2$")
    )
    for (f in setdiff(siblings, keep)) {
      unlink(file.path(asset_dir, f))
      message("sync: pruned superseded staging asset ", f)
    }
  }

  list(assets = assets, catalog = catalog)
}

# asset lockfile

read_lockfile <- function() {
  if (!file.exists(LOCKFILE)) {
    return(NULL)
  }
  tryCatch(readRDS(LOCKFILE), error = function(e) NULL)
}

# clock ids whose weights live in a pack
external_clock_ids <- function(catalog) {
  names(catalog[["clocks"]])[
    vapply(
      catalog[["clocks"]],
      function(e) isTRUE(e[["external_group"]]),
      logical(1L)
    )
  ]
}

# per-clock bundle_hash for the external clocks, the pack staleness key.
external_bundle_hashes <- function(catalog) {
  ids <- external_clock_ids(catalog)
  hashes <- vapply(
    catalog[["clocks"]][ids],
    function(e) as.character(e[["bundle_hash"]] %||% NA_character_),
    character(1L)
  )
  hashes[order(names(hashes))]
}

# lockfile hit when encoder code, every external bundle_hash, and packs match.
lockfile_hit <- function(lock, hashes, fingerprint) {
  if (is.null(lock)) {
    return(FALSE)
  }
  if (is.na(fingerprint) || !identical(lock$code_fingerprint, fingerprint)) {
    return(FALSE)
  }
  prev <- lock$bundle_hashes
  if (is.null(prev) || !identical(prev, hashes) || anyNA(hashes)) {
    return(FALSE)
  }
  files <- vapply(
    lock$assets %||% list(),
    function(a) as.character(a[["file"]] %||% ""),
    character(1L)
  )
  length(files) > 0L &&
    all(nzchar(files)) &&
    all(file.exists(file.path(asset_dir, files)))
}

# store pre-trim external catalog entries for lockfile reuse.
write_lockfile <- function(sha, hashes, assets, ext_clocks, fingerprint) {
  saveRDS(
    list(
      source_git_sha = sha,
      bundle_hashes = hashes,
      code_fingerprint = fingerprint,
      assets = assets,
      ext_clocks = ext_clocks
    ),
    LOCKFILE
  )
  message("sync: wrote ", LOCKFILE, " (asset lockfile @ ", sha, ")")
}

# main

sync <- function(source_git_sha = NULL, upload = FALSE, force = FALSE) {
  src <- resolve_source(source_git_sha = source_git_sha)
  current_sha <- src$source_git_sha

  message("sync: building catalog @ ", current_sha)
  catalog <- build_catalog(src$path, read_manifest(src$path))
  catalog[["source_git_sha"]] <- current_sha

  # a group with both bundled and external members is in both buckets
  external <- split_group_ids(catalog, TRUE)
  ship <- split_group_ids(catalog, FALSE)
  message(
    "sync: ",
    catalog[["n_clocks"]],
    " clocks; ship groups=",
    length(ship),
    "; external=",
    paste(external, collapse = ", ")
  )

  ext_hashes <- external_bundle_hashes(catalog)
  code_fp <- sync_code_fingerprint()
  lock <- if (isTRUE(force)) NULL else read_lockfile()
  # name the case the bundle_hashes cannot explain, or the rebuild looks arbitrary
  if (
    !is.null(lock) &&
      identical(lock$bundle_hashes, ext_hashes) &&
      !identical(lock$code_fingerprint, code_fp)
  ) {
    message(
      "sync: upstream unchanged but ",
      SYNC_SCRIPT,
      " differs from the lockfile -- rebuilding external assets."
    )
  }
  if (lockfile_hit(lock, ext_hashes, code_fp)) {
    message(
      "sync: external assets unchanged (bundle_hash match) -- reusing ",
      length(lock$assets),
      " cached pack(s), skipping rebuild. Run sync(force = TRUE) to reconcile drift."
    )
    assets <- lock$assets
    # restore external resolved probe sets from the lockfile.
    for (cid in names(lock$ext_clocks)) {
      catalog[["clocks"]][[cid]] <- lock$ext_clocks[[cid]]
    }
  } else {
    ext <- build_external_assets(src$path, catalog, external)
    assets <- ext$assets
    # adopt catalog with resolved external probe sets.
    catalog <- ext$catalog
    write_lockfile(
      current_sha,
      ext_hashes,
      assets,
      catalog[["clocks"]][external_clock_ids(catalog)],
      code_fp
    )
  }

  sys <- build_sysdata(src$path, catalog, ship, external_assets = assets)
  bib <- vendor_bibliography(src$path)

  if (isTRUE(upload)) {
    upload_external_assets(assets)
  } else {
    message(
      "sync: staged ",
      length(assets),
      " external asset(s) under data-raw/assets/ (upload=FALSE)"
    )
  }

  message("sync: done @ ", current_sha)
  invisible(list(
    source_git_sha = current_sha,
    n_clocks = catalog[["n_clocks"]],
    ship_groups = ship,
    external_groups = external,
    sysdata = sys,
    bibliography = bib,
    assets = assets,
    uploaded = isTRUE(upload)
  ))
}

if (interactive()) {
  message(
    "sync.R loaded. Clones/fetches ",
    META_REMOTE,
    "\n",
    "  into data-raw/methylCIPHER-meta, then materializes package data.\n",
    "  sync()\n",
    "  sync(upload = TRUE)\n",
    "  sync(force = TRUE)   # rebuild external assets, ignore the lockfile\n",
    "  sync(source_git_sha = \"dc543a7b...\")"
  )
}
