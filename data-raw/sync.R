# Vendor methylCIPHER-meta into the package. source("data-raw/sync.R"); sync().

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

asset_dir <- file.path("data-raw", "assets")
meta_dir <- file.path("data-raw", "methylCIPHER-meta")

# Gitignored lockfile for external pack rebuild skip.
LOCKFILE <- file.path(asset_dir, "lockfile.rds")

META_REMOTE <- "https://github.com/hhp94/methylCIPHER-meta.git"

# External families as release assets; rest in sysdata.
EXTERNAL_GROUPS <- c("SystemsAge", "PCClocks", "PCBrainAge")

# Bump when pack layout changes (new payload_hash).
EXTERNAL_ENCODING_VERSION <- 3L

# Pin-only fields; excluded from the content-addressed pack.
EXTERNAL_PIN_FIELDS <- c("source_git_sha", "manifest_generated_at_sha")

# Anchored weights/ tensor path pattern.
WEIGHTS_REF_RE <- "^weights/.+\\.(csv\\.gz|csv|[Rr])$"

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
  "external",
  "probe_sets",
  "code_ref",
  "definition",
  "depends_on_clocks",
  "n_cpgs_normalization",
  "covers",
  "license"
)

# Build-time only; stripped after resolution.
CATALOG_BUILD_ONLY_FIELDS <- c("covers", "shared")

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
EXTERNAL_FIELDS <- c(
  "r_package",
  "github",
  "commit",
  "function",
  "model_key",
  "depends"
)
RECIPE_STEP_DROP <- c("note")

# Keep named fields; preserves explicit JSON nulls.
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

prune_fixture <- function(fx) {
  if (is.null(fx) || !is.list(fx)) {
    return(NULL)
  }
  parity <- fx[["parity"]] %||% list()
  stub <- list(
    expected = fx[["expected"]] %||% NULL,
    oracle = fx[["oracle"]] %||% NULL,
    parity_policy = parity[["policy"]] %||% NULL,
    parity_metric = parity[["metric"]] %||% NULL
  )
  stub <- stub[!vapply(stub, is.null, logical(1L))]
  if (!length(stub)) NULL else stub
}

# Declared probe_sets are kept for every weights_format. Dropping them off
# external_package silently discarded DunedinPACE's QN background when it
# migrated to cpg_coefficient, leaving the accessor to guess at the tensor.
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
  if ("external" %in% names(out)) {
    out[["external"]] <- keep_fields(out[["external"]], EXTERNAL_FIELDS)
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

read_papers_csv <- function(repo_path) {
  path <- file.path(repo_path, "bibliography", "papers.csv")
  df <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  list(
    bib_key = stats::setNames(trimws(df[["bib_key"]]), trimws(df[["pmid"]])),
    n = nrow(df)
  )
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

# Clone/fetch meta_dir and checkout a commit (NULL = origin tip).
resolve_source <- function(source_git_sha = NULL) {
  fs::dir_create(dirname(meta_dir))

  # Discard a non-repo meta_dir (interrupted clone) and re-clone.
  is_repo <- dir.exists(file.path(meta_dir, ".git"))
  if (dir.exists(meta_dir) && !is_repo) {
    unlink(meta_dir, recursive = TRUE, force = TRUE)
  }

  if (!is_repo) {
    # Clone to a temp sibling, then rename so a partial clone cannot poison meta_dir.
    message("sync: cloning ", META_REMOTE, " -> ", meta_dir)
    tmp <- paste0(meta_dir, ".tmp-", Sys.getpid())
    unlink(tmp, recursive = TRUE, force = TRUE)
    git_exec("clone", "--filter=blob:none", META_REMOTE, tmp)
    file.rename(tmp, meta_dir)
  } else {
    message("sync: fetching ", META_REMOTE, " into ", meta_dir)
    git_exec("remote", "set-url", "origin", META_REMOTE, dir = meta_dir)
    git_exec("fetch", "origin", "--tags", "--prune", dir = meta_dir)
  }

  if (!is.null(source_git_sha) && nzchar(source_git_sha)) {
    ref <- source_git_sha
  } else {
    ref <- tryCatch(
      git_value(
        "symbolic-ref",
        "--short",
        "refs/remotes/origin/HEAD",
        dir = meta_dir
      ),
      error = function(e) NA_character_
    )
    if (is.na(ref) || !nzchar(ref)) {
      ref <- "origin/master"
    }
  }
  message("sync: checkout ", ref)
  # Force detach checkout; mirror is disposable.
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

# catalog crawl

# List meta files under weights/; clock vs group by basename.
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

# Covariate names from recipe, top-level field, and sex-keyed impute refs.
extract_covariates <- function(meta) {
  flatten_names <- function(x) {
    if (is.null(x)) {
      return(character())
    }
    if (is.character(x)) {
      return(x)
    }
    if (is.list(x)) {
      nms <- names(x)
      if (!is.null(nms) && any(nzchar(nms))) {
        return(nms[nzchar(nms)])
      }
      return(as.character(unlist(x, use.names = FALSE)))
    }
    character()
  }

  covariate_keys <- function(x) {
    nms <- names(x)
    if (is.null(nms)) {
      return(character())
    }
    nms[grepl("covariates$", nms)]
  }

  covs <- flatten_names(meta[["covariates"]])
  for (step in meta[["recipe"]] %||% list()) {
    for (k in covariate_keys(step)) {
      covs <- c(covs, flatten_names(step[[k]]))
    }
    if (identical(step[["op"]], "linear_sex")) {
      covs <- c(covs, "Female")
    }
  }

  for (sp in meta[["sex_params"]] %||% list()) {
    covs <- c(covs, flatten_names(sp[["covariates"]]))
  }
  if (isTRUE(meta[["sex_stratified"]])) {
    covs <- c(covs, "Female")
  }

  ref <- meta[["imputation"]][["ref"]]
  if (is.list(ref) && any(c("female", "male") %in% names(ref))) {
    covs <- c(covs, "Female")
  }

  unique(covs[nzchar(covs) & !is.na(covs)])
}

# Ops where subset scores != subset of full-cohort scores.
extract_batch_ops <- function(meta) {
  ops <- vapply(
    meta[["recipe"]] %||% list(),
    function(s) as.character(s[["op"]] %||% NA_character_),
    character(1L)
  )
  intersect(ops[!is.na(ops)], c("cohort_zscore", "sample_scale"))
}

# Relative path from repo root (forward slashes).
rel_from_repo <- function(abs_path, repo_path) {
  as.character(fs::path_rel(abs_path, start = repo_path))
}

# Repo-relative weights/ path -> absolute.
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

  papers <- read_papers_csv(repo_path)

  clocks <- list()
  for (mp in files$clock) {
    meta <- jsonlite::fromJSON(mp, simplifyVector = FALSE)
    cid <- as.character(meta[["clock_id"]] %||% NA_character_)
    gid <- as.character(meta[["group_id"]] %||% NA_character_)

    wf <- as.character(meta[["weights_format"]] %||% NA_character_)
    batch_ops <- extract_batch_ops(meta)
    covs <- extract_covariates(meta)

    default_coef <- file.path("weights", gid, paste0(cid, ".csv.gz"))
    has_default_coef <- file.exists(file.path(repo_path, default_coef))

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
    if (!is.null(entry[["depends_on_clocks"]])) {
      entry[["depends_on_clocks"]] <- unlist(
        entry[["depends_on_clocks"]],
        use.names = FALSE
      )
    }
    entry[["pmid"]] <- as.character(entry[["pmid"]] %||% NA_character_)
    entry[["bib_key"]] <- unname(papers$bib_key[entry[["pmid"]]])

    entry[["imputation_policy"]] <- as.character(
      entry[["imputation"]][["policy"]] %||%
        meta[["imputation"]][["policy"]] %||%
        NA_character_
    )
    entry[["covariates_required"]] <- covs
    entry[["batch_ops"]] <- batch_ops
    entry[["batch_dependent"]] <- length(batch_ops) > 0L
    entry[["external_group"]] <- gid %in% EXTERNAL_GROUPS
    entry[["fixture"]] <- prune_fixture(meta[["fixture"]])
    entry[["meta_path"]] <- rel_from_repo(mp, repo_path)
    entry[["coef_path"]] <- if (
      identical(wf, "cpg_coefficient") && has_default_coef
    ) {
      default_coef
    } else {
      NULL
    }

    clocks[[cid]] <- entry
  }

  list(
    clocks = clocks,
    groups = groups,
    source_git_sha = NA_character_,
    schema_version = manifest[["schema_version"]],
    n_clocks = length(clocks)
  )
}

# tensor IO

# Load weights/*.csv.gz as named numeric, character vector, or data.frame.
read_tensor_csv <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) {
    return(df)
  }
  if (ncol(df) == 1L) {
    return(as.character(df[[1L]]))
  }
  # Two-col key/value only when col2 is numeric; else keep as data.frame.
  if (ncol(df) == 2L && is.numeric(df[[2L]])) {
    key <- as.character(df[[1L]])
    if (anyDuplicated(key)) {
      stop(
        "Duplicate keys in ",
        path,
        " -- a named vector would silently collapse them",
        call. = FALSE
      )
    }
    return(stats::setNames(df[[2L]], key))
  }
  df
}

# Collect every weights/ path referenced anywhere in a meta JSON document.
collect_weights_refs <- function(x) {
  found <- character()
  walk <- function(node) {
    if (is.character(node)) {
      found <<- c(found, node[grepl(WEIGHTS_REF_RE, node)])
    } else if (is.list(node)) {
      lapply(node, walk)
    }
    invisible(NULL)
  }
  walk(x)
  unique(found[nzchar(found) & !is.na(found)])
}

collect_file_refs <- function(entry, repo_path) {
  meta_abs <- resolve_repo_rel(repo_path, entry[["meta_path"]])
  raw <- jsonlite::fromJSON(meta_abs, simplifyVector = FALSE)
  paths <- collect_weights_refs(raw)
  if (!is.null(entry[["coef_path"]])) {
    paths <- c(paths, entry[["coef_path"]])
  }
  unique(paths[nzchar(paths) & !is.na(paths)])
}

# custom-format payloads

# A `custom` group ships frozen author artifacts -- the code is the declared
# definition, so no meta `weights/` ref points at the parameter blob. One loader
# per custom group materializes it as an ordinary bundle tensor; the paired
# component declaration is grafted onto the clock entry so probe-set resolution
# and the accessors read it like any other cpg-keyed tensor.

MIAGE_PANEL_FILE <- "weights/MiAge/Additional_File1.csv.gz"
MIAGE_PARAMS_FILE <- "weights/MiAge/site_specific_parameters.Rdata"

# MiAge site-specific (b, c, d), keyed by the 268-CpG panel in panel order.
miage_site_parameters <- function(repo_path) {
  af1 <- utils::read.csv(
    gzfile(resolve_repo_rel(repo_path, MIAGE_PANEL_FILE)),
    stringsAsFactors = FALSE
  )
  cpgs <- as.character(af1[[1L]])
  env <- new.env(parent = emptyenv())
  load(resolve_repo_rel(repo_path, MIAGE_PARAMS_FILE), envir = env)

  # (b, c, d) align positionally with the panel and nothing downstream can
  # re-check it, so this is the one place a silent mis-key would get caught.
  bcd <- lapply(env$methyl.age[1:3], as.numeric)
  if (any(lengths(bcd) != length(cpgs))) {
    stop(
      "MiAge: site-specific (b, c, d) do not align 1:1 with the ",
      length(cpgs),
      "-CpG panel",
      call. = FALSE
    )
  }

  stats::setNames(
    list(data.frame(
      cpg = cpgs,
      b = bcd[[1L]],
      c = bcd[[2L]],
      d = bcd[[3L]],
      stringsAsFactors = FALSE
    )),
    MIAGE_PARAMS_FILE
  )
}

CUSTOM_GROUPS <- list(
  MiAge = list(
    tensors = miage_site_parameters,
    components = list(
      MiAge = list(list(
        name = "site_parameters",
        file = MIAGE_PARAMS_FILE,
        row_key = "cpg"
      ))
    )
  )
)

# Extra tensors a custom group contributes beyond its meta-referenced files.
custom_group_tensors <- function(gid, repo_path) {
  spec <- CUSTOM_GROUPS[[gid]]
  if (is.null(spec)) {
    return(list())
  }
  spec$tensors(repo_path)
}

# Graft custom-group component declarations onto their clock entries.
attach_custom_components <- function(clocks) {
  for (spec in CUSTOM_GROUPS) {
    for (cid in names(spec$components)) {
      clocks[[cid]][["components"]] <- c(
        clocks[[cid]][["components"]],
        spec$components[[cid]]
      )
    }
  }
  clocks
}

# sex-routed aliases

# Upstream resolves sex in the clock_id and records the pairing as
# `_group.meta.json` routing.sex. The un-suffixed stem is the callable a user
# actually wants, so it is minted here rather than asked for upstream -- adding
# it back to the meta would re-open the split upstream deliberately made. The
# alias owns no weights and no panel: it selects a member per sample.
SEX_SUFFIX <- c(female = "_Female", male = "_Male")

# {stem: {female = id, male = id}} for one group; empty when it does not route.
group_sex_routes <- function(gside) {
  routing <- gside[["routing"]][["sex"]]
  if (is.null(routing)) {
    return(list())
  }
  stems <- list()
  for (sx in names(SEX_SUFFIX)) {
    for (id in as.character(unlist(routing[[sx]] %||% character()))) {
      stem <- sub(paste0(SEX_SUFFIX[[sx]], "$"), "", id)
      if (identical(stem, id)) {
        stop(
          "routing.sex: member '",
          id,
          "' lacks the expected ",
          SEX_SUFFIX[[sx]],
          " suffix.",
          call. = FALSE
        )
      }
      stems[[stem]][[sx]] <- id
    }
  }
  # A stem missing a sex could not route every sample, so it is an error here
  # rather than an NA column at scoring time.
  partial <- names(Filter(
    function(r) !setequal(names(r), names(SEX_SUFFIX)),
    stems
  ))
  if (length(partial)) {
    stop(
      "routing.sex: stem(s) without both sexes: ",
      paste(partial, collapse = ", "),
      call. = FALSE
    )
  }
  stems
}

# Mint one alias clock per routed stem. Runs after probe-set resolution so an
# alias never acquires a panel -- no panel means no coverage gate and no
# coverage record, which is the point: its members' panels are disjoint over
# samples, so any union would describe no sample.
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
      donor <- catalog[["clocks"]][[routes[[stem]]$female]]
      catalog[["clocks"]][[stem]] <- list(
        clock_id = stem,
        group_id = gid,
        weights_format = "routed",
        computation_type = "sex_routed",
        output_transform = "identity",
        normalization = donor[["normalization"]],
        routing = routes[[stem]],
        depends_on_clocks = c(routes[[stem]]$female, routes[[stem]]$male),
        covariates_required = "Female",
        imputation_policy = donor[["imputation_policy"]],
        pmid = donor[["pmid"]],
        bib_key = donor[["bib_key"]],
        license = donor[["license"]],
        batch_dependent = FALSE,
        external_group = FALSE
      )
    }
  }
  catalog[["n_clocks"]] <- length(catalog[["clocks"]])
  catalog
}

# materialize

# Per-group tensor payload for the given group_ids.
build_group_bundles <- function(repo_path, catalog, group_ids) {
  bundles <- list()
  for (gid in group_ids) {
    member_ids <- names(catalog[["clocks"]])[
      vapply(
        catalog[["clocks"]],
        function(c) identical(c[["group_id"]], gid),
        logical(1L)
      )
    ]
    members <- catalog[["clocks"]][member_ids]
    rels <- character()
    for (entry in members) {
      rels <- c(rels, collect_file_refs(entry, repo_path))
    }
    gside <- catalog[["groups"]][[gid]]
    if (!is.null(gside[["shared_tensors"]])) {
      rels <- c(rels, collect_weights_refs(gside[["shared_tensors"]]))
    }

    tensors <- list()
    for (rel in unique(rels)) {
      abs <- resolve_repo_rel(repo_path, rel)
      # Missing refs would ship intercept-only scores; hard-error instead.
      if (!file.exists(abs)) {
        stop(
          "Referenced tensor missing from snapshot: ",
          rel,
          " (group ",
          gid,
          ")",
          call. = FALSE
        )
      }
      if (grepl("\\.[Rr]$", rel)) {
        tensors[[rel]] <- list(
          type = "r_source",
          text = readLines(abs, warn = FALSE)
        )
      } else {
        tensors[[rel]] <- read_tensor_csv(abs)
      }
    }

    custom <- custom_group_tensors(gid, repo_path)
    tensors[names(custom)] <- custom

    bundles[[gid]] <- list(
      group_id = gid,
      clocks = member_ids,
      tensors = tensors
    )
  }
  bundles
}

# scoring CpG resolution
# Materialize each clock's scoring CpGs as probe_sets {name, role, cpgs}.

# Row labels of one loaded tensor (cpg col or names).
tensor_row_keys <- function(t) {
  if (is.null(t)) {
    return(character())
  }
  if (is.data.frame(t)) {
    col <- if ("cpg" %in% names(t)) "cpg" else 1L
    return(as.character(t[[col]]))
  }
  if (is.numeric(t) && !is.null(names(t))) {
    return(names(t))
  }
  if (is.character(t) && is.null(dim(t))) {
    return(t)
  }
  character()
}

# Materialize one probe_set from a file pointer.
materialize_probe_set <- function(ps, tensors, cid = NA_character_) {
  where <- paste0(
    "probe_set '",
    ps[["name"]] %||% "?",
    "' (role ",
    ps[["role"]] %||% "?",
    ") of clock ",
    cid
  )
  if (is.null(ps[["file"]]) || !nzchar(as.character(ps[["file"]]))) {
    stop(where, " has no `file` pointer.", call. = FALSE)
  }
  if (is.null(tensors[[ps[["file"]]]])) {
    stop(
      where,
      " points at a tensor absent from the bundle: ",
      ps[["file"]],
      call. = FALSE
    )
  }
  cpgs <- tensor_row_keys(tensors[[ps[["file"]]]])
  if (!length(cpgs)) {
    stop(where, " resolved to zero CpGs from ", ps[["file"]], call. = FALSE)
  }
  if (anyDuplicated(cpgs)) {
    stop(
      where,
      " contains duplicate CpG id(s): ",
      ps[["file"]],
      call. = FALSE
    )
  }
  # `file` rides along so accessors can fetch the tensor's values, not just its
  # row keys (Dunedin's QN target means).
  list(
    name = ps[["name"]],
    role = ps[["role"]],
    file = ps[["file"]],
    cpgs = cpgs
  )
}

# Union of a clock's own cpg-keyed components.
own_component_cpgs <- function(entry, tensors) {
  cpgs <- character()
  for (comp in entry[["components"]] %||% list()) {
    if (!identical(comp[["row_key"]], "cpg")) {
      next
    }
    cpgs <- c(cpgs, tensor_row_keys(tensors[[comp[["file"]]]]))
  }
  unique(cpgs[nzchar(cpgs) & !is.na(cpgs)])
}

# Group shared_tensors entry that is a bare one-column probe list.
group_shared_cpg_list <- function(gside, tensors) {
  for (rel in gside[["shared_tensors"]] %||% character()) {
    t <- tensors[[rel]]
    if (is.character(t) && is.null(dim(t)) && length(t) > 0L) {
      return(as.character(t))
    }
  }
  character()
}

# Scoring CpGs already on an entry, or character(0).
resolved_scoring_cpgs <- function(entry) {
  hits <- Filter(
    function(p) identical(p[["role"]], "scoring") && length(p[["cpgs"]]) > 0L,
    entry[["probe_sets"]] %||% list()
  )
  unique(unlist(
    lapply(hits, function(p) p[["cpgs"]]),
    use.names = FALSE
  )) %||%
    character()
}

add_scoring_probe_set <- function(entry, cpgs) {
  c(
    entry[["probe_sets"]] %||% list(),
    list(list(name = "scoring_derived", role = "scoring", cpgs = unique(cpgs)))
  )
}

# Tiers a clock can resolve without any other clock being resolved first.
own_scoring_cpgs <- function(entry, tensors) {
  if (!is.null(entry[["coef_path"]])) {
    cpgs <- tensor_row_keys(tensors[[entry[["coef_path"]]]])
    if (length(cpgs)) {
      return(cpgs)
    }
  }
  own_component_cpgs(entry, tensors)
}

# Materialize probe_sets, then fill missing scoring panels in tiers: the clock's
# own tensors, then the union over the clocks it consumes, then the group's
# shared bare CpG list. Dependencies outrank the shared list -- taking the
# shared list first gave the DNAmFitAge composites the family-wide 627-probe
# prep panel instead of what they actually consume.
resolve_group_scoring_probe_sets <- function(catalog, bundles) {
  for (gid in names(bundles)) {
    tensors <- bundles[[gid]][["tensors"]]
    for (cid in bundles[[gid]][["clocks"]]) {
      entry <- catalog[["clocks"]][[cid]]
      if (length(entry[["probe_sets"]])) {
        catalog[["clocks"]][[cid]][["probe_sets"]] <- lapply(
          entry[["probe_sets"]],
          materialize_probe_set,
          tensors = tensors,
          cid = cid
        )
      }
      if (length(resolved_scoring_cpgs(catalog[["clocks"]][[cid]]))) {
        next
      }
      cpgs <- own_scoring_cpgs(catalog[["clocks"]][[cid]], tensors)
      if (length(cpgs)) {
        catalog[["clocks"]][[cid]][["probe_sets"]] <- add_scoring_probe_set(
          catalog[["clocks"]][[cid]],
          cpgs
        )
      }
    }
  }

  # Score-assembled clocks take the union of the in-group clocks they consume,
  # repeated until nothing new resolves (in-group chains). Out-of-group inputs
  # are deliberately excluded: DNAmFitAge_{Sex} reads GrimAgeV1, but GrimAgeV1
  # is scored as its own column with its own coverage row, and folding its 1030
  # probes in here would both double-count and swamp the family's own 172/190.
  ids <- unlist(
    lapply(bundles, function(b) b[["clocks"]]),
    use.names = FALSE
  )
  repeat {
    progressed <- FALSE
    for (cid in ids) {
      entry <- catalog[["clocks"]][[cid]]
      if (length(resolved_scoring_cpgs(entry))) {
        next
      }
      # `covers` is a family-wide label (every DNAmFitAge member "covers" all
      # 14), so it describes the group, not what this clock reads. Where a
      # clock declares real inputs, those win; `covers` stays the fallback for
      # clocks that declare none.
      inputs <- as.character(entry[["depends_on_clocks"]] %||% character())
      if (!length(inputs)) {
        inputs <- as.character(entry[["covers"]] %||% character())
      }
      inputs <- setdiff(unique(inputs), cid)
      inputs <- Filter(
        function(d) {
          identical(
            catalog[["clocks"]][[d]][["group_id"]],
            entry[["group_id"]]
          )
        },
        inputs
      )
      if (!length(inputs)) {
        next
      }
      cpgs <- unique(unlist(
        lapply(inputs, function(d) {
          resolved_scoring_cpgs(catalog[["clocks"]][[d]] %||% list())
        }),
        use.names = FALSE
      ))
      if (length(cpgs)) {
        catalog[["clocks"]][[cid]][["probe_sets"]] <- add_scoring_probe_set(
          entry,
          cpgs
        )
        progressed <- TRUE
      }
    }
    if (!progressed) {
      break
    }
  }

  for (gid in names(bundles)) {
    gside <- catalog[["groups"]][[gid]]
    for (cid in bundles[[gid]][["clocks"]]) {
      entry <- catalog[["clocks"]][[cid]]
      if (length(resolved_scoring_cpgs(entry))) {
        next
      }
      cpgs <- group_shared_cpg_list(gside, bundles[[gid]][["tensors"]])
      if (length(cpgs)) {
        catalog[["clocks"]][[cid]][["probe_sets"]] <- add_scoring_probe_set(
          entry,
          cpgs
        )
      }
    }
  }
  catalog
}

# external asset encoding

# Named numeric -> double[n] in cpgs order.
align_double <- function(x, cpgs, label) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(label, ": expected a named numeric vector", call. = FALSE)
  }
  if (is.null(names(x))) {
    if (length(x) != length(cpgs)) {
      stop(label, ": unnamed vector length != n_cpgs", call. = FALSE)
    }
    return(as.double(x))
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

# Cbind named numeric tensors into a double matrix.
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

# Canonical probe order for an external pack, from the longest named vec. The pack
# stores this order once and its matrices are aligned to it, so the per-tensor CpG
# names are dropped. Bundled groups keep names on each tensor and never need this.
pack_canonical_cpgs <- function(tensors, group_id) {
  is_named_num <- vapply(
    tensors,
    function(x) {
      is.numeric(x) && !is.null(names(x)) && length(x) > 0L && is.null(dim(x))
    },
    logical(1L)
  )
  is_probe_list <- vapply(
    tensors,
    function(x) is.character(x) && length(x) > 0L && is.null(dim(x)),
    logical(1L)
  )
  named_rels <- names(tensors)[is_named_num]
  probe_list_rels <- names(tensors)[is_probe_list]
  if (!length(named_rels)) {
    stop(
      group_id,
      ": no named numeric tensors to resolve cpgs from",
      call. = FALSE
    )
  }
  lens <- vapply(named_rels, function(r) length(tensors[[r]]), integer(1L))
  ref <- names(tensors[[named_rels[which.max(lens)]]])
  for (r in named_rels[lens == max(lens)]) {
    if (!setequal(names(tensors[[r]]), ref)) {
      stop(
        group_id,
        ": probe set mismatch among cpg-aligned tensors (",
        r,
        ")",
        call. = FALSE
      )
    }
  }
  cpgs <- ref
  drop_lists <- character()
  for (r in probe_list_rels) {
    pl <- as.character(unname(tensors[[r]]))
    if (setequal(pl, ref)) {
      cpgs <- pl
      drop_lists <- c(drop_lists, r)
    }
  }
  list(cpgs = as.character(cpgs), drop_lists = unique(drop_lists))
}

# Leftover tensors that are not cpg-aligned bulk.
residual_tensors <- function(tensors, used_rels) {
  keep <- setdiff(names(tensors), used_rels)
  if (!length(keep)) {
    return(list())
  }
  tensors[keep]
}

encode_pcclocks <- function(bundle) {
  tensors <- bundle[["tensors"]]
  gid <- "PCClocks"
  resolved <- pack_canonical_cpgs(tensors, gid)
  cpgs <- resolved$cpgs

  coef_rels <- grep(
    "^weights/PCClocks/PC[^/]+\\.csv\\.gz$",
    names(tensors),
    value = TRUE
  )
  if (!length(coef_rels)) {
    stop(gid, ": no PC*.csv.gz coefficient tensors found", call. = FALSE)
  }
  col_names <- sub("\\.csv\\.gz$", "", basename(coef_rels))
  ord <- order(col_names)
  coef_rels <- coef_rels[ord]
  col_names <- col_names[ord]

  impute_rel <- "weights/PCClocks/_shared/imputeMissingCpGs.csv.gz"
  if (is.null(tensors[[impute_rel]])) {
    stop(gid, ": missing ", impute_rel, call. = FALSE)
  }

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

# 11 organ/system columns shared by organs + systems matrices.
SYSTEMSAGE_ORGANS <- c(
  "Blood",
  "Brain",
  "Heart",
  "Hormone",
  "Immune",
  "Inflammation",
  "Kidney",
  "Liver",
  "Lung",
  "Metabolic",
  "MusculoSkeletal"
)

encode_systemsage <- function(bundle) {
  tensors <- bundle[["tensors"]]
  gid <- "SystemsAge"
  resolved <- pack_canonical_cpgs(tensors, gid)
  cpgs <- resolved$cpgs

  organ_rels <- file.path(
    "weights",
    "SystemsAge",
    paste0(SYSTEMSAGE_ORGANS, ".csv.gz")
  )
  system_rels <- file.path(
    "weights",
    "SystemsAge",
    "systems",
    paste0(SYSTEMSAGE_ORGANS, ".csv.gz")
  )
  organ_rels <- gsub("\\\\", "/", organ_rels)
  system_rels <- gsub("\\\\", "/", system_rels)
  age_rel <- "weights/SystemsAge/age/age_pc_coef.csv.gz"
  impute_rel <- "weights/SystemsAge/_shared/imputeMissingCpGs.csv.gz"
  cpgs_rel <- "weights/SystemsAge/_shared/CpGs.csv.gz"

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
    cpgs_rel,
    resolved$drop_lists
  )
  bundle[["cpgs"]] <- cpgs
  bundle[["organs"]] <- cbind_aligned(
    tensors,
    organ_rels,
    cpgs,
    SYSTEMSAGE_ORGANS
  )
  bundle[["systems"]] <- cbind_aligned(
    tensors,
    system_rels,
    cpgs,
    SYSTEMSAGE_ORGANS
  )
  bundle[["age"]] <- align_double(tensors[[age_rel]], cpgs, age_rel)
  bundle[["impute"]] <- align_double(tensors[[impute_rel]], cpgs, impute_rel)
  bundle[["tensors"]] <- residual_tensors(tensors, used)
  bundle[["encoding"]] <- "canonical_matrices"
  bundle
}

encode_pcbrainage <- function(bundle) {
  tensors <- bundle[["tensors"]]
  gid <- "PCBrainAge"
  resolved <- pack_canonical_cpgs(tensors, gid)
  cpgs <- resolved$cpgs

  coef_rel <- "weights/PCBrainAge/PCBrainAge.csv.gz"
  impute_rel <- "weights/PCBrainAge/imputeMissingCpGs.csv.gz"
  for (r in c(coef_rel, impute_rel)) {
    if (is.null(tensors[[r]])) {
      stop(gid, ": missing ", r, call. = FALSE)
    }
  }

  used <- c(coef_rel, impute_rel, resolved$drop_lists)
  bundle[["cpgs"]] <- cpgs
  bundle[["coefficient_matrix"]] <- cbind_aligned(
    tensors,
    coef_rel,
    cpgs,
    "PCBrainAge"
  )
  bundle[["impute"]] <- align_double(tensors[[impute_rel]], cpgs, impute_rel)
  bundle[["tensors"]] <- residual_tensors(tensors, used)
  bundle[["encoding"]] <- "canonical_matrices"
  bundle
}

encode_external_asset <- function(bundle) {
  gid <- bundle[["group_id"]] %||% NA_character_
  if (length(bundle[["tensors"]])) {
    names(bundle[["tensors"]]) <- gsub("\\\\", "/", names(bundle[["tensors"]]))
  }
  if (identical(gid, "PCClocks")) {
    encode_pcclocks(bundle)
  } else if (identical(gid, "SystemsAge")) {
    encode_systemsage(bundle)
  } else if (identical(gid, "PCBrainAge")) {
    encode_pcbrainage(bundle)
  } else {
    stop("No external encoding for group_id=", gid, call. = FALSE)
  }
}

# Runtime registry row for mc_provenance. The content address lives in the
# filename (<group>-<payload_hash>.qs2), so it is not repeated as a field.
external_asset_registry_row <- function(a) {
  list(
    group_id = a[["group_id"]],
    # release_tag is <group>-<hash>; bare 64-hex tags are rejected by GitHub.
    release_tag = a[["release_tag"]] %||%
      sub("\\.qs2$", "", a[["file"]] %||% ""),
    file = a[["file"]],
    size_bytes = a[["size_bytes"]],
    encoding = a[["encoding"]],
    encoding_version = a[["encoding_version"]],
    n_clocks = a[["n_clocks"]],
    n_cpgs = a[["n_cpgs"]]
  )
}

# sysdata

# Flat per-clock index for list_clocks() filters.
build_index <- function(catalog) {
  clocks <- catalog[["clocks"]]
  ids <- names(clocks)

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

  idx <- data.frame(
    clock_id = ids,
    group_id = scal("group_id"),
    weights_format = scal("weights_format"),
    computation_type = scal("computation_type"),
    output_transform = scal("output_transform", "identity"),
    imputation_policy = scal("imputation_policy"),
    batch_dependent = lgl("batch_dependent"),
    external_group = lgl("external_group"),
    bib_key = scal("bib_key"),
    pmid = scal("pmid"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  idx$covariates_required <- unname(lapply(clocks, function(e) {
    as.character(e[["covariates_required"]] %||% character())
  }))
  idx$n_covariates <- lengths(idx$covariates_required)

  idx
}

# External groups keep their scoring panel in the pack (clock_scoring_cpgs()
# reads it from there), so the bundled catalog drops its copy: 92% of all panel
# elements, and the only other place they could drift from the weights.
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
  catalog[["clocks"]] <- attach_custom_components(catalog[["clocks"]])
  bundles <- build_group_bundles(repo_path, catalog, ship_groups)
  catalog <- resolve_group_scoring_probe_sets(catalog, bundles)
  catalog <- attach_sex_routed_aliases(catalog)
  catalog[["clocks"]] <- trim_build_only_fields(catalog[["clocks"]])

  mc_catalog <- drop_external_probe_cpgs(catalog[["clocks"]])
  mc_groups <- catalog[["groups"]]
  mc_bundles <- bundles
  mc_index <- build_index(catalog)
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
    external_groups = EXTERNAL_GROUPS,
    external_assets = ext_reg
  )

  usethis::use_data(
    mc_catalog,
    mc_groups,
    mc_bundles,
    mc_index,
    mc_provenance,
    internal = TRUE,
    overwrite = TRUE,
    compress = "gzip"
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
      "mc_provenance"
    )
  ))
}

# content-addressed external packs

# Low ZSTD, no shuffle.
QS2_COMPRESS_LEVEL <- 1L
QS2_SHUFFLE <- FALSE

# Canonical pack for hash + qs_save. Carries weights only: the catalog and group
# metadata live in sysdata (mc_catalog/mc_groups), which is what the runtime reads.
stable_external_payload <- function(bundle) {
  for (f in EXTERNAL_PIN_FIELDS) {
    bundle[[f]] <- NULL
  }

  tensors <- bundle[["tensors"]] %||% list()
  if (length(tensors) && !is.null(names(tensors))) {
    tensors <- tensors[sort(names(tensors))]
  }

  clocks <- as.character(bundle[["clocks"]] %||% character())
  if (length(clocks)) {
    clocks <- sort(unique(clocks))
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

# Stable serialize (version=2L, xdr=TRUE).
payload_hash_of <- function(payload) {
  digest::digest(
    serialize(payload, connection = NULL, version = 2L, xdr = TRUE),
    algo = "sha256",
    serialize = FALSE
  )
}

# GitHub release target

parse_github_owner_repo <- function(url) {
  url <- trimws(as.character(url %||% ""))
  if (!nzchar(url)) {
    return(NULL)
  }
  url <- sub("\\.git$", "", url)
  m <- regexec("(?:github\\.com[:/])([^/]+)/([^/]+)$", url, perl = TRUE)
  parts <- regmatches(url, m)[[1L]]
  if (length(parts) != 3L) {
    return(NULL)
  }
  list(
    owner = parts[[2L]],
    repo = parts[[3L]],
    slug = paste0(parts[[2L]], "/", parts[[3L]])
  )
}

# Release target: env MC_RELEASE_REPO, else package origin remote.
package_release_repo <- function() {
  env <- Sys.getenv("MC_RELEASE_REPO", unset = "")
  if (nzchar(env)) {
    if (grepl("/", env) && !grepl("github\\.com", env)) {
      return(list(
        owner = sub("/.*", "", env),
        repo = sub(".*/", "", env),
        slug = env
      ))
    }
    parsed <- parse_github_owner_repo(env)
    if (!is.null(parsed)) {
      return(parsed)
    }
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

# Upload via `uv run python data-raw/gh_upload.py` (idempotent by release_tag).
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

# Prefer MC_UPLOAD_PAT over GITHUB_PAT.
upload_pat <- function() {
  for (v in c("MC_UPLOAD_PAT", "GITHUB_TOKEN", "GH_TOKEN")) {
    pat <- Sys.getenv(v, unset = "")
    if (nzchar(pat)) return(pat)
  }
  ""
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
      "(a fine-grained PAT with Contents:read/write on the package repo). ",
      "Set it in ~/.Renviron, and keep it OUT of GITHUB_PAT/GITHUB_TOKEN so it ",
      "does not shadow the broad token remotes::install_github() reads. ",
      "(GITHUB_TOKEN/GH_TOKEN remain honored as a fallback.)",
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
  # Feed PAT only to the child env; stage JSON for processx stdin.
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

  res <- tryCatch(
    jsonlite::fromJSON(proc$stdout, simplifyVector = FALSE),
    error = function(e) NULL
  )
  for (r in res$results %||% list()) {
    message(
      "sync: ",
      r$action,
      " ",
      r$name,
      " -> ",
      repo$slug,
      " @ tag ",
      r$tag
    )
  }
  invisible(assets)
}

# Build content-addressed external packs: <group>-<payload_hash>.qs2.
build_external_assets <- function(repo_path, catalog, external_groups) {
  fs::dir_create(asset_dir)
  assets <- list()

  for (gid in external_groups) {
    message("sync: building external asset for ", gid, "...")
    raw_bundle <- build_group_bundles(repo_path, catalog, gid)[[gid]]
    # Resolve probe sets before encoding; adopt catalog so pack and sysdata match.
    catalog <- resolve_group_scoring_probe_sets(
      catalog,
      stats::setNames(list(raw_bundle), gid)
    )
    bundle <- encode_external_asset(raw_bundle)
    bundle[["schema_version"]] <- catalog[["schema_version"]]
    bundle[["encoding_version"]] <- EXTERNAL_ENCODING_VERSION

    payload <- stable_external_payload(bundle)
    phash <- payload_hash_of(payload)
    fname <- sprintf("%s-%s.qs2", tolower(gid), phash)
    fpath <- file.path(asset_dir, fname)
    # Tag = filename stem (<group>-<hash>); bare hex tags are rejected by GitHub.
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

  # Prune superseded local staging packs; published releases are untouched.
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

# Hit when source sha matches and every recorded pack exists on disk.
lockfile_hit <- function(lock, current_sha) {
  if (
    is.null(lock) ||
      !identical(as.character(lock$source_git_sha %||% ""), current_sha)
  ) {
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

# Store pre-trim external catalog entries for lockfile reuse.
write_lockfile <- function(sha, assets, ext_clocks) {
  saveRDS(
    list(source_git_sha = sha, assets = assets, ext_clocks = ext_clocks),
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

  gids <- unique(vapply(
    catalog[["clocks"]],
    function(c) c[["group_id"]],
    character(1L)
  ))
  external <- sort(intersect(gids, EXTERNAL_GROUPS))
  ship <- sort(setdiff(gids, EXTERNAL_GROUPS))
  message(
    "sync: ",
    catalog[["n_clocks"]],
    " clocks; ship groups=",
    length(ship),
    "; external=",
    paste(external, collapse = ", ")
  )

  lock <- if (isTRUE(force)) NULL else read_lockfile()
  if (lockfile_hit(lock, current_sha)) {
    message(
      "sync: external assets unchanged (git sha match) -- reusing ",
      length(lock$assets),
      " cached pack(s), skipping rebuild. Run sync(force = TRUE) to reconcile drift."
    )
    assets <- lock$assets
    # Restore external resolved probe sets from the lockfile.
    for (cid in names(lock$ext_clocks)) {
      catalog[["clocks"]][[cid]] <- lock$ext_clocks[[cid]]
    }
  } else {
    ext <- build_external_assets(src$path, catalog, external)
    assets <- ext$assets
    # Adopt catalog with resolved external probe sets.
    catalog <- ext$catalog
    ext_ids <- names(catalog[["clocks"]])[
      vapply(
        catalog[["clocks"]],
        function(e) as.character(e[["group_id"]] %||% "") %in% external,
        logical(1L)
      )
    ]
    write_lockfile(current_sha, assets, catalog[["clocks"]][ext_ids])
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
