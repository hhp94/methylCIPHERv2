# Descriptive catalog verbs: codebook() and bibliography(). Plain functions (not
# S3 methods) -- they accept an mc_result, clock ids / groups / tags / "all", and
# describe those clocks from the bundled catalog.

# clock ids to describe: a result's returned clocks, or resolve the tokens
codebook_ids <- function(x) {
  if (inherits(x, "mc_result")) {
    return(x[["provenance"]][["clocks"]])
  }
  resolve_clocks(x)
}

# CpG count for the codebook. A sex-routed alias has no scoring panel of its own,
# so report the union of its routed members' panels rather than a misleading 0.
codebook_n_cpgs <- function(id) {
  if (identical(clock_kind(id), "sex_routed_alias")) {
    members <- unlist(clock_routing(id), use.names = FALSE)
    cpgs <- unique(unlist(lapply(
      members,
      function(m) tryCatch(clock_scoring_cpgs(m), error = function(e) character(0))
    )))
    return(length(cpgs))
  }
  tryCatch(length(clock_scoring_cpgs(id)), error = function(e) NA_integer_)
}

# one descriptive row per clock: what it is, what it needs, and its reference.
# Presents the scoring-contract metadata the catalog carries; training-population
# / phenotype fields are not synced into the package yet (see dev notes).
#' @export
codebook <- function(x = "all") {
  ids <- codebook_ids(x)
  idx <- match(ids, mc_index[["clock_id"]])

  n_cpgs <- vapply(ids, codebook_n_cpgs, integer(1L))
  covariates <- vapply(
    ids,
    function(id) {
      cv <- clock_covariates_required(id)
      if (is.null(cv) || !length(cv)) "" else paste(cv, collapse = ", ")
    },
    character(1L)
  )

  data.frame(
    clock_id = ids,
    group_id = mc_index[["group_id"]][idx],
    n_cpgs = unname(n_cpgs),
    computation = mc_index[["computation_type"]][idx],
    covariates = unname(covariates),
    external = mc_index[["external_group"]][idx],
    batch_dependent = mc_index[["batch_dependent"]][idx],
    pmid = mc_index[["pmid"]][idx],
    reference = mc_index[["bib_key"]][idx],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# short "Author et al. Year" from a bib_key of the form Author_Year_PMID
bib_citation <- function(key) {
  vapply(
    key,
    function(k) {
      parts <- strsplit(k, "_", fixed = TRUE)[[1L]]
      if (length(parts) < 2L) {
        return(k)
      }
      sprintf("%s et al. %s", parts[[1L]], parts[[2L]])
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

pubmed_url <- function(pmid) {
  pmid <- as.character(pmid)
  ifelse(
    !is.na(pmid) & nzchar(pmid) & pmid != "NA",
    paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid),
    NA_character_
  )
}

# read the shipped clocks.bib into a named list: cite key -> raw BibTeX entry.
# The file is sync-vendored to inst/bibliography/ (data-raw/sync.R,
# vendor_bibliography()). NULL if absent (older installs).
mc_read_bib <- function() {
  path <- system.file("bibliography", "clocks.bib", package = "methylCIPHERv2")
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^@", lines)
  if (!length(starts)) {
    return(NULL)
  }
  ends <- c(starts[-1L] - 1L, length(lines))
  out <- lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    sub("\\s+$", "", paste(block, collapse = "\n"))
  })
  names(out) <- sub("^@\\w+\\{([^,]+),.*$", "\\1", lines[starts])
  out
}

# one {field} value from a raw BibTeX entry (tolerates single-level nesting)
bib_field <- function(raw, field) {
  m <- regmatches(
    raw,
    regexec(paste0(field, "\\s*=\\s*\\{((?:[^{}]|\\{[^{}]*\\})*)\\}"), raw, perl = TRUE)
  )[[1L]]
  if (length(m) < 2L) NA_character_ else trimws(gsub("[{}]", "", m[[2L]]))
}

bib_first_author <- function(raw) {
  a <- bib_field(raw, "author")
  if (is.na(a)) {
    return(NA_character_)
  }
  trimws(sub(",.*$", "", strsplit(a, " and ", fixed = TRUE)[[1L]][1L]))
}

# references for the requested clocks, one row per unique publication, enriched
# from the shipped clocks.bib. format = "bibtex" prints the full BibTeX entries
# (a stub when a key is missing from the .bib).
#' @export
bibliography <- function(x = "all", format = c("data.frame", "bibtex")) {
  format <- match.arg(format)
  ids <- codebook_ids(x)
  idx <- match(ids, mc_index[["clock_id"]])

  refs <- data.frame(
    reference = mc_index[["bib_key"]][idx],
    pmid = as.character(mc_index[["pmid"]][idx]),
    stringsAsFactors = FALSE
  )
  refs <- unique(refs[nzchar(refs[["reference"]]), , drop = FALSE])
  refs <- refs[order(refs[["reference"]]), , drop = FALSE]
  bib <- mc_read_bib()

  if (identical(format, "bibtex")) {
    entries <- vapply(
      seq_len(nrow(refs)),
      function(i) {
        k <- refs[["reference"]][i]
        if (!is.null(bib) && !is.null(bib[[k]])) {
          bib[[k]]
        } else {
          sprintf(
            "@article{%s,\n  pmid = {%s},\n  url = {%s}\n}",
            k, refs[["pmid"]][i], pubmed_url(refs[["pmid"]][i])
          )
        }
      },
      character(1L)
    )
    cat(paste(entries, collapse = "\n\n"), "\n")
    return(invisible(refs[["reference"]]))
  }

  from_bib <- function(fn) {
    vapply(
      refs[["reference"]],
      function(k) if (is.null(bib) || is.null(bib[[k]])) NA_character_ else fn(bib[[k]]),
      character(1L),
      USE.NAMES = FALSE
    )
  }
  author1 <- from_bib(bib_first_author)
  year <- sub("^.*_(\\d{4})_.*$", "\\1", refs[["reference"]])
  data.frame(
    reference = refs[["reference"]],
    citation = ifelse(
      is.na(author1),
      bib_citation(refs[["reference"]]),
      sprintf("%s et al. %s", author1, year)
    ),
    title = from_bib(function(r) bib_field(r, "title")),
    journal = from_bib(function(r) bib_field(r, "journal")),
    year = year,
    doi = from_bib(function(r) bib_field(r, "doi")),
    pmid = refs[["pmid"]],
    url = pubmed_url(refs[["pmid"]]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
