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

# one descriptive row per clock: what it is, what it needs, and its reference.
# Presents the scoring-contract metadata the catalog carries; training-population
# / phenotype fields are not synced into the package yet (see dev notes).
#' @export
codebook <- function(x = "all") {
  ids <- codebook_ids(x)
  idx <- match(ids, mc_index[["clock_id"]])

  n_cpgs <- vapply(
    ids,
    function(id) tryCatch(length(clock_scoring_cpgs(id)), error = function(e) NA_integer_),
    integer(1L)
  )
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

# references for the requested clocks, one row per unique publication. Full
# BibTeX (titles/authors) needs the upstream clocks.bib, which is not bundled;
# until then this emits the citation key, PMID, and a PubMed URL.
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
  refs <- refs[nzchar(refs[["reference"]]), , drop = FALSE]
  refs <- unique(refs)
  refs <- refs[order(refs[["reference"]]), , drop = FALSE]
  refs[["citation"]] <- bib_citation(refs[["reference"]])
  refs[["url"]] <- ifelse(
    nzchar(refs[["pmid"]]) & refs[["pmid"]] != "NA",
    paste0("https://pubmed.ncbi.nlm.nih.gov/", refs[["pmid"]]),
    NA_character_
  )
  refs <- refs[, c("reference", "citation", "pmid", "url")]
  rownames(refs) <- NULL

  if (identical(format, "data.frame")) {
    return(refs)
  }

  # minimal but valid BibTeX stubs (key + pmid + url)
  entries <- vapply(
    seq_len(nrow(refs)),
    function(i) {
      sprintf(
        "@article{%s,\n  pmid = {%s},\n  url = {%s}\n}",
        refs[["reference"]][i], refs[["pmid"]][i], refs[["url"]][i]
      )
    },
    character(1L)
  )
  cat(paste(entries, collapse = "\n\n"), "\n")
  invisible(refs)
}
