# citations for scored clocks (mc_citations join + clocks.bib text)

#' @export
cite_clocks <- function(x, ...) {
  UseMethod("cite_clocks")
}

#' @export
cite_clocks.character <- function(x, ...) {
  new_mc_citation(resolve_clocks(x))
}

# score columns only -- cite_also is declared upstream, never walked
#' @export
cite_clocks.mc_result <- function(x, ...) {
  new_mc_citation(colnames(x[["scores"]]))
}

#' @export
cite_clocks.default <- function(x, ...) {
  cli::cli_abort(
    c(
      "Can't cite {.obj_type_friendly {x}}.",
      "i" = "Pass clock or group ids, or an {.cls mc_result}."
    ),
    call = NULL
  )
}

#' @export
print.mc_citation <- function(x, ...) {
  links <- x[["links"]]
  n_clocks <- length(unique(links[["clock_id"]]))
  n_papers <- length(unique(links[["bib_key"]]))
  cli::cli_text("{n_clocks} clock{?s}, {n_papers} paper{?s}.")
  cat("\n")
  # bibtex is pre-aligned (cli reflows non-verbatim)
  cli::cli_verbatim(x[["bibtex"]])
  cat("\n")
  cli::cli_alert_info(
    "Export the citations with {.code writeLines(toBibtex(x), \"refs.bib\")}.
     To cite the package itself see {.code citation(\"methylCIPHERv2\")}."
  )
  invisible(x)
}

#' @export
as.data.frame.mc_citation <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  x[["links"]]
}

#' @export
toBibtex.mc_citation <- function(object, ...) {
  structure(object[["bibtex"]], class = "Bibtex")
}

# record: links (clock -> paper) + bibtex text
new_mc_citation <- function(ids) {
  links <- citation_links(ids)
  entries <- bib_entries()
  keys <- unique(links[["bib_key"]])

  absent <- setdiff(keys, names(entries))
  if (length(absent)) {
    stop(
      sprintf(
        paste0(
          "%d bib key(s) missing from the vendored bibliography: %s. ",
          "The catalog and clocks.bib are out of step -- please report it."
        ),
        length(absent),
        paste(absent, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  bibtex <- unlist(
    lapply(keys, function(k) c(entries[[k]], "")),
    use.names = FALSE
  )
  structure(
    list(links = links, bibtex = utils::head(bibtex, -1L)),
    class = "mc_citation"
  )
}

# citation rows per clock, in request order, primary first
citation_links <- function(ids) {
  rows <- lapply(ids, function(id) {
    # sex-routed alias cites through its donor
    key <- as.character(clock_entry(id)[["donor_clock_id"]] %||% id)
    hit <- mc_citations[mc_citations[["clock_id"]] == key, , drop = FALSE]
    if (!nrow(hit)) {
      stop(
        sprintf(
          paste0(
            "%s has no citation on record. ",
            "The catalog is out of step with its bibliography -- please ",
            "report it."
          ),
          id
        ),
        call. = FALSE
      )
    }
    hit[["clock_id"]] <- id
    hit
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

# bib_key -> entry lines from the vendored .bib
bib_entries <- function() {
  path <- system.file(
    "bibliography",
    "clocks.bib",
    package = "methylCIPHERv2"
  )
  if (!nzchar(path)) {
    stop(
      "The vendored clocks.bib is missing from the package.",
      call. = FALSE
    )
  }
  con <- file(path, encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  starts <- grep("^@", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  entries <- Map(
    function(s, e) {
      block <- lines[seq.int(s, e)]
      block[seq_len(max(which(nzchar(trimws(block)))))]
    },
    starts,
    ends
  )
  stats::setNames(
    entries,
    trimws(sub("^@[^{]+\\{([^,]+),.*$", "\\1", lines[starts]))
  )
}
