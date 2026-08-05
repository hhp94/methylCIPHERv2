# citations for scored clocks (mc_citations join + clocks.bib text)

#' Clock Citations
#'
#' Builds the citations for a set of clocks, or for the clocks scored in an
#' `mc_result` object.
#'
#' @param x A character vector. The clock ids or group ids to cite, or an
#'   `mc_result` object from [calc_clocks()].
#' @param ... Not used.
#'
#' @details
#' A character vector cites the clocks it names, and a group id cites every
#' clock in the group. An `mc_result` object cites the clocks in its scores.
#' Any other class raises an error that names the two accepted types.
#'
#' @returns An `mc_citation` object. It holds the clock-to-paper links and
#'   the bibtex text for each paper.
#'
#' @examples
#' cite_clocks(c("Horvath1", "Hannum"))
#'
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 5)
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#' cite_clocks(res)
#'
#' @export
cite_clocks <- function(x, ...) {
  UseMethod("cite_clocks")
}

#' @rdname cite_clocks
#' @export
cite_clocks.character <- function(x, ...) {
  new_mc_citation(resolve_clocks(x))
}

# score columns only -- cite_also is declared upstream, never walked
#' @rdname cite_clocks
#' @export
cite_clocks.mc_result <- function(x, ...) {
  new_mc_citation(colnames(x[["scores"]]))
}

#' @rdname cite_clocks
#' @export
cite_clocks.default <- function(x, ...) {
  cli::cli_abort(
    c(
      "{.obj_type_friendly {x}} has no {.fn cite_clocks} method.",
      "i" = "Pass a character vector of clock ids or group ids.",
      "i" = "Or pass an {.cls mc_result} from {.fn calc_clocks}."
    ),
    call = NULL
  )
}

# header + bibtex in the shared printer grammar, via cli_verbatim so alignment holds
#' Print Method For An mc_citation Object
#'
#' Prints the clock and paper counts, then the bibtex text, for an
#' `mc_citation` object.
#'
#' @param x An `mc_citation` object. The value returned by [cite_clocks()].
#' @param ... Not used.
#'
#' @returns An `mc_citation` object. Returns `x`, invisibly, after printing
#'   it.
#'
#' @examples
#' cite_clocks(c("Horvath1", "Hannum"))
#'
#' @export
print.mc_citation <- function(x, ...) {
  links <- x[["links"]]
  n_clocks <- length(unique(links[["clock_id"]]))
  n_papers <- length(unique(links[["bib_key"]]))

  # cli_verbatim drops an empty string, so blank lines are cat()
  cli::cli_verbatim(
    fmt_header("mc_citation", n_clocks, "clock", n_papers, "paper")
  )
  cat("\n")
  cli::cli_verbatim(fmt_section("bibtex", plural_count(n_papers, "paper")))
  # bibliography is never capped. writeLines, not a cli template.
  writeLines(x[["bibtex"]])
  cat("\n")
  cli::cli_bullets(c(
    "i" = "{.code as.data.frame(x)} returns the clock-to-paper table.",
    "i" = "{.code writeLines(toBibtex(x), \"refs.bib\")} writes the bibtex to
           a file.",
    "i" = "{.code citation(\"methylCIPHERv2\")} cites the package itself."
  ))
  invisible(x)
}

#' Data Frame Method For An mc_citation Object
#'
#' Converts the [cite_clocks()] output to a data.frame.
#'
#' @param x An `mc_citation` object. The value returned by [cite_clocks()].
#' @param row.names A character vector. Not used by this method. Default is
#'   `NULL`.
#' @param optional A boolean. Not used by this method. Default is `FALSE`.
#' @param ... Not used.
#'
#' @returns A data.frame. One row for each clock, with the bib key of the
#'   paper it cites and publication details such as the title, the authors,
#'   and the DOI.
#'
#' @examples
#' cites <- cite_clocks(c("Horvath1", "Hannum"))
#' as.data.frame(cites)
#'
#' @export
as.data.frame.mc_citation <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  x[["links"]]
}

#' Bibtex Method For An mc_citation Object
#'
#' Returns the bibtex text for the papers cited in an `mc_citation` object.
#'
#' @param object An `mc_citation` object. The value returned by
#'   [cite_clocks()].
#' @param ... Not used.
#'
#' @returns A character vector of class `Bibtex`. One bibtex entry for each
#'   cited paper, as lines of text.
#'
#' @examples
#' cites <- cite_clocks(c("Horvath1", "Hannum"))
#' toBibtex(cites)
#'
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
          "Bib key(s) missing from the vendored bibliography: %s. ",
          "The catalog and clocks.bib are out of step -- please report it."
        ),
        capped(absent)
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
            "No citation on record for %s. ",
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
