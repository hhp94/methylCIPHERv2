# result-record methods for class "mc_result"

# compact view of a scored result: dims, clock ids, a corner of the scores
#' @export
print.mc_result <- function(x, n = 6L, k = 4L, ...) {
  s <- x[["scores"]]
  cat(sprintf("<mc_result> %d sample(s) x %d clock(s)\n", nrow(s), ncol(s)))
  ids <- colnames(s)
  shown <- utils::head(ids, 8L)
  cat(sprintf(
    "clocks: %s%s\n\n",
    paste(shown, collapse = ", "),
    if (length(ids) > length(shown)) sprintf(", ... (%d total)", length(ids)) else ""
  ))
  nr <- min(n, nrow(s))
  nc <- min(k, ncol(s))
  print(round(s[seq_len(nr), seq_len(nc), drop = FALSE], 3L))
  if (nr < nrow(s) || nc < ncol(s)) {
    cat(sprintf("... %d more row(s), %d more col(s)\n", nrow(s) - nr, ncol(s) - nc))
  }
  cat("\nExtract scores with $scores; per-clock QC with clocks_coverage().\n")
  invisible(x)
}

# scores as a plain n x k double matrix
#' @export
as.matrix.mc_result <- function(x, ...) {
  x[["scores"]]
}

# scores as a data.frame: the id column plus one column per clock. No pheno --
# that stays on $pheno / augment(), so it cannot leak into an analysis by accident.
#' @export
as.data.frame.mc_result <- function(x, ...) {
  s <- x[["scores"]]
  id <- x[["provenance"]][["pheno_id"]]
  out <- data.frame(
    .id = rownames(s),
    s,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  names(out)[1L] <- id
  out
}

# analysis-ready table: scores + id + the aligned covariates the record carries
# (Age/Female/...), plus an optional user data.frame joined by the id column.
# This is the covariate-appended view for downstream modelling -- e.g.
# lm(Horvath1 ~ Age + Female, data = augment(res)).
#' @export
augment <- function(x, data = NULL, ...) {
  check_mc_result(x)
  id <- x[["provenance"]][["pheno_id"]]
  out <- as.data.frame(x)

  pheno <- x[["pheno"]]
  if (!is.null(pheno) && ncol(pheno) > 1L) {
    out <- merge(out, pheno, by = id, all.x = TRUE, sort = FALSE)
  }
  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (!id %in% names(data)) {
      cli::cli_abort(
        c(
          "{.arg data} must have the id column {.field {id}} to join on.",
          "i" = "It is joined to the scores by sample id."
        ),
        call = NULL
      )
    }
    out <- merge(out, data, by = id, all.x = TRUE, sort = FALSE)
  }
  out
}

# rbind is refused -- re-run calc_clocks() on the combined DNAm instead
#' @export
rbind.mc_result <- function(...) {
  cli::cli_abort(
    c(
      "{.cls mc_result} records cannot be {.fn rbind}-ed.",
      "i" = "Re-run {.fn calc_clocks} on the combined DNAm -- batch-dependent
             clocks must see all samples at once."
    ),
    call = NULL
  )
}
