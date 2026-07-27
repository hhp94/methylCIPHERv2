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
