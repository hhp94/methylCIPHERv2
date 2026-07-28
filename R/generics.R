# result-record methods for class "mc_result"

# print mc_result: every component is labelled with the name it is reached
# by, so the record reads as the list it is
#' @export
print.mc_result <- function(x, n = 6, p = 6, ...) {
  scores <- x[["scores"]]
  pheno <- x[["pheno"]]
  nr <- nrow(scores)
  nc <- ncol(scores)

  cat(sprintf("<mc_result> %d sample(s) x %d clock(s)\n", nr, nc))

  cat("\n<$pheno>")
  if (is.null(pheno)) {
    cat(" NULL\n")
  } else {
    pn <- min(n, nrow(pheno))
    cat(sprintf(" [%d of %d row(s)]\n", pn, nrow(pheno)))
    print(pheno[seq_len(pn), , drop = FALSE])
    if (pn < nrow(pheno)) {
      cat("...\n")
    }
  }

  ni <- min(n, nr)
  pi <- min(p, nc)
  cat(sprintf(
    "\n<$scores> [%d of %d row(s), %d of %d clock(s)]\n",
    ni,
    nr,
    pi,
    nc
  ))
  print(scores[seq_len(ni), seq_len(pi), drop = FALSE])
  if (ni < nr || pi < nc) {
    cat("...\n")
  }

  invisible(x)
}

# naked scores; coverage and provenance stay on the record
#' @export
as.matrix.mc_result <- function(x, ...) {
  x[["scores"]]
}

# rbind is refused -- re-run calc_clocks() on the combined DNAm instead
#' @export
rbind.mc_result <- function(..., deparse.level = 1) {
  cli::cli_abort(
    c(
      "{.cls mc_result} records cannot be {.fn rbind}-ed.",
      "i" = "Re-run {.fn calc_clocks} on the combined DNAm -- batch-dependent
             clocks must see all samples at once."
    ),
    call = NULL
  )
}
