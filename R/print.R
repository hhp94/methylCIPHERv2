# one visual grammar for every print.mc_* method. these records are lists, so
# each prints as a "<class> A x B" header, then one "$component [what is shown]"
# section per element, then a "... N more" tail on the axes that were cut.
# builders return strings, so a cli printer emits the same text verbatim.

# "3 clock(s)" -- one plural form for every count we print
plural_count <- function(n, noun) {
  sprintf("%d %s(s)", n, noun)
}

# "6 of 10 row(s)" -- the "of" form means the axis can be cut
shown_count <- function(i, n, noun) {
  sprintf("%d of %s", i, plural_count(n, noun))
}

# "4 more row(s)", or nothing when the axis is whole
more_count <- function(i, n, noun) {
  if (i >= n) {
    return(character(0))
  }
  sprintf("%d more %s(s)", n - i, noun)
}

# "<mc_result> 10 sample(s) x 3 clock(s)"
fmt_header <- function(cls, n, n_noun, k, k_noun) {
  sprintf(
    "<%s> %s x %s",
    cls,
    plural_count(n, n_noun),
    plural_count(k, k_noun)
  )
}

# "$scores [6 of 10 row(s), 3 of 3 clock(s)]"
fmt_section <- function(name, ...) {
  sprintf("$%s [%s]", name, paste(c(...), collapse = ", "))
}

# one component: its section line, an ni x pi head, then the axes it cut.
# cut_cols = FALSE for a component we never narrow (pheno keeps its columns)
print_block <- function(name, x, ni, pi, col_noun, cut_cols = TRUE) {
  nr <- nrow(x)
  nc <- ncol(x)
  cols <- if (cut_cols) {
    shown_count(pi, nc, col_noun)
  } else {
    plural_count(nc, col_noun)
  }
  cat("\n", fmt_section(name, shown_count(ni, nr, "row"), cols), "\n", sep = "")
  print(x[seq_len(ni), seq_len(pi), drop = FALSE])

  tail <- c(more_count(ni, nr, "row"), more_count(pi, nc, col_noun))
  if (length(tail)) {
    cat("... ", paste(tail, collapse = ", "), "\n", sep = "")
  }
  invisible(NULL)
}
