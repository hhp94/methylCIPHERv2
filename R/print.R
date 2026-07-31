# shared print grammar for every print.mc_* method. builders return strings so cli and cat match

# one plural form for every count we print. suffix covers nouns like batch(es)
plural_count <- function(n, noun, suffix = "s") {
  sprintf("%d %s(%s)", n, noun, suffix)
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

# one component block. cut_cols = false when columns stay whole (pheno)
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
