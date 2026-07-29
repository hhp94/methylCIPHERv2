# shared cross-cutting helpers

# name a character vector as cli "*" bullets
bullets <- function(x) {
  stats::setNames(x, rep("*", length(x)))
}

# cli "*" bullets, capped at 10 with a "... and n more" tail
capped_bullets <- function(lines) {
  shown <- utils::head(lines, 10L)
  if (length(lines) > length(shown)) {
    shown <- c(shown, sprintf("... and %d more", length(lines) - length(shown)))
  }
  bullets(shown)
}

# polynomial eval, lowest degree first (horner-style, 1-row safe)
poly_eval <- function(x, coef) {
  out <- rep(0, length(x))
  for (k in seq_along(coef)) {
    out <- out + coef[[k]] * x^(k - 1L)
  }
  out
}
