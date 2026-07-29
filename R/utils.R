# shared cross-cutting helpers

# name a character vector as cli "*" bullets
bullets <- function(x) {
  stats::setNames(x, rep("*", length(x)))
}

# polynomial eval, lowest degree first (horner-style, 1-row safe)
poly_eval <- function(x, coef) {
  out <- rep(0, length(x))
  for (k in seq_along(coef)) {
    out <- out + coef[[k]] * x^(k - 1L)
  }
  out
}
