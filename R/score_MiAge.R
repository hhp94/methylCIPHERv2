# MiAge: multi-start L-BFGS-B mitotic age over n in [10, 10000]

MIAGE_LOWER <- 10
MIAGE_UPPER <- 10000

# four interior starts plus author default 500
MIAGE_STARTS <- c(
  MIAGE_LOWER + seq_len(4L) * (MIAGE_UPPER - MIAGE_LOWER) / 5,
  500
)

# best multi-start fit for one sample.
#
# `log(b)` does not depend on n, and L-BFGS-B asks for the objective and the
# gradient at the same n -- so hoist the log out of the search and bind
# `b^(n - 1)` and the residual once per n, letting the second call read them.
# The arithmetic is left exactly as written (`d` is not folded into `log(b)`),
# so both values are bit-identical to recomputing them.
miage_fit <- function(betaj, b, c, d) {
  logb <- log(b)
  at <- NA_real_
  bn <- NULL
  res <- NULL
  bind <- function(n) {
    if (!isTRUE(n == at)) {
      bn <<- b^(n - 1)
      res <<- c + bn * d - betaj
      at <<- n
    }
  }

  objective <- function(n) {
    bind(n)
    sum(res^2)
  }
  gradient <- function(n) {
    bind(n)
    2 * sum(res * bn * logb * d)
  }

  fits <- lapply(MIAGE_STARTS, function(start) {
    stats::optim(
      par = start,
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      lower = MIAGE_LOWER,
      upper = MIAGE_UPPER,
      control = list(factr = 1)
    )
  })
  fits[[which.min(vapply(fits, function(f) f[["value"]], numeric(1)))]][["par"]]
}

score_MiAge <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)

  params <- miage_params(id)
  present <- cpgs[["score_present"]]
  cache <- block[["partial_cache"]]
  cached <- cached_cols(present, cache)
  betas <- block[["DNAm"]][, present, drop = FALSE]
  if (length(cached)) {
    betas[, cached] <- cache[, cached]
  }

  b <- params[["b"]][present]
  cc <- params[["c"]][present]
  d <- params[["d"]][present]
  score_matrix(
    vapply(seq_len(n), function(i) miage_fit(betas[i, ], b, cc, d), numeric(1)),
    sample_id,
    id
  )
}

# MiAge site-specific params: named b, c, d vectors in panel order
miage_params <- function(id) {
  tab <- component_tensor(id, "cpg")
  lapply(
    list(b = tab[["b"]], c = tab[["c"]], d = tab[["d"]]),
    function(x) stats::setNames(as.numeric(x), tab[["cpg"]])
  )
}
