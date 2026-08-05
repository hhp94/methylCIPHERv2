# miAge: multi-start L-BFGS-B mitotic age over n in [10, 10000]

MIAGE_BOUNDS <- c(10, 10000)

# four interior starts plus author default 500
MIAGE_STARTS <- c(
  MIAGE_BOUNDS[[1L]] + seq_len(4L) * diff(MIAGE_BOUNDS) / 5,
  500
)

miage_fit <- function(betaj, b, cc, d) {
  logb <- log(b)
  st <- new.env(parent = emptyenv())
  st[["at"]] <- NA_real_
  bind <- function(n) {
    if (!isTRUE(n == st[["at"]])) {
      bn <- b^(n - 1)
      st[["bn"]] <- bn
      st[["res"]] <- cc + bn * d - betaj
      st[["at"]] <- n
    }
    st
  }

  objective <- function(n) {
    sum(bind(n)[["res"]]^2)
  }
  gradient <- function(n) {
    s <- bind(n)
    2 * sum(s[["res"]] * s[["bn"]] * logb * d)
  }

  fits <- lapply(MIAGE_STARTS, function(start) {
    stats::optim(
      par = start,
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      lower = MIAGE_BOUNDS[[1L]],
      upper = MIAGE_BOUNDS[[2L]],
      control = list(factr = 1)
    )
  })
  fits[[which.min(vapply(fits, function(f) f[["value"]], numeric(1)))]][["par"]]
}

score_MiAge <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  params <- miage_params(id)
  obs <- observed_panel(
    cpgs[["score_present"]],
    cpgs[["score_present_idx"]],
    block
  )
  betas <- obs[["values"]]
  panel <- obs[["cols"]]

  b <- params[["b"]][panel]
  cc <- params[["c"]][panel]
  d <- params[["d"]][panel]
  score_matrix(
    vapply(
      seq_along(sample_id),
      function(i) miage_fit(betas[i, ], b, cc, d),
      numeric(1)
    ),
    sample_id,
    id
  )
}

# miAge site-specific params: named b, c, d vectors in panel order
miage_params <- function(id) {
  tab <- component_tensor(id, "cpg")
  lapply(
    list(b = tab[["b"]], c = tab[["c"]], d = tab[["d"]]),
    function(x) stats::setNames(as.numeric(x), tab[["cpg"]])
  )
}
