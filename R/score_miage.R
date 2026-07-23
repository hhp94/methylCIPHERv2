# MiAge (mitotic age): per-sample lifetime cell divisions
# n = argmin over n in [10, 10000] of sum_i (c_i + b_i^(n-1) * d_i - beta_i)^2,
# with site-specific (b, c, d) fixed from training. No closed form -- the author
# minimizes with bounded L-BFGS-B from 5 starts and keeps the best fit.
# Absent CpGs drop out of the sum (policy "omit"), as they do in the author code.

MIAGE_LOWER <- 10
MIAGE_UPPER <- 10000

# Author start grid: 4 evenly spaced interior points, then the 500 default last.
MIAGE_STARTS <- c(
  MIAGE_LOWER + seq_len(4L) * (MIAGE_UPPER - MIAGE_LOWER) / 5,
  500
)

# Divisions for one sample: best of the multi-start fits, earliest start on ties.
# betaj is fully observed (the engine drops absent CpGs and fills partial NA), so
# the author's na.rm on the residual sum has nothing left to do here.
miage_fit <- function(betaj, b, c, d) {
  objective <- function(n) sum((c + b^(n - 1) * d - betaj)^2)
  gradient <- function(n) {
    2 * sum((c + b^(n - 1) * d - betaj) * b^(n - 1) * log(b) * d)
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
  fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]$par
}

score_miage <- function(id, cpgs, DNAm, partial_cache = NULL) {
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)

  params <- miage_params(id)
  # Panel order is the author's summation order; keep it through the subset so
  # the residual sum accumulates in the same sequence.
  present <- cpgs$score_present
  cached <- if (is.null(partial_cache)) {
    character(0)
  } else {
    intersect(present, colnames(partial_cache))
  }
  betas <- DNAm[, present, drop = FALSE]
  if (length(cached)) {
    betas[, cached] <- partial_cache[, cached]
  }

  b <- params$b[present]
  cc <- params$c[present]
  d <- params$d[present]
  score <- matrix(
    vapply(seq_len(n), function(i) miage_fit(betas[i, ], b, cc, d), numeric(1)),
    nrow = n,
    ncol = 1L,
    dimnames = list(sample_id, id)
  )

  sample_miss <- if (length(cached)) {
    slideimp::mat_miss(DNAm[, cached, drop = FALSE], col = FALSE)
  } else {
    integer(n)
  }
  names(sample_miss) <- sample_id

  coverage <- list(
    clock_id = id,
    policy = clock_impute(id)[["policy"]],
    score_needed = length(cpgs$score_needed),
    score_present = length(present),
    score_used = length(present),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = 0L,
    score_dropped = length(cpgs$score_absent),
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = cpgs$score_absent
  )

  list(score = score, coverage = coverage, sample_miss = sample_miss)
}
