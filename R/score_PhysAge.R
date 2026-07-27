# DNAmPhysAge: surrogate means, reverse-code, cohort z-score, row sum.
# Split at the cohort reduction (recipe step 11): physage_raws() is per-sample
# and runs in the scoring loop, finalize_PhysAge() reduces over samples and
# runs once, after every block is in hand.

# the chunk-safe half -- an n x n_surrogate matrix of reverse-coded raws.
# Each row depends only on itself, so a block's rows equal the same rows of a
# whole-cohort run.
physage_raws <- function(id, cpgs, DNAm, partial_cache = NULL) {
  n <- nrow(DNAm)
  surrogates <- physage_surrogates(id)
  ref <- clock_impute_ref(id)

  cols <- lapply(surrogates, function(s) {
    coef <- s$coef
    present <- intersect(names(coef), cpgs$score_present)
    absent <- setdiff(names(coef), present)
    absent_offset <- vendor_offset(
      coef,
      absent,
      ref,
      paste0(id, " surrogate ", s$name)
    )
    lp <- linear_predictor(
      coef = coef,
      intercept = 0,
      cov_coefs = numeric(0),
      score_present = present,
      DNAm = DNAm,
      partial_cache = partial_cache,
      id = s$name
    )
    raw <- (as.numeric(lp$cpg_contrib) + absent_offset) / length(coef)
    if (s$negate) -raw else raw
  })

  # built by hand rather than vapply: a 1-row block must still be a matrix,
  # and the rownames are what assembly reorders on
  matrix(
    unlist(cols, use.names = FALSE),
    nrow = n,
    dimnames = list(
      rownames(DNAm),
      vapply(surrogates, function(s) s$name, character(1))
    )
  )
}

# the cohort reduction. DNAmPhysAge reduces once (scale -> row_sum ->
# transform); DNAmPhysAge_years reduces twice (a second scale before the
# polynomial), so this follows the branch rather than a step index.
finalize_PhysAge <- function(id, raws) {
  n <- nrow(raws)
  if (n < 2L) {
    cli::cli_abort(
      c(
        "{.val {id}} needs at least 2 samples (cohort z-score), got {n}.",
        "i" = "Score it with a larger DNAm matrix."
      ),
      call = NULL
    )
  }

  phys <- rowSums(scale(raws))

  poly <- physage_poly_coef(id)
  score_vec <- if (is.null(poly)) {
    phys
  } else {
    poly_eval(as.numeric(scale(phys)), poly)
  }

  score_matrix(score_vec, rownames(raws), id)
}

# y = sum_k coef[k+1] * x^k (lowest degree first)
poly_eval <- function(x, coef) {
  powers <- vapply(seq_along(coef) - 1L, function(k) x^k, numeric(length(x)))
  as.numeric(powers %*% coef)
}
