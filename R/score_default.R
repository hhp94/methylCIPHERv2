# shared linear scorers

# present CpGs covered by the cohort-mean cache
cached_cols <- function(present, partial_cache) {
  if (is.null(partial_cache)) {
    character(0)
  } else {
    intersect(present, colnames(partial_cache))
  }
}

# column positions in the block matrix for CpGs the block declares usable
block_cols <- function(cols, block) {
  idx <- block[["usable_idx"]][match(cols, block[["usable"]])]
  # a name outside `usable` would index an NA column, not error -- say so loudly
  if (anyNA(idx)) {
    bad <- cols[is.na(idx)]
    stop(
      sprintf(
        paste0(
          "block_cols: %d CpG(s) outside the block's usable set: %s. ",
          "This is a package bug -- please report it."
        ),
        length(bad),
        paste(utils::head(bad, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  idx
}

# observed betas for `present`: cohort-mean-filled columns first, then raw
observed_panel <- function(present, block) {
  cache <- block[["partial_cache"]]
  cached <- cached_cols(present, cache)
  raw <- setdiff(present, cached)
  list(
    cols = c(cached, raw),
    values = cbind(
      cache[, cached, drop = FALSE],
      block[["DNAm"]][, block_cols(raw, block), drop = FALSE]
    )
  )
}

# vendor-mean fill for fully absent CpGs
vendor_offset <- function(coef, absent, ref, id) {
  miss_ref <- setdiff(absent, names(ref))
  if (length(miss_ref)) {
    stop(
      sprintf(
        "%s: no vendor mean for %d absent CpG(s): %s.",
        id,
        length(miss_ref),
        paste(utils::head(miss_ref, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sum(coef[absent] * ref[absent])
}

# absent-CpG contribution under the declared policy (vendor_mean or drop)
absent_fill <- function(id, coef, absent, ref = NULL, label = id) {
  no_fill <- list(offset = 0, filled = character(0))
  if (!length(absent)) {
    return(no_fill)
  }
  if (!identical(clock_impute(id)[["policy"]], "vendor_mean")) {
    return(no_fill)
  }
  if (is.null(ref)) {
    ref <- clock_impute_ref(id)
  }
  list(offset = vendor_offset(coef, absent, ref, label), filled = absent)
}

# n x 1 score matrix every branch returns
score_matrix <- function(values, sample_id, id) {
  matrix(
    as.numeric(values),
    nrow = length(sample_id),
    ncol = 1L,
    dimnames = list(sample_id, id)
  )
}

# horvath age back-transform (adult age 20)
ADULT_AGE <- 20

anti_trafo <- function(x) {
  ifelse(x < 0, (1 + ADULT_AGE) * exp(x) - 1, (1 + ADULT_AGE) * x + ADULT_AGE)
}

resolve_output_transform <- function(name) {
  switch(
    name,
    identity = function(x) x,
    anti.trafo = anti_trafo,
    stop(sprintf("Unknown output_transform %s.", name), call. = FALSE)
  )
}

# linpred = intercept + sum(coef * beta) + covariates
linear_predictor <- function(
  coef,
  intercept,
  cov_coefs,
  score_present,
  block,
  observed = NULL
) {
  # observed lets a pre-transform branch supply already-normalized betas
  obs <- if (is.null(observed)) {
    observed_panel(score_present, block)
  } else {
    observed
  }
  cpg_contrib <- obs[["values"]] %*% coef[obs[["cols"]]]

  cov_contrib <- 0
  if (length(cov_coefs)) {
    # presence is a front-door check (check_pheno)
    need <- names(cov_coefs)
    cov_mat <- as.matrix(block[["pheno"]][, need, drop = FALSE])
    cov_contrib <- cov_mat %*% cov_coefs[need]
  }

  linpred <- cpg_contrib + cov_contrib + intercept
  list(
    linpred = linpred,
    cpg_contrib = cpg_contrib,
    cov_contrib = cov_contrib
  )
}

# mean-reduced CpG term: filled offset included, over the terms actually used
mean_cpg_contrib <- function(lp, fill, n_present) {
  (lp[["cpg_contrib"]] + fill[["offset"]]) /
    (n_present + length(fill[["filled"]]))
}

# one component's linear predictor: absent fill, then sum or mean the CpG terms
component_linpred <- function(
  id,
  coef,
  present,
  block,
  label = id,
  intercept = 0,
  cov_coefs = numeric(0),
  reduction = c("sum", "mean")
) {
  reduction <- match.arg(reduction)
  fill <- absent_fill(id, coef, setdiff(names(coef), present), label = label)
  lp <- linear_predictor(
    coef = coef,
    intercept = intercept,
    cov_coefs = cov_coefs,
    score_present = present,
    block = block
  )
  if (identical(reduction, "mean")) {
    as.numeric(
      mean_cpg_contrib(lp, fill, length(present)) +
        lp[["cov_contrib"]] +
        intercept
    )
  } else {
    as.numeric(lp[["linpred"]] + fill[["offset"]])
  }
}

# linear engine for one cpg_coefficient clock
linear_score <- function(cpgs, block, observed = NULL) {
  id <- cpgs[["clock_id"]]
  reduction <- clock_reduction(id)
  coef <- clock_coefs(id)
  intercept <- clock_intercept(id)
  # the declared absent set, not setdiff(names(coef), present) -- see
  # component_linpred(), which is for callers holding only a coef vector
  fill <- absent_fill(id, coef, cpgs[["score_absent"]])

  lp <- linear_predictor(
    coef = coef,
    intercept = intercept,
    cov_coefs = clock_covariates_coefs(id),
    score_present = cpgs[["score_present"]],
    block = block,
    observed = observed
  )

  if (identical(reduction, "mean")) {
    linpred <- mean_cpg_contrib(lp, fill, length(cpgs[["score_present"]])) +
      lp[["cov_contrib"]] +
      intercept
  } else {
    linpred <- lp[["linpred"]] + fill[["offset"]]
  }

  transform <- resolve_output_transform(clock_output_transform(id))
  score_matrix(transform(linpred), block[["sample_id"]], id)
}
