# shared linear scorers

# present CpGs covered by the cohort-mean cache, as names and usable positions
cached_cols <- function(present, present_idx, block) {
  mask <- block[["cached_mask"]]
  # null when nothing was partly missing -- no panel has a cached column
  if (is.null(mask)) {
    return(list(cols = character(0), idx = integer(0)))
  }
  keep <- mask[present_idx]
  list(cols = present[keep], idx = present_idx[keep])
}

# column positions in the block matrix for usable positions the caller resolved
block_cols <- function(pos, block) {
  # 0 silently drops an element and a high position yields NA -- say so loudly
  bad <- is.na(pos) | pos < 1L | pos > length(block[["usable_idx"]])
  if (any(bad)) {
    stop(
      sprintf(
        paste0(
          "block_cols: %d position(s) outside the block's usable set: %s. ",
          "This is a package bug -- please report it."
        ),
        sum(bad),
        paste(utils::head(pos[bad], 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  block[["usable_idx"]][pos]
}

# observed betas for present: one subset, then cohort-mean columns on top
observed_panel <- function(present, present_idx, block) {
  cached <- cached_cols(present, present_idx, block)
  values <- block[["DNAm"]][, block_cols(present_idx, block), drop = FALSE]
  if (length(cached[["cols"]])) {
    cache <- block[["partial_cache"]]
    values[, cached[["cols"]]] <- cache[, cached[["cols"]], drop = FALSE]
  }
  list(cols = present, values = values)
}

# present from the clock's declared panel, never the block's usable set
component_present <- function(coef, cpgs, label) {
  extra <- setdiff(names(coef), cpgs[["score_needed"]])
  if (length(extra)) {
    catalog_bug(
      "%s: %d coefficient CpG(s) outside the declared scoring panel: %s.",
      label,
      length(extra),
      paste(utils::head(extra, 5L), collapse = ", ")
    )
  }
  # hash the component, not the panel. cols follows the panel's order, which
  # is what keeps it aligned with idx.
  keep <- cpgs[["score_present"]] %in% names(coef)
  list(
    cols = cpgs[["score_present"]][keep],
    idx = cpgs[["score_present_idx"]][keep]
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

# horvath age back-transform, adult age 20
anti_trafo <- function(x) {
  ifelse(x < 0, 21 * exp(x) - 1, 21 * x + 20)
}

# retroelement pan-mammalian back-transform: trained on log(age + 2)
log_offset_anti_trafo <- function(x) {
  exp(x) - 2
}

resolve_output_transform <- function(name) {
  switch(
    name,
    identity = function(x) x,
    anti.trafo = anti_trafo,
    log_offset_anti_trafo = log_offset_anti_trafo,
    stop(sprintf("Unknown output_transform %s.", name), call. = FALSE)
  )
}

# linpred = intercept + sum(coef * beta) + covariates
linear_predictor <- function(
  coef,
  intercept,
  cov_coefs,
  score_present,
  score_idx,
  block,
  observed = NULL
) {
  # observed lets a pre-transform branch supply already-normalized betas
  obs <- if (is.null(observed)) {
    observed_panel(score_present, score_idx, block)
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
  # present is the component's resolved panel: names plus usable positions
  cols <- present[["cols"]]
  fill <- absent_fill(id, coef, setdiff(names(coef), cols), label = label)
  lp <- linear_predictor(
    coef = coef,
    intercept = intercept,
    cov_coefs = cov_coefs,
    score_present = cols,
    score_idx = present[["idx"]],
    block = block
  )
  if (identical(reduction, "mean")) {
    as.numeric(
      mean_cpg_contrib(lp, fill, length(cols)) +
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
  # declared absent set, not setdiff(names(coef), present)
  fill <- absent_fill(id, coef, cpgs[["score_absent"]])

  lp <- linear_predictor(
    coef = coef,
    intercept = intercept,
    cov_coefs = clock_covariates_coefs(id),
    score_present = cpgs[["score_present"]],
    score_idx = cpgs[["score_present_idx"]],
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
