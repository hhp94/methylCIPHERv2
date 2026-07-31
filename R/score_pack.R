# batched scorers for external packs (PCClocks, PCBrainAge, SystemsAge)

# pack cpg panel from the upstream-resolved split (the clock's declared panel)
pack_design <- function(id, cpgs, block) {
  present <- cpgs[["score_present"]]
  absent <- cpgs[["score_absent"]]
  obs <- observed_panel(present, block)
  list(
    present = present,
    absent = absent,
    used = obs[["cols"]],
    X = obs[["values"]],
    ref = clock_impute_ref(id, block[["packs"]])
  )
}

# vendor-mean-filled linear predictors over cols
pack_linpred <- function(design, M, cols) {
  contrib <- design[["X"]] %*% M[design[["used"]], cols, drop = FALSE]
  if (length(design[["absent"]])) {
    off <- as.numeric(
      design[["ref"]][design[["absent"]]] %*%
        M[design[["absent"]], cols, drop = FALSE]
    )
    contrib <- sweep(contrib, 2L, off, "+")
  }
  colnames(contrib) <- cols
  contrib
}

# per-clock covariate contributions (n x k)
pack_cov_contrib <- function(ids, pheno, n) {
  cc <- lapply(ids, clock_covariates_coefs)
  need <- unique(unlist(lapply(cc, names), use.names = FALSE))
  if (!length(need)) {
    return(matrix(0, nrow = n, ncol = length(ids)))
  }
  # presence is a front-door check (check_pheno)
  Cmat <- matrix(0, length(need), length(ids), dimnames = list(need, ids))
  for (j in seq_along(ids)) {
    v <- cc[[j]]
    if (length(v)) {
      Cmat[names(v), ids[j]] <- v
    }
  }
  as.matrix(pheno[, need, drop = FALSE]) %*% Cmat
}

# dispatch a pack group to its batched scorer (cpgs: the group's shared panel)
score_pack_group <- function(ids, cpgs, block) {
  ty <- score_type(ids[[1]])
  switch(
    ty,
    pack_systemsage = score_systemsage_group(ids, cpgs, block),
    pack_linear = score_linear_pack(ids, cpgs, block),
    stop(
      sprintf("No batched scorer for score_type %s.", ty),
      call. = FALSE
    )
  )
}

# coefficient_matrix packs (PCClocks, PCBrainAge)
score_linear_pack <- function(ids, cpgs, block) {
  # every clock here is declared vendor_mean + sum
  pack <- clock_pack(ids[[1]], block[["packs"]])
  M <- pack[["coefficient_matrix"]]
  rownames(M) <- pack[["cpgs"]]
  design <- pack_design(ids[[1]], cpgs, block)
  sample_id <- block[["sample_id"]]

  linpred <- sweep(
    pack_linpred(design, M, ids),
    2L,
    vapply(ids, clock_intercept, numeric(1)),
    "+"
  ) +
    pack_cov_contrib(ids, block[["pheno"]], length(sample_id))

  out <- vector("list", length(ids))
  names(out) <- ids
  for (id in ids) {
    tf <- resolve_output_transform(clock_output_transform(id))
    out[[id]] <- score_matrix(tf(linpred[, id]), sample_id, id)
  }
  out
}
