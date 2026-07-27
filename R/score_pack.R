# batched scorers for external packs (PCClocks, PCBrainAge, SystemsAge)

# pack CpG panel: subset, present/absent, vendor ref
pack_design <- function(pack, usable, DNAm, partial_cache) {
  panel <- pack[["cpgs"]]
  hit <- match(panel, usable, 0L) > 0L
  present <- panel[hit]
  absent <- panel[!hit]
  obs <- observed_panel(present, DNAm, partial_cache)
  list(
    present = present,
    absent = absent,
    used = obs$cols,
    X = obs$values,
    ref = stats::setNames(as.numeric(pack[["impute"]]), pack[["cpgs"]])
  )
}

# vendor-mean-filled linear predictors over cols
pack_linpred <- function(design, M, cols) {
  contrib <- design$X %*% M[design$used, cols, drop = FALSE]
  if (length(design$absent)) {
    off <- as.numeric(
      design$ref[design$absent] %*% M[design$absent, cols, drop = FALSE]
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
  if (is.null(pheno) || !all(need %in% names(pheno))) {
    cli::cli_abort(
      c(
        "These pack clocks need {cli::qty(need)} pheno column{?s}
         {.field {need}}.",
        "i" = "Add {cli::qty(need)}{?it/them} to {.arg pheno}."
      ),
      call = NULL
    )
  }
  Cmat <- matrix(0, length(need), length(ids), dimnames = list(need, ids))
  for (j in seq_along(ids)) {
    v <- cc[[j]]
    if (length(v)) {
      Cmat[names(v), ids[j]] <- v
    }
  }
  as.matrix(pheno[, need, drop = FALSE]) %*% Cmat
}

# dispatch a pack group to its batched scorer
score_pack_group <- function(
  group_id,
  ids,
  usable,
  DNAm,
  partial_cache,
  pheno,
  packs
) {
  switch(
    score_type(ids[[1]]),
    pack_systemsage = score_systemsage_group(
      ids,
      usable,
      DNAm,
      partial_cache,
      pheno,
      packs
    ),
    pack_linear = score_linear_pack(
      ids,
      usable,
      DNAm,
      partial_cache,
      pheno,
      packs
    ),
    stop(
      "score_pack_group(): group '",
      group_id,
      "' has no batched scorer.",
      call. = FALSE
    )
  )
}

# coefficient_matrix packs (PCClocks, PCBrainAge)
score_linear_pack <- function(
  ids,
  usable,
  DNAm,
  partial_cache,
  pheno,
  packs
) {
  pack <- clock_pack(ids[[1]], packs)
  for (id in ids) {
    if (!identical(clock_impute(id)[["policy"]], "vendor_mean")) {
      stop(
        "score_linear_pack(): '",
        id,
        "' policy != vendor_mean.",
        call. = FALSE
      )
    }
    if (!identical(clock_reduction(id), "sum")) {
      stop("score_linear_pack(): '", id, "' reduction != sum.", call. = FALSE)
    }
  }

  M <- pack[["coefficient_matrix"]]
  rownames(M) <- pack[["cpgs"]]
  design <- pack_design(pack, usable, DNAm, partial_cache)

  linpred <- sweep(
    pack_linpred(design, M, ids),
    2L,
    vapply(ids, clock_intercept, numeric(1)),
    "+"
  ) +
    pack_cov_contrib(ids, pheno, nrow(DNAm))

  sample_id <- rownames(DNAm)
  out <- vector("list", length(ids))
  names(out) <- ids
  for (id in ids) {
    tf <- resolve_output_transform(clock_output_transform(id))
    out[[id]] <- score_matrix(tf(linpred[, id]), sample_id, id)
  }
  out
}
