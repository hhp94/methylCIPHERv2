# Batched scorers for external packs (PCClocks, PCBrainAge, SystemsAge).

# Shared design over a pack CpG panel: subset matrix, present/absent, vendor ref.
pack_design <- function(pack, usable, DNAm, partial_cache) {
  panel <- pack[["cpgs"]]
  hit <- match(panel, usable, 0L) > 0L
  present <- panel[hit]
  absent <- panel[!hit]
  cached_hit <- if (is.null(partial_cache)) {
    logical(length(present))
  } else {
    match(present, colnames(partial_cache), 0L) > 0L
  }
  cached <- present[cached_hit]
  raw <- present[!cached_hit]
  X <- cbind(
    if (length(cached)) partial_cache[, cached, drop = FALSE] else NULL,
    if (length(raw)) DNAm[, raw, drop = FALSE] else NULL
  )
  if (is.null(X)) {
    X <- matrix(0, nrow = nrow(DNAm), ncol = 0L)
  }
  sample_miss <- if (length(cached)) {
    slideimp::mat_miss(DNAm[, cached, drop = FALSE], col = FALSE)
  } else {
    integer(nrow(DNAm))
  }
  names(sample_miss) <- rownames(DNAm)
  list(
    present = present,
    absent = absent,
    cached = cached,
    used = c(cached, raw),
    X = X,
    ref = stats::setNames(as.numeric(pack[["impute"]]), pack[["cpgs"]]),
    sample_miss = sample_miss
  )
}

# Vendor-mean-filled linear predictors over cols, reusing pack_design.
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

# Per-clock covariate contributions as one n x k matrix.
pack_cov_contrib <- function(ids, pheno, n) {
  cc <- lapply(ids, clock_covariate_coefs)
  need <- unique(unlist(lapply(cc, names), use.names = FALSE))
  if (!length(need)) {
    return(matrix(0, nrow = n, ncol = length(ids)))
  }
  if (is.null(pheno) || !all(need %in% names(pheno))) {
    stop(
      "pack scorer: clock(s) need covariate(s) ",
      paste(need, collapse = ", "),
      " but they are absent from `pheno`.",
      call. = FALSE
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

# Coverage record for a vendor-mean linear pack member.
pack_linear_coverage <- function(cpgs, sample_miss) {
  list(
    clock_id = cpgs$clock_id,
    policy = "vendor_mean",
    score_needed = length(cpgs$score_needed),
    score_present = length(cpgs$score_present),
    score_used = length(cpgs$score_present) + length(cpgs$score_absent),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = length(cpgs$score_absent),
    score_dropped = 0L,
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = cpgs$score_absent
  )
}

# Dispatch a pack group to its batched scorer. Members of a group share a tag.
score_pack_group <- function(
  group_id,
  ids,
  cpg_list,
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
      cpg_list,
      usable,
      DNAm,
      partial_cache,
      pheno,
      packs
    ),
    pack_linear = score_linear_pack(
      ids,
      cpg_list,
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

# Batched scorer for coefficient_matrix packs (PCClocks, PCBrainAge).
score_linear_pack <- function(
  ids,
  cpg_list,
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
    score <- matrix(
      as.numeric(tf(linpred[, id])),
      ncol = 1L,
      dimnames = list(sample_id, id)
    )
    out[[id]] <- list(
      score = score,
      coverage = pack_linear_coverage(
        cpg_list$per_clock[[id]],
        design$sample_miss
      ),
      sample_miss = design$sample_miss
    )
  }
  out
}
