# SystemsAge: organ sub-clocks linear; Age_prediction and SystemsAge composites.

# c0 + c1*L + c2*L^2.
sa_poly <- function(L, coef) {
  out <- rep(0, length(L))
  for (k in seq_along(coef)) {
    out <- out + coef[[k]] * L^(k - 1L)
  }
  out
}

# Coverage record for a SystemsAge composite.
systemsage_composite_coverage <- function(cpgs, sample_miss) {
  list(
    clock_id = cpgs$clock_id,
    policy = "vendor_mean",
    score_needed = length(cpgs$score_needed),
    score_present = length(cpgs$score_present),
    score_used = length(cpgs$score_present),
    score_imputed_partial = sum(sample_miss),
    score_imputed_full = length(cpgs$score_absent),
    score_dropped = 0L,
    norm_needed = length(cpgs$norm_needed),
    norm_present = length(cpgs$norm_present),
    missing_cpgs = cpgs$score_absent
  )
}

# Batched scorer for the SystemsAge group.
score_systemsage_group <- function(
  ids,
  cpg_list,
  usable,
  DNAm,
  partial_cache,
  pheno,
  packs
) {
  pack <- clock_pack(ids[[1]], packs)
  design <- pack_design(pack, usable, DNAm, partial_cache)
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)

  composites <- intersect(ids, c("Age_prediction", "SystemsAge"))
  organs_req <- setdiff(ids, composites)

  record <- function(id, score_vec, coverage) {
    list(
      score = matrix(
        as.numeric(score_vec),
        ncol = 1L,
        dimnames = list(sample_id, id)
      ),
      coverage = coverage,
      sample_miss = design$sample_miss
    )
  }

  out <- vector("list", length(ids))
  names(out) <- ids

  # Organ sub-clocks: plain linear over the pack's `organs` matrix.
  if (length(organs_req)) {
    Mo <- pack[["organs"]]
    rownames(Mo) <- pack[["cpgs"]]
    O <- pack_linpred(design, Mo, organs_req)
    for (org in organs_req) {
      out[[org]] <- record(
        org,
        O[, org] + clock_intercept(org),
        pack_linear_coverage(cpg_list$per_clock[[org]], design$sample_miss)
      )
    }
  }

  # Composites share the age-linear front L.
  if (length(composites)) {
    Ma <- matrix(
      as.numeric(pack[["age"]]),
      ncol = 1L,
      dimnames = list(pack[["cpgs"]], "age")
    )
    age_matmul <- as.numeric(pack_linpred(design, Ma, "age"))

    if ("Age_prediction" %in% composites) {
      L <- age_matmul + systemsage_age_intercept("Age_prediction")
      out[["Age_prediction"]] <- record(
        "Age_prediction",
        sa_poly(L, systemsage_poly("Age_prediction", "score")),
        systemsage_composite_coverage(
          cpg_list$per_clock[["Age_prediction"]],
          design$sample_miss
        )
      )
    }

    if ("SystemsAge" %in% composites) {
      id <- "SystemsAge"
      L <- age_matmul + systemsage_age_intercept(id)
      ap_scaled <- sa_poly(L, systemsage_poly(id, "ap_scaled"))

      order <- systemsage_stack_order(id)
      organs_pca <- setdiff(order, "Age_prediction")
      Ms <- pack[["systems"]]
      rownames(Ms) <- pack[["cpgs"]]
      S <- sweep(
        pack_linpred(design, Ms, organs_pca),
        2L,
        systemsage_raw_intercepts(id)[organs_pca],
        "+"
      )

      sysscores <- matrix(0, n, length(order), dimnames = list(NULL, order))
      sysscores[, "Age_prediction"] <- ap_scaled
      sysscores[, organs_pca] <- S[, organs_pca]

      pca <- systemsage_pca(id, packs, order)
      cs <- sweep(sweep(sysscores, 2L, pca$center, "-"), 2L, pca$scale, "/")
      pcs <- cs %*% pca$rotation
      out[[id]] <- record(
        id,
        as.numeric(systemsage_final_intercept(id) + pcs %*% pca$model),
        systemsage_composite_coverage(
          cpg_list$per_clock[[id]],
          design$sample_miss
        )
      )
    }
  }
  out
}
