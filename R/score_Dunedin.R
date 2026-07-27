# Dunedin pace-of-aging family (PoAm linear, PACE QN then linear)
score_Dunedin <- function(id, cpgs, DNAm, partial_cache = NULL) {
  sample_id <- rownames(DNAm)
  n <- nrow(DNAm)

  coef <- clock_coefs(id)
  intercept <- clock_intercept(id)
  model_needed <- cpgs$score_needed
  model_present <- cpgs$score_present
  model_absent <- cpgs$score_absent

  # PACE uses the gold QN panel, PoAm the model CpGs
  qn <- isTRUE(cpgs$normalizes)
  if (qn) {
    fill_ref <- clock_norm_target(id)
    panel_needed <- cpgs$norm_needed
    panel_present <- cpgs$norm_present
    panel_absent <- cpgs$norm_absent
  } else {
    fill_ref <- clock_impute_ref(id)
    panel_needed <- model_needed
    panel_present <- model_present
    panel_absent <- model_absent
  }

  score <- score_matrix(NA_real_, sample_id, id)

  obs <- observed_panel(panel_present, DNAm, partial_cache)

  panel <- matrix(
    0,
    nrow = n,
    ncol = length(panel_needed),
    dimnames = list(sample_id, panel_needed)
  )
  if (length(obs$cols)) {
    panel[, obs$cols] <- obs$values
  }
  if (length(panel_absent)) {
    panel[, panel_absent] <- rep(fill_ref[panel_absent], each = n)
  }

  scored <- if (qn) {
    require_betanorm(id)
    norm <- betanorm::quantile_norm(
      panel,
      target = as.numeric(fill_ref[panel_needed])
    )
    dimnames(norm) <- dimnames(panel)
    norm
  } else {
    panel
  }
  score[, 1] <- as.numeric(
    intercept + scored[, model_needed, drop = FALSE] %*% coef[model_needed]
  )
  score
}
