# dunedin pace-of-aging family (PoAm linear, PACE QN then linear)
score_Dunedin <- function(id, cpgs, block, results) {
  # poAm is vendor-filled linear. rest is PACE quantile normalization.
  if (!cpgs[["normalizes"]]) {
    return(linear_score(cpgs, block))
  }
  require_betanorm(id)

  sample_id <- block[["sample_id"]]
  n <- length(sample_id)
  needed <- cpgs[["norm_needed"]]
  absent <- cpgs[["norm_absent"]]
  target <- clock_norm_target(id)

  # absent background CpG takes the gold target's value (QN needs full panel)
  panel <- matrix(0, nrow = n, ncol = length(needed), dimnames = list(sample_id, needed))
  obs <- observed_panel(cpgs[["norm_present"]], block)
  if (length(obs[["cols"]])) {
    panel[, obs[["cols"]]] <- obs[["values"]]
  }
  if (length(absent)) {
    panel[, absent] <- rep(target[absent], each = n)
  }

  norm <- betanorm::quantile_norm(panel, target = as.numeric(target[needed]))
  dimnames(norm) <- dimnames(panel)

  model <- cpgs[["score_needed"]]
  score_matrix(
    clock_intercept(id) +
      norm[, model, drop = FALSE] %*% clock_coefs(id)[model],
    sample_id,
    id
  )
}
