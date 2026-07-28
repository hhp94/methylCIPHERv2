# DNAmFitAge_{Sex}: Klemera-Doubal mix of upstream member scores
score_DNAmFitAge <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)

  kdm <- fitage_kdm_params(id)

  # components are scored first -- resolve_clocks_sequence() guarantees the order
  score_vec <- numeric(n)
  for (i in seq_len(nrow(kdm))) {
    cv <- as.numeric(results[[kdm[["component"]][i]]])
    score_vec <- score_vec +
      kdm[["weight"]][i] * (cv - kdm[["center"]][i]) / kdm[["scale"]][i]
  }

  score_matrix(score_vec, sample_id, id)
}

# klemera-Doubal mixing table (component, weight, center, scale)
fitage_kdm_params <- function(id) {
  require_fields(
    component_tensor(id, "component"),
    c("component", "weight", "center", "scale"),
    "KDM table",
    id
  )
}
