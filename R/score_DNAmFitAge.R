# DNAmFitAge_{Sex}: klemera-doubal mix of upstream member scores
score_DNAmFitAge <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)

  # KDM mixing table (component, weight, center, scale) -- sync asserts the header
  kdm <- component_tensor(id, "component")

  # components are scored first (sequence order)
  score_vec <- numeric(n)
  for (i in seq_len(nrow(kdm))) {
    cv <- as.numeric(results[[kdm[["component"]][i]]])
    score_vec <- score_vec +
      kdm[["weight"]][i] * (cv - kdm[["center"]][i]) / kdm[["scale"]][i]
  }

  score_matrix(score_vec, sample_id, id)
}
