# sex-routed alias: pick the member score matching each sample's sex
score_sex_routed <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)
  route <- clock_routing(id)

  score_vec <- rep(NA_real_, n)
  rows <- sex_rows(block[["pheno"]][["Female"]], n)
  for (key in names(rows)) {
    i <- rows[[key]]
    if (!any(i)) {
      next
    }
    score_vec[i] <- as.numeric(results[[route[[key]]]])[i]
  }

  score_matrix(score_vec, sample_id, id)
}

# routed members: scored for coverage, never a score column
drop_routed_members <- function(ids) {
  setdiff(ids, names(sex_routed_members()[["sex"]]))
}
