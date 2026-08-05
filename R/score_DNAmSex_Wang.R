# dnamSex_Wang: z-score over the autosomal ref, centre sex panel, project onto PC1.

score_DNAmSex_Wang <- function(id, cpgs, block, results) {
  step_center <- recipe_step_op(id, "center_scale")
  # a declared scale would be silently ignored by the centring below
  if (!is.null(step_center[["scale"]])) {
    catalog_bug(
      "%s: center_scale declares a scale this branch does not apply.",
      id
    )
  }
  center <- component_tensor_named(id, step_center[["center"]])
  rotation <- component_tensor_named(
    id,
    recipe_step_op(id, "project")[["rotation"]]
  )
  # both PCA tensors are keyed on the panel -- an absent name would index NA
  if (!setequal(names(center), names(rotation))) {
    catalog_bug("%s: center and rotation cover different CpG sets.", id)
  }

  # the declared panel, never the block's usable set
  panel <- component_present(rotation, cpgs, id)
  present <- panel[["cols"]]
  obs <- observed_panel(present, panel[["idx"]], block)

  # per-sample moments over the declared ref, banked by mc_cohort()
  mom <- block_domain_moments(block, id)

  # matmul plus two scalar reductions over present (same as score_Zhang2019).
  r <- rotation[present]
  lp <- obs[["values"]] %*% r
  # mean/sd are per-sample, so they recycle down the single column
  score <- (lp - mom[["mean"]] * sum(r)) /
    mom[["sd"]] -
    sum(center[present] * r)

  # n < 2 on the ref leaves moments NA. sample is unscorable.
  failed <- block[["sample_id"]][is.na(mom[["sd"]])]
  if (length(failed)) {
    note_scoring_failure(block, id, failed)
    warning(
      sprintf(
        paste0(
          "%s: %d sample(s) have fewer than 2 observed CpGs on the ",
          "z-score reference. Scored NA. ",
          "Also recorded in $provenance$scoring_failures."
        ),
        id,
        length(failed)
      ),
      call. = FALSE
    )
  }

  score_matrix(score, block[["sample_id"]], id)
}
