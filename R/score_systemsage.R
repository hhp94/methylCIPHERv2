# systemsAge: organ sub-clocks plus Age_prediction / SystemsAge composites

# batched scorer for the SystemsAge group
score_systemsage_group <- function(ids, block) {
  packs <- block[["packs"]]
  pack <- clock_pack(ids[[1]], packs)
  design <- pack_design(ids[[1]], pack, block)
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)

  composites <- intersect(ids, c("Age_prediction", "SystemsAge"))
  organs_req <- setdiff(ids, composites)

  record <- function(id, score_vec) {
    score_matrix(score_vec, sample_id, id)
  }

  out <- vector("list", length(ids))
  names(out) <- ids

  if (length(organs_req)) {
    Mo <- pack[["organs"]]
    rownames(Mo) <- pack[["cpgs"]]
    O <- pack_linpred(design, Mo, organs_req)
    for (org in organs_req) {
      out[[org]] <- record(org, O[, org] + clock_intercept(org))
    }
  }

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
        poly_eval(L, systemsage_poly("Age_prediction", "score"))
      )
    }

    if ("SystemsAge" %in% composites) {
      id <- "SystemsAge"
      L <- age_matmul + systemsage_age_intercept(id)
      ap_scaled <- poly_eval(L, systemsage_poly(id, "ap_scaled"))

      # stack step says which column the scaled age front lands in
      map <- systemsage_stack_map(id)
      order <- unname(map)
      ap_col <- unname(map[["ap_scaled"]])
      organs_pca <- setdiff(order, ap_col)
      Ms <- pack[["systems"]]
      rownames(Ms) <- pack[["cpgs"]]
      S <- sweep(
        pack_linpred(design, Ms, organs_pca),
        2L,
        systemsage_raw_intercepts(id)[organs_pca],
        "+"
      )

      sysscores <- matrix(0, n, length(order), dimnames = list(NULL, order))
      sysscores[, ap_col] <- ap_scaled
      sysscores[, organs_pca] <- S[, organs_pca]

      pca <- systemsage_pca(id, packs, order)
      cs <- sweep(sweep(sysscores, 2L, pca[["center"]], "-"), 2L, pca[["scale"]], "/")
      pcs <- cs %*% pca[["rotation"]]
      out[[id]] <- record(
        id,
        as.numeric(systemsage_final_intercept(id) + pcs %*% pca[["model"]])
      )
    }
  }
  out
}

# unique recipe step producing out, or error (>1 already stopped upstream)
systemsage_step <- function(id, out) {
  step <- recipe_step_out(id, out)
  if (is.null(step)) {
    cli::cli_abort(
      c(
        "{.val {id}} has no recipe step with out {.field {out}}.",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  step
}

# intercept of the age-linear front
systemsage_age_intercept <- function(id) {
  as.numeric(systemsage_step(id, "L")[["intercept"]])
}

# quadratic poly coef for the poly step producing out
systemsage_poly <- function(id, out) {
  as.numeric(unlist(systemsage_step(id, out)[["coef"]]))
}

# organ labels come from the stack step (not from stripping raw_)
systemsage_stack_map <- function(id) {
  stack_label_map(systemsage_step(id, "sysscores"), id)
}

# the linear steps producing one stack column each
systemsage_raw_steps <- function(id) {
  operands <- names(systemsage_stack_map(id))
  Filter(
    function(s) {
      identical(s[["op"]], "linear") &&
        !is.null(s[["out"]]) &&
        as.character(s[["out"]]) %in% operands
    },
    clock_entry(id)[["recipe"]]
  )
}

# system linear intercepts, named by stack column label
systemsage_raw_intercepts <- function(id) {
  map <- systemsage_stack_map(id)
  steps <- systemsage_raw_steps(id)
  stats::setNames(
    vapply(steps, function(s) as.numeric(s[["intercept"]]), numeric(1)),
    vapply(
      steps,
      function(s) unname(map[[as.character(s[["out"]])]]),
      character(1)
    )
  )
}

# stack column labels, in column order
systemsage_stack_order <- function(id) {
  unname(systemsage_stack_map(id))
}

# intercept of the final systems_model linear head
systemsage_final_intercept <- function(id) {
  as.numeric(systemsage_step(id, "score")[["intercept"]])
}

# systems_PCA tensors from the pack, ordered to stack rows/PC cols
systemsage_pca <- function(id, packs, order) {
  pack <- clock_pack(id, packs)
  comps <- clock_components(id)
  tensor_by_component <- function(name) {
    comp <- component_named(comps, name, id)
    t <- pack[["tensors"]][[comp[["file"]]]]
    if (is.null(t)) {
      cli::cli_abort(
        c(
          "Pack for {.val {id}} has no tensor {.file {comp[['file']]}}
           (component {.field {name}}).",
          CATALOG_BUG
        ),
        call = NULL
      )
    }
    t
  }

  center <- tensor_by_component("systems_pca_center")
  scale <- tensor_by_component("systems_pca_scale")
  model <- tensor_by_component("systems_model")
  rot_df <- tensor_by_component("systems_pca_rotation")

  # rotation row key is declared -- no first-column fallback
  sys_col <- component_row_key(component_named(comps, "systems_pca_rotation", id))
  if (!sys_col %in% names(rot_df)) {
    cli::cli_abort(
      c(
        "{.val {id}}: rotation tensor has no declared row_key column
         {.field {sys_col}}.",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  pc_cols <- setdiff(names(rot_df), sys_col)
  rot <- as.matrix(rot_df[, pc_cols, drop = FALSE])
  rownames(rot) <- as.character(rot_df[[sys_col]])

  list(
    center = center[order],
    scale = scale[order],
    rotation = rot[order, pc_cols, drop = FALSE],
    model = model[pc_cols]
  )
}
