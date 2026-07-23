# Catalog accessors for scorers; never read raw mc_catalog lists.
# Every read here uses `[[` -- `$` partial-matches on lists, so `entry$covariates`
# silently resolves to `covariates_required` on entries that lack the exact field.

# One catalog entry, or error.
clock_entry <- function(id) {
  if (length(id) != 1L) {
    stop(
      "clock_entry() takes a single clock id, got ",
      length(id),
      call. = FALSE
    )
  }
  entry <- mc_catalog[[id]]
  if (is.null(entry)) {
    stop("Unknown clock id: ", id, call. = FALSE)
  }
  entry
}

# Named numeric tensor for a weights/ path.
bundle_tensor <- function(group_id, path) {
  bundle <- mc_bundles[[group_id]]
  if (is.null(bundle)) {
    stop(
      "No shipped bundle for group: ",
      group_id,
      " (external/unshipped group?)",
      call. = FALSE
    )
  }
  tensor <- bundle[["tensors"]][[path]]
  if (is.null(tensor)) {
    stop("Tensor not in bundle ", group_id, ": ", path, call. = FALSE)
  }
  tensor
}

probe_sets_cpgs <- function(entry, role) {
  hits <- Filter(
    function(p) identical(p[["role"]], role),
    entry[["probe_sets"]]
  )
  if (!length(hits)) {
    return(character(0))
  }
  unique(unlist(lapply(hits, function(p) p[["cpgs"]]), use.names = FALSE))
}

# Scoring CpGs for one clock. External groups keep the panel in their pack, so
# the pack is the single source there and the catalog never carries a copy.
clock_scoring_cpgs <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  if (isTRUE(entry[["external_group"]])) {
    return(clock_pack(id, packs)[["cpgs"]])
  }
  probe_sets_cpgs(entry, "scoring")
}

# Normalization/background panel; character(0) when none.
clock_norm_cpgs <- function(id) {
  probe_sets_cpgs(clock_entry(id), "quantile_normalization_background")
}

# The declared quantile-normalization background probe_set, or error.
qn_background_probe_set <- function(id) {
  entry <- clock_entry(id)
  hits <- Filter(
    function(p) identical(p[["role"]], "quantile_normalization_background"),
    entry[["probe_sets"]]
  )
  if (length(hits) != 1L || !length(hits[[1]][["file"]] %||% NULL)) {
    stop(
      "qn_background_probe_set(): '",
      id,
      "' declares normalization 'quantile' but has ",
      length(hits),
      " quantile_normalization_background probe_set(s) with a file pointer ",
      "(expected 1).",
      call. = FALSE
    )
  }
  hits[[1]]
}

# Gold-standard QN target means, or NULL when the clock needs no QN.
dunedin_gold_means <- function(id) {
  if (!identical(clock_norm_scheme(id), "quantile")) {
    return(NULL)
  }
  ps <- qn_background_probe_set(id)
  bundle_tensor(clock_group_id(id), ps[["file"]])
}

# {female, male} member ids for a sex-routed alias, or NULL.
clock_routing <- function(id) {
  clock_entry(id)[["routing"]]
}

# Routing tables flattened: member clock_id -> its sex, and -> its alias. The
# single source for the callable pool, the not-callable error and its
# suggestion, so those three can never disagree.
sex_routed_members <- function() {
  sex <- character(0)
  alias <- character(0)
  aliases <- mc_index[["clock_id"]][
    mc_index[["computation_type"]] == "sex_routed"
  ]
  for (cid in aliases) {
    route <- clock_routing(cid)
    for (sx in names(route)) {
      member <- as.character(route[[sx]])
      sex[[member]] <- sx
      alias[[member]] <- cid
    }
  }
  list(sex = sex, alias = alias)
}

# Paper-assumed array-normalization scheme (annotation only).
clock_norm_scheme <- function(id) {
  scheme <- clock_entry(id)[["normalization"]]
  if (is.null(scheme)) {
    return(character(0))
  }
  as.character(scheme)
}

# Parity fixture stub for tests, or NULL.
clock_fixture <- function(id) {
  clock_entry(id)[["fixture"]]
}

# Imputation policy and ref.
clock_impute <- function(id) {
  clock_entry(id)[["imputation"]]
}

# Vendor-mean ref for absent-CpG fill.
clock_impute_ref <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  if (isTRUE(entry[["external_group"]])) {
    pack <- clock_pack(id, packs)
    ref <- pack[["impute"]]
    if (is.null(ref)) {
      stop(
        "clock_impute_ref(): external group '",
        entry[["group_id"]],
        "' pack has no $impute vector.",
        call. = FALSE
      )
    }
    ref <- as.numeric(ref)
    names(ref) <- pack[["cpgs"]]
    return(ref)
  }
  imp <- entry[["imputation"]]
  ref <- imp[["ref"]]
  if (is.null(ref) || !is.character(ref) || length(ref) != 1L || !nzchar(ref)) {
    stop(
      "clock_impute_ref(): '",
      id,
      "' has no scalar vendor-mean ref path (policy '",
      if (is.null(imp[["policy"]])) NA else imp[["policy"]],
      "').",
      call. = FALSE
    )
  }
  bundle_tensor(entry[["group_id"]], ref)
}

# Linear reduction: "mean" if recipe has linear_mean, else "sum".
clock_reduction <- function(id) {
  ops <- vapply(
    clock_entry(id)[["recipe"]],
    function(s) as.character(s[["op"]]),
    character(1)
  )
  if ("linear_mean" %in% ops) "mean" else "sum"
}

# TRUE when scores depend on which samples are scored together.
clock_batch_dependent <- function(id) {
  isTRUE(clock_entry(id)[["batch_dependent"]])
}

# Named cpg->coef for single-vector clocks.
clock_coefs <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  wf <- entry[["weights_format"]]
  if (identical(wf, "cpg_coefficient")) {
    if (isTRUE(entry[["external_group"]])) {
      pack <- clock_pack(id, packs)
      m <- if (identical(entry[["group_id"]], "SystemsAge")) {
        pack[["organs"]]
      } else {
        pack[["coefficient_matrix"]]
      }
      if (is.null(m) || is.null(colnames(m)) || !id %in% colnames(m)) {
        stop(
          "clock_coefs(): external clock '",
          id,
          "' is not a column of group '",
          entry[["group_id"]],
          "' coefficient_matrix.",
          call. = FALSE
        )
      }
      coef <- as.numeric(m[, id])
      names(coef) <- pack[["cpgs"]]
      return(coef)
    }
    path <- entry[["coef_path"]]
    if (length(path) != 1L || !nzchar(path)) {
      stop(
        "clock_coefs(): ",
        id,
        " is cpg_coefficient but has no coef tensor path.",
        call. = FALSE
      )
    }
    return(bundle_tensor(entry[["group_id"]], path))
  }
  if (identical(wf, "component_matrices")) {
    cpg_comps <- Filter(
      function(c) identical(c[["row_key"]], "cpg"),
      entry[["components"]]
    )
    if (length(cpg_comps) == 1L) {
      return(bundle_tensor(entry[["group_id"]], cpg_comps[[1]][["file"]]))
    }
  }
  stop(
    "clock_coefs(): ",
    id,
    " is weights_format='",
    wf,
    "' and does not reduce to a single cpg->coef vector; ",
    "use clock_group_bundle() + bundle_tensor() (or its family orchestrator).",
    call. = FALSE
  )
}

# {name: coef} list -> named numeric; numeric(0) if empty.
covariate_coefs_from <- function(cov) {
  empty <- stats::setNames(numeric(0), character(0))
  if (is.null(cov) || !length(cov)) {
    return(empty)
  }
  nms <- names(cov)
  if (is.null(nms) || !all(nzchar(nms))) {
    return(empty)
  }
  stats::setNames(vapply(cov, as.numeric, numeric(1)), nms)
}

# Unique recipe step producing `out`; NULL when the clock has no such step.
recipe_step_out <- function(id, out) {
  step <- Filter(
    function(s) identical(s[["out"]], out),
    clock_entry(id)[["recipe"]]
  )
  if (length(step) == 1L) step[[1]] else NULL
}

# Covariate weights; numeric(0) when none. Single-tensor clocks carry them at
# top level, recipe-borne ones on the step that produces the score.
clock_covariate_coefs <- function(id) {
  cov <- clock_entry(id)[["covariates"]]
  if (is.null(cov)) {
    step <- recipe_step_out(id, "score")
    cov <- if (is.null(step)) NULL else step[["covariates"]]
  }
  covariate_coefs_from(cov)
}

# Covariate names required for pheno checks.
clock_covariates_required <- function(id) {
  covs <- clock_entry(id)[["covariates_required"]]
  if (is.null(covs) || length(covs) == 0) character(0) else as.character(covs)
}

# Clock ids this clock consumes as inputs.
clock_depends_on <- function(id) {
  deps <- clock_entry(id)[["depends_on_clocks"]]
  if (is.null(deps) || length(deps) == 0) character(0) else as.character(deps)
}

# Model intercept; 0 when unset.
clock_intercept <- function(id) {
  intercept <- clock_entry(id)[["intercept"]]
  if (is.null(intercept)) 0 else intercept
}

clock_type <- function(id) {
  computation_type <- clock_entry(id)[["computation_type"]]
  if (is.null(computation_type)) {
    stop("Unexpected Error: all `computation_type` should be not `NULL`")
  } else {
    computation_type
  }
}

# Catalog output_transform name; default "identity".
clock_output_transform <- function(id) {
  ot <- clock_entry(id)[["output_transform"]]
  if (is.null(ot)) "identity" else as.character(ot)
}

# weights_format from the catalog.
clock_weights_format <- function(id) {
  weights_format <- clock_entry(id)[["weights_format"]]
  if (is.null(weights_format)) {
    stop("Unexpected Error: all `weights_format` should be not `NULL`")
  } else {
    weights_format
  }
}

# TRUE when weights are not in mc_bundles.
clock_is_external <- function(id) {
  isTRUE(clock_entry(id)[["external_group"]])
}

# Loaded pack for an external clock's group.
clock_pack <- function(id, packs) {
  gid <- clock_group_id(id)
  pack <- if (is.null(packs)) NULL else packs[[gid]]
  if (is.null(pack)) {
    stop(
      "clock '",
      id,
      "': pack for external group '",
      gid,
      "' is not loaded; pass the registry from load_mc_assets() via `packs`.",
      call. = FALSE
    )
  }
  pack
}

# Shipped group bundle: $group_id, $clocks, $tensors.
clock_group_bundle <- function(id) {
  group_id <- clock_entry(id)[["group_id"]]
  bundle <- mc_bundles[[group_id]]
  if (is.null(bundle)) {
    stop("No shipped bundle for group: ", group_id, call. = FALSE)
  }
  bundle
}

# Family/group label for pack dispatch.
clock_group_id <- function(id) {
  gid <- clock_entry(id)[["group_id"]]
  if (is.null(gid)) {
    stop("Unexpected Error: clock '", id, "' has no group_id", call. = FALSE)
  }
  gid
}

# GrimAge components for score_grimage().
clock_components <- function(id) {
  comps <- clock_entry(id)[["components"]]
  if (is.null(comps)) list() else comps
}

# GrimAge Cox coef vector.
grimage_cox_coef <- function(id) {
  entry <- clock_entry(id)
  model <- Filter(
    function(c) identical(c[["row_key"]], "component"),
    entry[["components"]]
  )
  if (length(model) != 1L) {
    stop(
      "grimage_cox_coef(): ",
      id,
      " has ",
      length(model),
      " component-keyed model matrices (expected exactly 1).",
      call. = FALSE
    )
  }
  bundle_tensor(entry[["group_id"]], model[[1]][["file"]])
}

# grimage_rescale params: m_cox, sd_cox, m_age, sd_age.
grimage_rescale_params <- function(id) {
  recipe <- clock_entry(id)[["recipe"]]
  step <- Filter(
    function(s) {
      identical(s[["op"]], "transform") &&
        identical(s[["name"]], "grimage_rescale")
    },
    recipe
  )
  if (length(step) != 1L) {
    stop(
      "grimage_rescale_params(): ",
      id,
      " has no unique grimage_rescale transform step.",
      call. = FALSE
    )
  }
  p <- step[[1]][["params"]]
  need <- c("m_cox", "sd_cox", "m_age", "sd_age")
  miss <- setdiff(need, names(p))
  if (length(miss)) {
    stop(
      "grimage_rescale_params(): ",
      id,
      " transform is missing param(s) ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  vapply(p[need], as.numeric, numeric(1))
}

# Klemera-Doubal mixing table for a DNAmFitAge_{Sex} composite: one row per
# input clock, with weight/center/scale.
fitage_kdm_params <- function(id) {
  entry <- clock_entry(id)
  comp <- Filter(
    function(c) identical(c[["row_key"]], "component"),
    entry[["components"]]
  )
  if (length(comp) != 1L) {
    stop(
      "fitage_kdm_params(): ",
      id,
      " has ",
      length(comp),
      " component-keyed KDM table(s) (expected 1).",
      call. = FALSE
    )
  }
  tab <- bundle_tensor(entry[["group_id"]], comp[[1]][["file"]])
  miss <- setdiff(c("component", "weight", "center", "scale"), names(tab))
  if (length(miss)) {
    stop(
      "fitage_kdm_params(): ",
      id,
      " KDM table lacks column(s) ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  tab
}

# Ordered surrogates: each {name, coef, negate}.
physage_surrogates <- function(id) {
  entry <- clock_entry(id)
  recipe <- entry[["recipe"]]

  stack_step <- Filter(function(s) identical(s[["op"]], "stack"), recipe)
  if (length(stack_step) != 1L) {
    stop(
      "physage_surrogates(): ",
      id,
      " has ",
      length(stack_step),
      " stack op(s) (expected 1).",
      call. = FALSE
    )
  }
  order <- as.character(unlist(stack_step[[1]][["inputs"]]))

  zs <- Filter(
    function(s) {
      identical(s[["op"]], "cohort_zscore") && identical(s[["in"]], "raws")
    },
    recipe
  )
  if (length(zs) != 1L) {
    stop(
      "physage_surrogates(): ",
      id,
      " has ",
      length(zs),
      " cohort_zscore op(s) over 'raws' (expected 1).",
      call. = FALSE
    )
  }
  negate_set <- as.character(unlist(zs[[1]][["negate"]]))

  lm_ops <- Filter(function(s) identical(s[["op"]], "linear_mean"), recipe)
  by_out <- stats::setNames(
    lm_ops,
    vapply(lm_ops, function(s) s[["out"]], character(1))
  )

  lapply(order, function(raw_name) {
    op <- by_out[[raw_name]]
    if (is.null(op)) {
      stop(
        "physage_surrogates(): ",
        id,
        " stack input '",
        raw_name,
        "' has no matching linear_mean op.",
        call. = FALSE
      )
    }
    comp <- Filter(
      function(c) identical(c[["name"]], op[["coef"]]),
      entry[["components"]]
    )
    if (length(comp) != 1L) {
      stop(
        "physage_surrogates(): ",
        id,
        " component '",
        op[["coef"]],
        "' resolves to ",
        length(comp),
        " tensor(s) (expected 1).",
        call. = FALSE
      )
    }
    list(
      name = raw_name,
      coef = bundle_tensor(entry[["group_id"]], comp[[1]][["file"]]),
      negate = raw_name %in% negate_set
    )
  })
}

# Poly coef for DNAmPhysAge_years, or NULL.
physage_poly_coef <- function(id) {
  step <- Filter(
    function(s) identical(s[["op"]], "poly"),
    clock_entry(id)[["recipe"]]
  )
  if (!length(step)) {
    return(NULL)
  }
  if (length(step) != 1L) {
    stop(
      "physage_poly_coef(): ",
      id,
      " has ",
      length(step),
      " poly op(s) (expected 0 or 1).",
      call. = FALSE
    )
  }
  as.numeric(unlist(step[[1]][["coef"]]))
}

# EpiTOC2 per-CpG ground state: named delta (de-novo rate) and beta0 vectors.
epitoc2_params <- function(id) {
  entry <- clock_entry(id)
  comp <- Filter(
    function(c) identical(c[["row_key"]], "cpg"),
    entry[["components"]]
  )
  if (length(comp) != 1L) {
    stop(
      "epitoc2_params(): ",
      id,
      " has ",
      length(comp),
      " cpg-keyed component(s) (expected 1).",
      call. = FALSE
    )
  }
  tab <- bundle_tensor(entry[["group_id"]], comp[[1]][["file"]])
  miss <- setdiff(c("cpg", "delta", "beta0"), names(tab))
  if (length(miss)) {
    stop(
      "epitoc2_params(): ",
      id,
      " params tensor lacks column(s) ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  list(
    delta = stats::setNames(as.numeric(tab[["delta"]]), tab[["cpg"]]),
    beta0 = stats::setNames(as.numeric(tab[["beta0"]]), tab[["cpg"]])
  )
}

# MiAge site-specific params: named b, c, d vectors in panel order. The tensor
# and its component entry are both written by sync's custom-group registry.
miage_params <- function(id) {
  entry <- clock_entry(id)
  comp <- Filter(
    function(c) identical(c[["row_key"]], "cpg"),
    entry[["components"]]
  )
  tab <- bundle_tensor(entry[["group_id"]], comp[[1]][["file"]])
  lapply(
    list(b = tab[["b"]], c = tab[["c"]], d = tab[["d"]]),
    function(x) stats::setNames(as.numeric(x), tab[["cpg"]])
  )
}

# Unique recipe step producing out, or error.
systemsage_step <- function(id, out) {
  step <- recipe_step_out(id, out)
  if (is.null(step)) {
    stop(
      "systemsage_step(): ",
      id,
      " has no unique recipe step with out='",
      out,
      "'.",
      call. = FALSE
    )
  }
  step
}

# Intercept of the age-linear front.
systemsage_age_intercept <- function(id) {
  as.numeric(systemsage_step(id, "L")[["intercept"]])
}

# Quadratic poly coef for the poly step producing out.
systemsage_poly <- function(id, out) {
  as.numeric(unlist(systemsage_step(id, out)[["coef"]]))
}

# Raw-system linear intercepts, named by organ.
systemsage_raw_intercepts <- function(id) {
  steps <- Filter(
    function(s) {
      identical(s[["op"]], "linear") &&
        !is.null(s[["out"]]) &&
        startsWith(s[["out"]], "raw_")
    },
    clock_entry(id)[["recipe"]]
  )
  ints <- vapply(steps, function(s) as.numeric(s[["intercept"]]), numeric(1))
  names(ints) <- sub(
    "^raw_",
    "",
    vapply(steps, function(s) s[["out"]], character(1))
  )
  ints
}

# Stack column order as system labels.
systemsage_stack_order <- function(id) {
  inputs <- as.character(unlist(systemsage_step(id, "sysscores")[["inputs"]]))
  vapply(
    inputs,
    function(x) {
      if (identical(x, "ap_scaled")) "Age_prediction" else sub("^raw_", "", x)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# Intercept of the final systems_model linear head.
systemsage_final_intercept <- function(id) {
  as.numeric(systemsage_step(id, "score")[["intercept"]])
}

# systems_PCA tensors from the pack, ordered to stack rows/PC cols.
systemsage_pca <- function(id, packs, order) {
  pack <- clock_pack(id, packs)
  comps <- clock_components(id)
  tensor_by_component <- function(name) {
    comp <- Filter(function(c) identical(c[["name"]], name), comps)
    if (length(comp) != 1L) {
      stop(
        "systemsage_pca(): ",
        id,
        " lacks component '",
        name,
        "'.",
        call. = FALSE
      )
    }
    t <- pack[["tensors"]][[comp[[1]][["file"]]]]
    if (is.null(t)) {
      stop(
        "systemsage_pca(): pack for '",
        id,
        "' has no tensor ",
        comp[[1]][["file"]],
        " (component '",
        name,
        "').",
        call. = FALSE
      )
    }
    t
  }

  center <- tensor_by_component("systems_pca_center")
  scale <- tensor_by_component("systems_pca_scale")
  model <- tensor_by_component("systems_model")
  rot_df <- tensor_by_component("systems_pca_rotation")

  sys_col <- if ("system" %in% names(rot_df)) "system" else names(rot_df)[1]
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
