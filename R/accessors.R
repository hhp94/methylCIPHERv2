# catalog accessors -- always `[[`, never `$`

# shared suffix for catalog/sync faults
CATALOG_BUG <- c("i" = "Catalog/sync bug -- please report it.")

# one catalog entry, or error
clock_entry <- function(id) {
  if (length(id) != 1L) {
    cli::cli_abort(
      "{.fn clock_entry} takes a single clock id, got {length(id)}.",
      call = NULL
    )
  }
  entry <- mc_catalog[[id]]
  if (is.null(entry)) {
    cli::cli_abort("Unknown clock id: {.val {id}}.", call = NULL)
  }
  entry
}

# a catalog field that must be declared, or a stop naming it
required_field <- function(id, field) {
  value <- clock_entry(id)[[field]]
  if (is.null(value)) {
    cli::cli_abort(
      c("Catalog entry {.val {id}} is missing {.field {field}}.", CATALOG_BUG),
      call = NULL
    )
  }
  value
}

# a catalog field that may be absent, with the declared default
optional_field <- function(id, field, default) {
  value <- clock_entry(id)[[field]]
  if (is.null(value)) default else value
}

# `x` must carry every name in `need`, or a stop listing what is absent
require_fields <- function(x, need, what, id) {
  miss <- setdiff(need, names(x))
  if (length(miss)) {
    cli::cli_abort(
      c(
        "{.val {id}}: {what} lacks {cli::qty(miss)} field{?s} {.field {miss}}.",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  x
}

# named numeric tensor for a weights/ path
bundle_tensor <- function(group_id, path) {
  bundle <- mc_bundles[[group_id]]
  if (is.null(bundle)) {
    cli::cli_abort(
      c(
        "No shipped bundle for group {.val {group_id}}.",
        "i" = "External or unshipped group?"
      ),
      call = NULL
    )
  }
  tensor <- bundle[["tensors"]][[path]]
  if (is.null(tensor)) {
    cli::cli_abort(
      c(
        "Bundle {.val {group_id}} has no tensor {.file {path}}.",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  tensor
}

# the single element of `items` matching `pred`, or a cli stop naming `what`
pick_one <- function(items, pred, what, id) {
  hits <- Filter(pred, items)
  if (length(hits) != 1L) {
    cli::cli_abort(
      c("{.val {id}} has {length(hits)} {what} (expected 1).", CATALOG_BUG),
      call = NULL
    )
  }
  hits[[1]]
}

# the single component with this row_key, resolved to its bundled tensor
component_tensor <- function(id, row_key) {
  entry <- clock_entry(id)
  comp <- pick_one(
    entry[["components"]],
    function(c) identical(c[["row_key"]], row_key),
    sprintf("%s-keyed components", row_key),
    id
  )
  bundle_tensor(entry[["group_id"]], comp[["file"]])
}

# the single component declaration with this name
component_named <- function(comps, name, id) {
  pick_one(
    comps,
    function(c) identical(c[["name"]], name),
    sprintf("components named '%s'", name),
    id
  )
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

# scoring CpGs for one clock (pack panel for external groups)
clock_scoring_cpgs <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  if (isTRUE(entry[["external_group"]])) {
    return(clock_pack(id, packs)[["cpgs"]])
  }
  probe_sets_cpgs(entry, "scoring")
}

# probe_set roles carrying a background panel plus the target it calibrates onto
NORM_ROLES <- c("quantile_normalization_background", "bmiq_gold_standard")

# schemes expressible as a declared panel + a vendored target
NORM_SCHEMES <- c("quantile", "bmiq")

# schemes that are part of the clock's definition and cannot be declined
NORM_CONSTITUTIVE <- "quantile"

# normalization panel for one clock, character(0) unless it normalizes
clock_norm_cpgs <- function(id, normalize = FALSE) {
  if (!isTRUE(normalize)) {
    return(character(0))
  }
  entry <- clock_entry(id)
  cpgs <- unlist(
    lapply(NORM_ROLES, function(role) probe_sets_cpgs(entry, role)),
    use.names = FALSE
  )
  if (is.null(cpgs)) character(0) else unique(cpgs)
}

# declared background probe_set (any NORM_ROLES role), or error
norm_background_probe_set <- function(id) {
  ps <- pick_one(
    clock_entry(id)[["probe_sets"]],
    function(p) p[["role"]] %in% NORM_ROLES,
    "normalization background probe_sets",
    id
  )
  if (!length(ps[["file"]])) {
    cli::cli_abort(
      c(
        "{.val {id}}: normalization background probe_set has no file pointer.",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  ps
}

# vendored normalization target, or NULL when the scheme is not expressible
clock_norm_target <- function(id) {
  if (!(clock_norm_scheme(id) %in% NORM_SCHEMES)) {
    return(NULL)
  }
  ps <- norm_background_probe_set(id)
  bundle_tensor(clock_group_id(id), ps[["file"]])
}

# package-minted classification: "clock", or "sex_routed_alias"
clock_kind <- function(id) {
  as.character(optional_field(id, "kind", "clock"))
}

# {female, male} member ids for a sex-routed alias, or NULL
clock_routing <- function(id) {
  clock_entry(id)[["routing"]]
}

# member clock_id -> sex and alias (callable-pool source)
sex_routed_members <- function() {
  sex <- character(0)
  alias <- character(0)
  aliases <- mc_index[["clock_id"]][
    mc_index[["kind"]] == "sex_routed_alias"
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

# declared array-normalization scheme, lowercased; "none" when absent
clock_norm_scheme <- function(id) {
  scheme <- tolower(as.character(optional_field(id, "normalization", "none")))
  if (!length(scheme)) {
    return("none")
  }
  if (length(scheme) > 1L) {
    cli::cli_abort(
      c(
        "{.val {id}} declares {length(scheme)} normalization schemes
         ({.val {scheme}}); exactly one is supported.",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  scheme
}

# per-cohort parity fixture blocks (one per registry cohort), or NULL
clock_fixtures <- function(id) {
  clock_entry(id)[["fixtures"]]
}

# the fixture block for one cohort, or NULL when the clock has none there
clock_fixture <- function(id, cohort) {
  hits <- Filter(
    function(fx) identical(as.character(fx[["cohort"]]), cohort),
    clock_fixtures(id) %||% list()
  )
  if (!length(hits)) NULL else hits[[1L]]
}

# imputation policy and ref
clock_impute <- function(id) {
  clock_entry(id)[["imputation"]]
}

# vendor-mean ref for absent-CpG fill
clock_impute_ref <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  if (isTRUE(entry[["external_group"]])) {
    pack <- clock_pack(id, packs)
    ref <- pack[["impute"]]
    if (is.null(ref)) {
      cli::cli_abort(
        c(
          "External group {.val {entry[['group_id']]}} pack has no
           impute vector.",
          CATALOG_BUG
        ),
        call = NULL
      )
    }
    ref <- as.numeric(ref)
    names(ref) <- pack[["cpgs"]]
    return(ref)
  }
  imp <- entry[["imputation"]]
  ref <- imp[["ref"]]
  if (is.null(ref) || !is.character(ref) || length(ref) != 1L || !nzchar(ref)) {
    cli::cli_abort(
      c(
        "{.val {id}} has no scalar vendor-mean ref path
         (policy {.val {imp[['policy']]}}).",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  bundle_tensor(entry[["group_id"]], ref)
}

# linear reduction: "mean" if recipe has linear_mean, else "sum"
clock_reduction <- function(id) {
  ops <- vapply(
    clock_entry(id)[["recipe"]],
    function(s) as.character(s[["op"]]),
    character(1)
  )
  if ("linear_mean" %in% ops) "mean" else "sum"
}

# named cpg->coef for single-vector clocks
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
        cli::cli_abort(
          c(
            "External clock {.val {id}} is not a column of group
             {.val {entry[['group_id']]}}'s coefficient_matrix.",
            CATALOG_BUG
          ),
          call = NULL
        )
      }
      coef <- as.numeric(m[, id])
      names(coef) <- pack[["cpgs"]]
      return(coef)
    }
    path <- entry[["coef_path"]]
    if (length(path) != 1L || !nzchar(path)) {
      cli::cli_abort(
        c(
          "{.val {id}} is cpg_coefficient but has no coef tensor path.",
          CATALOG_BUG
        ),
        call = NULL
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
  cli::cli_abort(
    c(
      "{.val {id}} is weights_format {.val {wf}} and does not reduce to a
       single cpg->coef vector.",
      "i" = "Use {.fn clock_group_bundle} + {.fn bundle_tensor}, or the
             clock's family orchestrator."
    ),
    call = NULL
  )
}

# {name: coef} list -> named numeric, numeric(0) if empty
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

# unique recipe step producing `out`, or NULL when the clock has no such step
recipe_step_out <- function(id, out) {
  step <- Filter(
    function(s) identical(s[["out"]], out),
    clock_entry(id)[["recipe"]]
  )
  if (length(step) == 1L) step[[1]] else NULL
}

# covariate weights, numeric(0) when none
clock_covariates_coefs <- function(id) {
  cov <- clock_entry(id)[["covariates"]]
  if (is.null(cov)) {
    step <- recipe_step_out(id, "score")
    cov <- if (is.null(step)) NULL else step[["covariates"]]
  }
  covariate_coefs_from(cov)
}

# covariate names required for pheno checks
clock_covariates_required <- function(id) {
  as.character(optional_field(id, "covariates_required", character(0)))
}

# clock ids this clock consumes as inputs (the recipe's `inputs`, or a
# sex-routed alias's two members)
clock_depends_on <- function(id) {
  as.character(optional_field(id, "clock_inputs", character(0)))
}

# model intercept, 0 when unset
clock_intercept <- function(id) {
  optional_field(id, "intercept", 0)
}

clock_type <- function(id) {
  required_field(id, "computation_type")
}

# catalog output_transform name, default "identity"
clock_output_transform <- function(id) {
  as.character(optional_field(id, "output_transform", "identity"))
}

# weights_format from the catalog
clock_weights_format <- function(id) {
  required_field(id, "weights_format")
}

# true when weights are not in mc_bundles
clock_is_external <- function(id) {
  isTRUE(clock_entry(id)[["external_group"]])
}

# loaded pack for an external clock's group
clock_pack <- function(id, packs) {
  gid <- clock_group_id(id)
  pack <- if (is.null(packs)) NULL else packs[[gid]]
  if (is.null(pack)) {
    cli::cli_abort(
      c(
        "{.val {id}}: pack for external group {.val {gid}} is not loaded.",
        "i" = "Pass {.fn load_mc_assets} output via {.arg packs}."
      ),
      call = NULL
    )
  }
  pack
}

# shipped group bundle: $group_id, $clocks, $tensors
clock_group_bundle <- function(id) {
  group_id <- clock_group_id(id)
  bundle <- mc_bundles[[group_id]]
  if (is.null(bundle)) {
    cli::cli_abort(
      c(
        "No shipped bundle for group {.val {group_id}}.",
        "i" = "External or unshipped group?"
      ),
      call = NULL
    )
  }
  bundle
}

# family/group label for pack dispatch
clock_group_id <- function(id) {
  required_field(id, "group_id")
}

# GrimAge components for score_GrimAge()
clock_components <- function(id) {
  optional_field(id, "components", list())
}

# GrimAge Cox coef vector
grimage_cox_coef <- function(id) {
  component_tensor(id, "component")
}

# grimage_rescale params: m_cox, sd_cox, m_age, sd_age
grimage_rescale_params <- function(id) {
  step <- pick_one(
    clock_entry(id)[["recipe"]],
    function(s) {
      identical(s[["op"]], "transform") &&
        identical(s[["name"]], "grimage_rescale")
    },
    "grimage_rescale transform steps",
    id
  )
  need <- c("m_cox", "sd_cox", "m_age", "sd_age")
  p <- require_fields(
    step[["params"]],
    need,
    "grimage_rescale transform",
    id
  )
  vapply(p[need], as.numeric, numeric(1))
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

# stack column order: inputs, then internal, then covariates
stack_operands <- function(step) {
  c(
    as.character(unlist(step[["inputs"]] %||% character())),
    as.character(unlist(step[["internal"]] %||% character())),
    as.character(unlist(step[["covariates"]] %||% character()))
  )
}

# ordered surrogates: each {name, coef, negate}
physage_surrogates <- function(id) {
  entry <- clock_entry(id)
  recipe <- entry[["recipe"]]

  stack_step <- pick_one(
    recipe,
    function(s) identical(s[["op"]], "stack"),
    "stack ops",
    id
  )
  order <- stack_operands(stack_step)

  zs <- pick_one(
    recipe,
    function(s) {
      identical(s[["op"]], "cohort_zscore") && identical(s[["in"]], "raws")
    },
    "cohort_zscore ops over 'raws'",
    id
  )
  negate_set <- as.character(unlist(zs[["negate"]]))

  lm_ops <- Filter(function(s) identical(s[["op"]], "linear_mean"), recipe)
  by_out <- stats::setNames(
    lm_ops,
    vapply(lm_ops, function(s) s[["out"]], character(1))
  )

  lapply(order, function(raw_name) {
    op <- by_out[[raw_name]]
    if (is.null(op)) {
      cli::cli_abort(
        c(
          "{.val {id}}: stack input {.field {raw_name}} has no matching
           linear_mean op.",
          CATALOG_BUG
        ),
        call = NULL
      )
    }
    comp <- component_named(entry[["components"]], op[["coef"]], id)
    list(
      name = raw_name,
      coef = bundle_tensor(entry[["group_id"]], comp[["file"]]),
      negate = raw_name %in% negate_set
    )
  })
}

# poly coef for DNAmPhysAge_years, or NULL
physage_poly_coef <- function(id) {
  step <- Filter(
    function(s) identical(s[["op"]], "poly"),
    clock_entry(id)[["recipe"]]
  )
  if (!length(step)) {
    return(NULL)
  }
  if (length(step) != 1L) {
    cli::cli_abort(
      c(
        "{.val {id}} has {length(step)} poly ops (expected 0 or 1).",
        CATALOG_BUG
      ),
      call = NULL
    )
  }
  as.numeric(unlist(step[[1]][["coef"]]))
}

# EpiTOC2 per-CpG ground state: named delta (de-novo rate) and beta0 vectors
epitoc2_params <- function(id) {
  tab <- require_fields(
    component_tensor(id, "cpg"),
    c("cpg", "delta", "beta0"),
    "params tensor",
    id
  )
  list(
    delta = stats::setNames(as.numeric(tab[["delta"]]), tab[["cpg"]]),
    beta0 = stats::setNames(as.numeric(tab[["beta0"]]), tab[["cpg"]])
  )
}

# MiAge site-specific params: named b, c, d vectors in panel order
miage_params <- function(id) {
  tab <- component_tensor(id, "cpg")
  lapply(
    list(b = tab[["b"]], c = tab[["c"]], d = tab[["d"]]),
    function(x) stats::setNames(as.numeric(x), tab[["cpg"]])
  )
}

# unique recipe step producing out, or error
systemsage_step <- function(id, out) {
  step <- recipe_step_out(id, out)
  if (is.null(step)) {
    cli::cli_abort(
      c(
        "{.val {id}} has no unique recipe step with out {.field {out}}.",
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

# raw-system linear intercepts, named by organ
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

# stack column order as system labels
systemsage_stack_order <- function(id) {
  inputs <- stack_operands(systemsage_step(id, "sysscores"))
  vapply(
    inputs,
    function(x) {
      if (identical(x, "ap_scaled")) "Age_prediction" else sub("^raw_", "", x)
    },
    character(1),
    USE.NAMES = FALSE
  )
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
