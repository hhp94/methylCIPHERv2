# a declaration the accessors require is absent or ambiguous
catalog_bug <- function(fmt, ...) {
  stop(
    sprintf(fmt, ...),
    " Catalog/sync bug -- please report it.",
    call. = FALSE
  )
}

# one catalog entry, or error
clock_entry <- function(id) {
  entry <- mc_catalog[[id]]
  if (is.null(entry)) {
    stop(sprintf("Unknown clock id: %s.", id), call. = FALSE)
  }
  entry
}

# one group entry, or error. never read mc_groups directly.
group_entry <- function(gid) {
  grp <- mc_groups[[gid]]
  if (is.null(grp)) {
    stop(sprintf("Unknown clock group: %s.", gid), call. = FALSE)
  }
  grp
}

# a catalog field that must be declared, or a stop naming it
required_field <- function(id, field) {
  value <- clock_entry(id)[[field]]
  if (is.null(value)) {
    catalog_bug("Catalog entry %s is missing %s.", id, field)
  }
  value
}

# a catalog field that may be absent, with the declared default
optional_field <- function(id, field, default) {
  value <- clock_entry(id)[[field]]
  if (is.null(value)) default else value
}

# named numeric tensor for a weights/ path
bundle_tensor <- function(group_id, path) {
  bundle <- mc_bundles[[group_id]]
  if (is.null(bundle)) {
    stop(
      sprintf(
        "No shipped bundle for group %s. External or unshipped group?",
        group_id
      ),
      call. = FALSE
    )
  }
  tensor <- bundle[["tensors"]][[path]]
  if (is.null(tensor)) {
    catalog_bug("Bundle %s has no tensor %s.", group_id, path)
  }
  tensor
}

# the single element of `items` matching `pred`, or a stop naming `what`
pick_one <- function(items, pred, what, id) {
  hits <- Filter(pred, items)
  # multiplicity matters: an unguarded Filter would silently score the first hit
  if (length(hits) != 1L) {
    catalog_bug("%s has %d %s (expected 1).", id, length(hits), what)
  }
  hits[[1]]
}

# recipe steps declaring this op, empty when none. `op` is not a key (scan).
recipe_steps_op <- function(id, op) {
  Filter(
    function(s) identical(as.character(s[["op"]]), op),
    clock_entry(id)[["recipe"]] %||% list()
  )
}

# the single recipe step declaring this op, or a stop naming it
recipe_step_op <- function(id, op) {
  pick_one(
    clock_entry(id)[["recipe"]] %||% list(),
    function(s) identical(as.character(s[["op"]]), op),
    sprintf("%s ops", op),
    id
  )
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

# the component declaration with this name (sync keys `components` by name)
component_named <- function(comps, name, id) {
  comp <- comps[[name]]
  if (is.null(comp)) {
    catalog_bug("%s has no component named '%s'.", id, name)
  }
  comp
}

# the named component a recipe operand points at, resolved to its tensor
component_tensor_named <- function(id, name) {
  entry <- clock_entry(id)
  comp <- component_named(entry[["components"]], as.character(name), id)
  bundle_tensor(entry[["group_id"]], comp[["file"]])
}

# shared-tensor declaration by name (name -> file for recipe operands).
shared_named <- function(entry, name, id) {
  sh <- entry[["shared"]][[name]]
  if (is.null(sh)) {
    catalog_bug("%s has no shared tensor named '%s'.", id, name)
  }
  sh
}

# clock specific accessors ----

# row-key column a component declares (no first-column fallback)
component_row_key <- function(comp) {
  as.character(comp[["row_key"]])
}

# cpGs of one probe_set role (sync keys `probe_sets` by role)
probe_sets_cpgs <- function(entry, role) {
  cpgs <- entry[["probe_sets"]][[role]][["cpgs"]]
  if (is.null(cpgs)) character(0) else unique(as.character(cpgs))
}

# scoring CpGs for one clock (pack panel for external groups)
clock_scoring_cpgs <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  if (isTRUE(entry[["external_group"]])) {
    return(clock_pack(id, packs)[["cpgs"]])
  }
  probe_sets_cpgs(entry, "scoring")
}

# probe_set roles with a background panel plus the target it calibrates onto.
NORM_ROLES <- c("quantile_normalization_background", "bmiq_gold_standard")

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
  pick_one(
    clock_entry(id)[["probe_sets"]],
    function(p) p[["role"]] %in% NORM_ROLES,
    "normalization background probe_sets",
    id
  )
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

# 1-based index of first cohort-reducing recipe step, or NA
clock_cross_sample_at <- function(id) {
  as.integer(optional_field(id, "cross_sample_at", NA_integer_))[[1L]]
}

# true when the score depends on other samples (not chunk-safe)
clock_is_cross_sample <- function(id) {
  !is.na(clock_cross_sample_at(id))
}

# split a compute sequence into per-sample vs cohort-reducing
split_cross_sample <- function(clock_ids) {
  hit <- vapply(clock_ids, clock_is_cross_sample, logical(1L))
  list(per_sample = clock_ids[!hit], cross_sample = clock_ids[hit])
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

# declared array-normalization scheme, lowercased ("none" when absent)
clock_norm_scheme <- function(id) {
  tolower(as.character(optional_field(id, "normalization", "none")))
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
    ref <- as.numeric(pack[["impute"]])
    names(ref) <- pack[["cpgs"]]
    return(ref)
  }
  bundle_tensor(entry[["group_id"]], entry[["imputation"]][["ref"]])
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

# named cpg->coef for one single-vector clock (pack tensor for external clocks)
clock_coefs <- function(id, packs = NULL) {
  entry <- clock_entry(id)
  wf <- entry[["weights_format"]]
  if (identical(wf, "cpg_coefficient")) {
    if (isTRUE(entry[["external_group"]])) {
      return(pack_tensor(id, packs, entry[["coef_path"]]))
    }
    return(bundle_tensor(entry[["group_id"]], entry[["coef_path"]]))
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
    sprintf(
      paste0(
        "%s is weights_format %s and does not reduce to a ",
        "single cpg->coef vector. Use bundle_tensor(), clock_pack(), ",
        "or the clock's family orchestrator."
      ),
      id,
      wf
    ),
    call. = FALSE
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

# the recipe step producing `out`, or NULL (sync keys `recipe` by out)
recipe_step_out <- function(id, out) {
  clock_entry(id)[["recipe"]][[out]]
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

# recipe steps taking a per-sample z-score, empty when the clock has none
sample_scale_steps <- function(id) {
  recipe_steps_op(id, "sample_scale")
}

# sample_scale step ref, or NULL. use clock_moment_key() to tell no-step from no-ref.
sample_scale_ref_name <- function(id) {
  refs <- unique(unlist(lapply(sample_scale_steps(id), function(s) s[["ref"]])))
  if (!length(refs)) {
    return(NULL)
  }
  if (length(refs) > 1L) {
    catalog_bug(
      "%s declares %d sample_scale refs (expected 1).",
      id,
      length(refs)
    )
  }
  as.character(refs)
}

# whole-matrix domain key when sample_scale has no declared `ref`.
FULL_MOMENT_KEY <- "full"

# moment domain key: NULL if no sample_scale, FULL_MOMENT_KEY if no `ref`, else from the declaration.
clock_moment_key <- function(id) {
  if (!length(sample_scale_steps(id))) {
    return(NULL)
  }
  nm <- sample_scale_ref_name(id)
  if (is.null(nm)) {
    return(FULL_MOMENT_KEY)
  }
  paste0(clock_group_id(id), ":", nm)
}

# true when the recipe z-scores over the full input matrix.
clock_needs_full_panel <- function(id) {
  identical(clock_moment_key(id), FULL_MOMENT_KEY)
}

# cpgs a sample_scale `ref` declares, or NULL for the whole-matrix case.
clock_sample_scale_ref <- function(id) {
  nm <- sample_scale_ref_name(id)
  if (is.null(nm)) {
    return(NULL)
  }
  entry <- clock_entry(id)
  sh <- shared_named(entry, nm, id)
  as.character(bundle_tensor(entry[["group_id"]], sh[["file"]]))
}

# moment domain: NULL if no sample_scale, else key and cpgs (NULL cpgs = every DNAm column).
clock_moment_domain <- function(id) {
  key <- clock_moment_key(id)
  if (is.null(key)) {
    return(NULL)
  }
  list(key = key, cpgs = clock_sample_scale_ref(id))
}

# covariate names required for pheno checks
clock_covariates_required <- function(id) {
  as.character(optional_field(id, "covariates_required", character(0)))
}

# clock ids this clock consumes as inputs
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
    stop(
      sprintf(
        paste0(
          "%s: pack for external group %s is not loaded. ",
          "Pass load_mc_assets() output via packs."
        ),
        id,
        gid
      ),
      call. = FALSE
    )
  }
  pack
}

# named numeric tensor a pack carries under its declared weights/ path
pack_tensor <- function(id, packs, path) {
  pack <- clock_pack(id, packs)
  tensor <- pack[["tensors"]][[path]]
  if (is.null(tensor)) {
    catalog_bug("Pack %s has no tensor %s.", clock_group_id(id), path)
  }
  tensor
}

# family/group label for pack dispatch
clock_group_id <- function(id) {
  required_field(id, "group_id")
}

# grimAge components for score_GrimAge()
clock_components <- function(id) {
  optional_field(id, "components", list())
}

# a stack step's three operand namespaces, in column order.
STACK_NAMESPACES <- c("inputs", "internal", "covariates")

# operand -> declaring namespace (inputs / internal / covariates)
stack_roles <- function(step) {
  ns <- stats::setNames(
    lapply(STACK_NAMESPACES, function(k) {
      as.character(unlist(step[[k]] %||% character()))
    }),
    STACK_NAMESPACES
  )
  stats::setNames(
    rep(names(ns), lengths(ns)),
    unlist(ns, use.names = FALSE)
  )
}

# stack column order: inputs, then internal, then covariates
stack_operands <- function(step) {
  names(stack_roles(step))
}

# operand -> column label (default: operand name, or declared `columns`)
stack_label_map <- function(step, id) {
  operands <- stack_operands(step)
  declared <- step[["columns"]]
  if (is.null(declared)) {
    return(stats::setNames(operands, operands))
  }
  labels <- as.character(unlist(declared))
  if (length(labels) != length(operands)) {
    catalog_bug(
      "%s: stack declares %d column label(s) for %d operand(s).",
      id,
      length(labels),
      length(operands)
    )
  }
  stats::setNames(labels, operands)
}

# stack column labels, in column order
stack_labels <- function(step, id) {
  unname(stack_label_map(step, id))
}

# the one stack step of a clock that has one
stack_step <- function(id) {
  recipe_step_op(id, "stack")
}
