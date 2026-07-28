# catalog accessors -- always `[[`, never `$`

# shared suffix for catalog/sync faults
CATALOG_BUG <- "Catalog/sync bug -- please report it."

# one catalog entry, or error
clock_entry <- function(id) {
  if (length(id) != 1L) {
    stop(
      sprintf("clock_entry() takes a single clock id, got %d.", length(id)),
      call. = FALSE
    )
  }
  entry <- mc_catalog[[id]]
  if (is.null(entry)) {
    stop(sprintf("Unknown clock id: %s.", id), call. = FALSE)
  }
  entry
}

# a catalog field that must be declared, or a stop naming it
required_field <- function(id, field) {
  value <- clock_entry(id)[[field]]
  if (is.null(value)) {
    stop(
      sprintf(
        "Catalog entry %s is missing %s. %s",
        id,
        field,
        CATALOG_BUG
      ),
      call. = FALSE
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
    stop(
      sprintf(
        "%s: %s lacks field(s) %s. %s",
        id,
        what,
        paste(miss, collapse = ", "),
        CATALOG_BUG
      ),
      call. = FALSE
    )
  }
  x
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
    stop(
      sprintf(
        "Bundle %s has no tensor %s. %s",
        group_id,
        path,
        CATALOG_BUG
      ),
      call. = FALSE
    )
  }
  tensor
}

# the single element of `items` matching `pred`, or a stop naming `what`
pick_one <- function(items, pred, what, id) {
  hits <- Filter(pred, items)
  if (length(hits) != 1L) {
    stop(
      sprintf(
        "%s has %d %s (expected 1). %s",
        id,
        length(hits),
        what,
        CATALOG_BUG
      ),
      call. = FALSE
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

# row-key column a component declares (no first-column fallback)
component_row_key <- function(comp) {
  key <- as.character(comp[["row_key"]])
  if (length(key) != 1L || !nzchar(key)) {
    stop(
      sprintf(
        "Component %s declares no row_key. %s",
        comp[["name"]],
        CATALOG_BUG
      ),
      call. = FALSE
    )
  }
  key
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

# schemes score_normalized() implements (quantile routes via Dunedin)
NORM_SCHEMES_ROUTED <- "bmiq"

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
    stop(
      sprintf(
        "%s: normalization background probe_set has no file pointer. %s",
        id,
        CATALOG_BUG
      ),
      call. = FALSE
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

# 1-based index of first cohort-reducing recipe step, or NA
clock_cross_sample_at <- function(id) {
  v <- optional_field(id, "cross_sample_at", NA_integer_)
  if (!length(v)) NA_integer_ else as.integer(v)[[1L]]
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
  scheme <- tolower(as.character(optional_field(id, "normalization", "none")))
  if (!length(scheme)) {
    return("none")
  }
  if (length(scheme) > 1L) {
    stop(
      sprintf(
        paste0(
          "%s declares %d normalization schemes (%s); ",
          "exactly one is supported. %s"
        ),
        id,
        length(scheme),
        paste(scheme, collapse = ", "),
        CATALOG_BUG
      ),
      call. = FALSE
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
      stop(
        sprintf(
          "External group %s pack has no impute vector. %s",
          entry[["group_id"]],
          CATALOG_BUG
        ),
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
      sprintf(
        "%s has no scalar vendor-mean ref path (policy %s). %s",
        id,
        imp[["policy"]],
        CATALOG_BUG
      ),
      call. = FALSE
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

# named cpg->coef for one bundled single-vector clock
clock_coefs <- function(id) {
  entry <- clock_entry(id)
  wf <- entry[["weights_format"]]
  if (identical(wf, "cpg_coefficient")) {
    path <- entry[["coef_path"]]
    if (length(path) != 1L || !nzchar(path)) {
      stop(
        sprintf(
          "%s is cpg_coefficient but has no coef tensor path. %s",
          id,
          CATALOG_BUG
        ),
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
    sprintf(
      paste0(
        "%s is weights_format %s and does not reduce to a ",
        "single cpg->coef vector. Use clock_group_bundle() + ",
        "bundle_tensor(), or the clock's family orchestrator."
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

# unique recipe step producing `out`, or NULL (stops if ambiguous)
recipe_step_out <- function(id, out) {
  step <- Filter(
    function(s) identical(s[["out"]], out),
    clock_entry(id)[["recipe"]]
  )
  if (!length(step)) {
    return(NULL)
  }
  if (length(step) > 1L) {
    stop(
      sprintf(
        "%s has %d recipe steps with out %s (expected at most 1). %s",
        id,
        length(step),
        out,
        CATALOG_BUG
      ),
      call. = FALSE
    )
  }
  step[[1]]
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

# true when the recipe z-scores over the full input matrix (sample_scale)
clock_needs_full_panel <- function(id) {
  any(vapply(
    clock_entry(id)[["recipe"]],
    function(s) identical(as.character(s[["op"]]), "sample_scale"),
    logical(1)
  ))
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

# family/group label for pack dispatch
clock_group_id <- function(id) {
  required_field(id, "group_id")
}

# grimAge components for score_GrimAge()
clock_components <- function(id) {
  optional_field(id, "components", list())
}

# a stack step's three operand namespaces, in column order
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
    stop(
      sprintf(
        "%s: stack declares %d column label(s) for %d operand(s). %s",
        id,
        length(labels),
        length(operands),
        CATALOG_BUG
      ),
      call. = FALSE
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
  pick_one(
    clock_entry(id)[["recipe"]],
    function(s) identical(s[["op"]], "stack"),
    "stack ops",
    id
  )
}
