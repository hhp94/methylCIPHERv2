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

# the column a data-frame tensor's rows are keyed by, as the component declares
# it. An absent declaration stops -- there is no first-column fallback.
component_row_key <- function(comp) {
  key <- as.character(comp[["row_key"]])
  if (length(key) != 1L || !nzchar(key)) {
    cli::cli_abort(
      c(
        "Component {.field {comp[['name']]}} declares no {.field row_key}.",
        CATALOG_BUG
      ),
      call = NULL
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

# of those, the ones score_normalized() actually implements. `quantile` is
# missing on purpose: today its only clock is DunedinPACE, which the Dunedin
# group tag claims first, so nothing reaches score_normalized() with it.
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

# 1-based position of the clock's first cohort-reducing recipe step, NA when
# every step is per-sample. Declared by sync from a closed op vocabulary.
clock_cross_sample_at <- function(id) {
  v <- optional_field(id, "cross_sample_at", NA_integer_)
  if (!length(v)) NA_integer_ else as.integer(v)[[1L]]
}

# does this clock's score depend on the other samples in the matrix? A TRUE
# clock is not chunk-safe: scoring a row subset does not reproduce the rows a
# whole-cohort run gives, because the reduction is still inside the branch.
clock_is_cross_sample <- function(id) {
  !is.na(clock_cross_sample_at(id))
}

# split a compute sequence on the sample axis. The one place the kind-1 /
# kind-2 partition is derived, so nothing downstream carries a clock list.
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

# named cpg->coef for one bundled single-vector clock. External clocks never
# reach here: score_type() sends every one of them to a pack_* tag, and the
# batched pack scorers multiply the whole coefficient matrix rather than one
# column at a time.
clock_coefs <- function(id) {
  entry <- clock_entry(id)
  wf <- entry[["weights_format"]]
  if (identical(wf, "cpg_coefficient")) {
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

# unique recipe step producing `out`, or NULL when the clock has no such step.
# Ambiguity is not absence: callers read the returned step's coefficients, so
# folding >1 into NULL would score the clock with the term dropped entirely
# rather than say the catalog is ambiguous.
recipe_step_out <- function(id, out) {
  step <- Filter(
    function(s) identical(s[["out"]], out),
    clock_entry(id)[["recipe"]]
  )
  if (!length(step)) {
    return(NULL)
  }
  if (length(step) > 1L) {
    cli::cli_abort(
      c(
        "{.val {id}} has {length(step)} recipe steps with out
         {.field {out}} (expected at most 1).",
        CATALOG_BUG
      ),
      call = NULL
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

# Does this clock's recipe take per-sample moments over the whole input matrix?
# A `sample_scale` op z-scores each sample over every probe it was handed, so
# the score moves with the width of the caller's matrix, not just with the
# scoring panel. Declared, never a clock list -- the parity tier reads the same
# op to decide which matrix to feed a fixture.
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

# family/group label for pack dispatch
clock_group_id <- function(id) {
  required_field(id, "group_id")
}

# GrimAge components for score_GrimAge()
clock_components <- function(id) {
  optional_field(id, "components", list())
}

# a stack step's three operand namespaces, in column order
STACK_NAMESPACES <- c("inputs", "internal", "covariates")

# operand -> the namespace that declares it. The namespace is what says where a
# stack column comes from: another clock's score, a component this clock
# computes itself, or a pheno column. Only the declaration says which -- the
# operand's name does not, and reading it off the name means a covariate that
# is not Age or Female silently falls through to the wrong source.
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

# operand -> column label, in column order. A label is a matrix index: the
# tensors a later step names (center, scale, rotation, coef) are keyed by it,
# by name. The default label is the operand's own name; an optional `columns`
# list, parallel to the operand concatenation, overrides it elementwise -- so
# a step that declares `columns` is saying its labels are not its operands.
stack_label_map <- function(step, id) {
  operands <- stack_operands(step)
  declared <- step[["columns"]]
  if (is.null(declared)) {
    return(stats::setNames(operands, operands))
  }
  labels <- as.character(unlist(declared))
  if (length(labels) != length(operands)) {
    cli::cli_abort(
      c(
        "{.val {id}}: stack declares {length(labels)} column label{?s} for
         {length(operands)} operand{?s}.",
        CATALOG_BUG
      ),
      call = NULL
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
