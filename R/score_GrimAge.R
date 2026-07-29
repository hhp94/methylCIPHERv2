# grimAgeV1/V2: Cox stack of surrogates + Age/Female, then rescale to years
score_GrimAge <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)
  pheno <- block[["pheno"]]

  group_id <- clock_group_id(id)
  cox <- grimage_cox_coef(id)
  comps <- clock_components(id)
  stack_names <- names(cox)

  X <- matrix(
    0,
    nrow = n,
    ncol = length(stack_names),
    dimnames = list(sample_id, stack_names)
  )
  # operand source is the declared namespace, not the name
  roles <- grimage_stack_roles(id, stack_names)
  for (nm in stack_names) {
    if (identical(roles[[nm]], "covariates")) {
      # presence is a front-door check (check_pheno)
      X[, nm] <- as.numeric(pheno[[nm]])
    } else if (identical(roles[[nm]], "internal")) {
      comp <- component_named(comps, nm, id)
      coef <- bundle_tensor(group_id, comp[["file"]])
      intercept <- if (is.null(comp[["intercept"]])) 0 else comp[["intercept"]]
      present <- intersect(names(coef), block[["usable"]])
      # surrogate follows the clock's declared absent-CpG policy
      X[, nm] <- component_linpred(
        id,
        coef,
        present,
        block,
        label = nm,
        intercept = intercept,
        cov_coefs = covariate_coefs_from(comp[["covariates"]])
      )
    } else {
      X[, nm] <- as.numeric(results[[nm]])
    }
  }

  cox_score <- X[, stack_names, drop = FALSE] %*% cox
  score_matrix(
    grimage_rescale(cox_score, grimage_rescale_params(id)),
    sample_id,
    id
  )
}

# cox scale -> years
grimage_rescale <- function(cox_score, params) {
  (cox_score - params[["m_cox"]]) /
    params[["sd_cox"]] *
    params[["sd_age"]] +
    params[["m_age"]]
}

# grimage cox coef vector
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
  p <- step[["params"]]
  miss <- setdiff(need, names(p))
  if (length(miss)) {
    catalog_bug(
      "%s: grimage_rescale transform lacks param(s) %s.",
      id,
      paste(miss, collapse = ", ")
    )
  }
  vapply(p[need], as.numeric, numeric(1))
}

# cox stack operands in coefficient order, tagged by namespace
grimage_stack_roles <- function(id, order) {
  stack_roles(stack_step(id))[order]
}
