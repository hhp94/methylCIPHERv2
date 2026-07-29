# DNAmPhysAge: physage_raws (per-sample) then finalize_PhysAge (cohort reduce)

# per-sample half: n x n_surrogate reverse-coded raws
physage_raws <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)
  surrogates <- physage_surrogates(id)

  cols <- lapply(surrogates, function(s) {
    coef <- s[["coef"]]
    present <- intersect(names(coef), cpgs[["score_present"]])
    # mean over present CpGs only
    raw <- component_linpred(
      id,
      coef,
      present,
      block,
      label = paste0(id, " surrogate ", s[["name"]]),
      reduction = "mean"
    )
    if (s[["negate"]]) -raw else raw
  })

  # hand-built matrix so a 1-row block keeps dim and rownames
  matrix(
    unlist(cols, use.names = FALSE),
    nrow = n,
    dimnames = list(
      sample_id,
      vapply(surrogates, function(s) s[["name"]], character(1))
    )
  )
}

# stop on scale() degeneracies (flat column or n < 2)
check_zscoreable <- function(id, raws) {
  n <- nrow(raws)
  if (n < 2L) {
    stop(
      sprintf(
        paste0(
          "%s needs at least 2 samples (cohort z-score), got %d. ",
          "Score it with a larger DNAm matrix."
        ),
        id,
        n
      ),
      call. = FALSE
    )
  }

  sds <- matrixStats::colSds(raws)
  flat <- colnames(raws)[!is.finite(sds) | sds == 0]
  if (!length(flat)) {
    return(invisible(TRUE))
  }

  # declared panel size per surrogate
  surrogates <- physage_surrogates(id)
  needed <- stats::setNames(
    vapply(surrogates, function(s) length(s[["coef"]]), integer(1L)),
    vapply(surrogates, function(s) s[["name"]], character(1L))
  )
  detail <- paste(
    vapply(
      flat,
      function(nm) sprintf("%s: %d declared CpG(s)", nm, needed[[nm]]),
      character(1L)
    ),
    collapse = "; "
  )
  stop(
    sprintf(
      paste0(
        "%s cannot be scored: %d surrogate(s) constant across the cohort, ",
        "so z-score is undefined (%s). A surrogate goes constant when none ",
        "of its CpGs were observed. clocks_coverage() reports panel counts."
      ),
      id,
      length(flat),
      detail
    ),
    call. = FALSE
  )
}

# cohort reduction (years branch reduces twice before the poly)
finalize_PhysAge <- function(id, raws) {
  check_zscoreable(id, raws)

  phys <- rowSums(scale(raws))

  poly <- physage_poly_coef(id)
  score_vec <- if (is.null(poly)) {
    phys
  } else {
    poly_eval(as.numeric(scale(phys)), poly)
  }

  score_matrix(score_vec, rownames(raws), id)
}

# ordered surrogates: each {name, coef, negate}
physage_surrogates <- function(id) {
  entry <- clock_entry(id)
  recipe <- entry[["recipe"]]

  order <- stack_operands(stack_step(id))

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
    # a stack input with no linear_mean op leaves component_named() 0 hits
    op <- by_out[[raw_name]]
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
  as.numeric(unlist(step[[1]][["coef"]]))
}
