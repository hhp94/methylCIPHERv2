# DNAmPhysAge: surrogate means, reverse-code, cohort z-score, row sum.
# Split at the cohort reduction (recipe step 11): physage_raws() is per-sample
# and runs in the scoring loop, finalize_PhysAge() reduces over samples and
# runs once, after every block is in hand.

# the chunk-safe half -- an n x n_surrogate matrix of reverse-coded raws.
# Each row depends only on itself, so a block's rows equal the same rows of a
# whole-cohort run.
physage_raws <- function(id, cpgs, block, results) {
  sample_id <- block[["sample_id"]]
  n <- length(sample_id)
  surrogates <- physage_surrogates(id)

  cols <- lapply(surrogates, function(s) {
    coef <- s[["coef"]]
    present <- intersect(names(coef), cpgs[["score_present"]])
    absent <- setdiff(names(coef), present)
    fill <- absent_fill(
      id,
      coef,
      absent,
      label = paste0(id, " surrogate ", s[["name"]])
    )
    lp <- linear_predictor(
      coef = coef,
      intercept = 0,
      cov_coefs = numeric(0),
      score_present = present,
      block = block,
      id = s[["name"]]
    )
    # the surrogate is a mean, so a dropped CpG leaves the denominator too
    n_terms <- length(present) + length(fill[["filled"]])
    raw <- (as.numeric(lp[["cpg_contrib"]]) + fill[["offset"]]) / n_terms
    if (s[["negate"]]) -raw else raw
  })

  # built by hand rather than vapply: a 1-row block must still be a matrix,
  # and the rownames are what assembly reorders on
  matrix(
    unlist(cols, use.names = FALSE),
    nrow = n,
    dimnames = list(
      sample_id,
      vapply(surrogates, function(s) s[["name"]], character(1))
    )
  )
}

# Both degeneracies that make scale() return NaN, in one place. rowSums()
# propagates a single NaN column across every sample, so either one silently
# turns the whole cohort into NaN rather than failing where it happened.
check_zscoreable <- function(id, raws) {
  n <- nrow(raws)
  if (n < 2L) {
    cli::cli_abort(
      c(
        "{.val {id}} needs at least 2 samples (cohort z-score), got {n}.",
        "i" = "Score it with a larger DNAm matrix."
      ),
      call = NULL
    )
  }

  sds <- matrixStats::colSds(raws)
  flat <- colnames(raws)[!is.finite(sds) | sds == 0]
  if (!length(flat)) {
    return(invisible(TRUE))
  }

  # declared panel size per surrogate -- the count that explains a flat one
  surrogates <- physage_surrogates(id)
  needed <- stats::setNames(
    vapply(surrogates, function(s) length(s[["coef"]]), integer(1L)),
    vapply(surrogates, function(s) s[["name"]], character(1L))
  )
  cli::cli_abort(
    c(
      "{.val {id}} cannot be scored: {length(flat)} surrogate{?s} {?is/are}
       constant across the cohort, so {?its/their} z-score is undefined.",
      bullets(vapply(
        flat,
        function(nm) {
          cli::format_inline("{.field {nm}}: {needed[[nm]]} declared CpG{?s}")
        },
        character(1L)
      )),
      "i" = "A surrogate goes constant when none of its CpGs were observed --
             every sample then gets the same vendor-filled value.",
      "i" = "{.fn clocks_coverage} reports the panel counts for {.val {id}}."
    ),
    call = NULL
  )
}

# the cohort reduction. DNAmPhysAge reduces once (scale -> row_sum ->
# transform); DNAmPhysAge_years reduces twice (a second scale before the
# polynomial), so this follows the branch rather than a step index.
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
