MC_ACCEL <- "accel_id"

finalized <- function(x) {
  # n_batches() first. a short-circuit would skip the count check.
  multi <- n_batches(x) > 1L
  pending <- x[["provenance"]][["pending"]]
  if (length(pending) && multi) {
    x <- refinalize_clocks(x)
  }
  x
}

shape_scores <- function(m, id_col, value_col, batch, long, label = NULL) {
  if (!long) {
    ids <- stats::setNames(
      data.frame(rownames(m), stringsAsFactors = FALSE),
      id_col
    )
    out <- cbind(ids, as.data.frame(m, optional = TRUE))
    if (!is.null(label)) {
      names(out) <- c(id_col, paste0(colnames(m), "_", label))
    }
    out[[MC_BATCH]] <- batch
    rownames(out) <- NULL
    return(drop_single_batch(out, batch))
  }
  out <- data.frame(
    id = rep(rownames(m), times = ncol(m)),
    clock_id = rep(colnames(m), each = nrow(m)),
    value = as.vector(m),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  names(out) <- c(id_col, "clock_id", value_col)
  if (is_multi_batch(batch)) {
    out[[MC_BATCH]] <- rep(batch, times = ncol(m))
  }
  if (!is.null(label)) {
    out[[MC_ACCEL]] <- label
    keep <- setdiff(names(out), MC_ACCEL)
    out <- out[append(keep, MC_ACCEL, after = match("clock_id", keep))]
  }
  drop_single_batch(out, batch)
}

#' Data Frame Method For An mc_result Object
#'
#' Converts the [calc_clocks()] output to a data.frame containing just the
#' clocks, in long or wide format.
#'
#' @inheritParams mc-params
#' @param row.names A character vector. Not used by this method. Default is
#'   `NULL`.
#' @param optional A boolean. Not used by this method. Default is
#'   `FALSE`.
#' @param ... Not used.
#'
#' @inheritSection mc-params Clocks that use all the samples
#'
#' @returns A data.frame. In long form, one row for each sample and clock,
#'   with the score and, when `x` holds more than one batch, an
#'   `mc_batch_id` column. In wide form, one row for each sample, with one
#'   column for each clock.
#'
#' @seealso
#' - [as.matrix.mc_result()] for the scores as a numeric matrix.
#' - [rbind.mc_result()] for two runs combined into one object.
#' - [refinalize_clocks()] for a cross-sample score recomputed after a bind.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 20)
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#'
#' head(as.data.frame(res))
#' head(as.data.frame(res, long = FALSE))
#'
#' @export
as.data.frame.mc_result <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  long = TRUE
) {
  check_mc_result(x)
  checkmate::assert_flag(long)
  x <- finalized(x)
  shape_scores(
    x[["scores"]],
    x[["provenance"]][["pheno_id"]],
    "score",
    x[["provenance"]][[MC_BATCH]],
    long
  )
}

type_family <- function(v) {
  if (is.numeric(v) || is.logical(v)) {
    "number"
  } else if (is.character(v) || is.factor(v)) {
    "string"
  } else {
    class(v)[[1L]]
  }
}

values_agree <- function(a, b) {
  if (is.numeric(a) || is.logical(a)) {
    a <- as.numeric(a)
    b <- as.numeric(b)
    same <- abs(a - b) <= 1e-8 * pmax(1, abs(a), abs(b))
  } else {
    same <- as.character(a) == as.character(b)
  }
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & same)
}

column_conflict <- function(a, b, nm) {
  if (type_family(a) != type_family(b)) {
    return(cli::format_inline(
      "{.field {nm}}: {.cls {class(a)[[1L]]}} in {.arg x},
       {.cls {class(b)[[1L]]}} in {.arg data}"
    ))
  }
  n_off <- sum(!values_agree(a, b))
  if (n_off == 0L) {
    return(NULL)
  }
  cli::format_inline("{.field {nm}}: {n_off} sample{?s} disagree")
}

merge_accel_data <- function(pheno, data, pheno_id) {
  if (is.null(data)) {
    return(pheno)
  }
  if (MC_BATCH %in% names(data)) {
    cli::cli_abort(
      c(
        "{.arg data} has a {.field {MC_BATCH}} column.",
        "i" = "That name is reserved for the batch label of the run.",
        "i" = "Rename the column in {.arg data}."
      ),
      call = NULL
    )
  }
  if (!pheno_id %in% names(data)) {
    cli::cli_abort(
      c(
        "{.arg data} has no {.field {pheno_id}} column.",
        "i" = "{.arg data} joins on the sample id column of {.arg x}."
      ),
      call = NULL
    )
  }
  ids <- as.character(data[[pheno_id]])
  dup <- unique(ids[duplicated(ids)])
  if (length(dup)) {
    cli::cli_abort(
      c(
        "{.arg data} has {length(dup)} duplicate id{?s} in
         {.field {pheno_id}}:",
        capped_bullets(dup, val_lines),
        "i" = "One row of {.arg data} per sample id."
      ),
      call = NULL
    )
  }

  # data says nothing about a sample it does not carry
  idx <- id_index(
    as.character(pheno[[pheno_id]]),
    ids,
    "merge_accel_data",
    unmatched = "na"
  )
  seen <- !is.na(idx)
  shared <- setdiff(intersect(names(pheno), names(data)), pheno_id)
  bad <- unlist(lapply(shared, function(cl) {
    column_conflict(pheno[[cl]][seen], data[[cl]][idx[seen]], cl)
  }))
  if (length(bad)) {
    cli::cli_abort(
      c(
        "{.arg data} changes {length(bad)} column{?s} the scores were
         computed from:",
        capped_bullets(bad),
        "i" = "{.arg data} may add a column. It may not change one.",
        "i" = "Call {.fn calc_clocks} again with the pheno you want."
      ),
      call = NULL
    )
  }

  if (any(!seen)) {
    absent <- as.character(pheno[[pheno_id]])[!seen]
    cli::cli_warn(
      c(
        "{.arg data} has no row for {sum(!seen)} sample{?s}:
         {.val {capped_vals(absent)}}.",
        "i" = "An added column is {.code NA} for
               {cli::qty(absent)}{?that sample/those samples}."
      ),
      call = NULL
    )
  }

  for (cl in setdiff(names(data), names(pheno))) {
    pheno[[cl]] <- data[[cl]][idx]
  }
  pheno
}

# the model's rhs, carried as language so terms like I(Age^2) survive
accel_formula <- function(formula, type) {
  if (is.null(formula)) {
    # accel with no formula is the classic age regression.
    if (type == "diff") {
      return(NULL)
    }
    formula <- stats::as.formula("~ Age")
  }
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    cli::cli_abort(
      c(
        "{.arg formula} must be a formula with no left side.",
        "i" = "Put the covariates after {.code ~}. For example,
               {.code ~ Age + Female}."
      ),
      call = NULL
    )
  }
  formula
}

accel_label <- function(formula, type) {
  rhs <- if (is.null(formula)) {
    NULL
  } else {
    attr(stats::terms(formula), "term.labels")
  }
  paste(c(rhs, type), collapse = "_")
}

pattern_residuals <- function(y, ph, formula) {
  mfr <- stats::model.frame(formula, data = ph, na.action = stats::na.fail)
  qrx <- qr(stats::model.matrix(formula, mfr))
  if (nrow(y) - qrx[["rank"]] < 1L) {
    return(NULL)
  }
  qr.resid(qrx, y)
}

residualize <- function(resp, ph, vars, formula) {
  keep <- rep(TRUE, nrow(ph))
  for (v in vars) {
    keep <- keep & !is.na(ph[[v]])
  }

  out <- resp
  out[] <- NA_real_
  okm <- !is.na(resp) & keep
  patt <- apply(okm, 2L, function(v) paste0(which(!v), collapse = ","))
  dead <- character(0)
  for (cols in split(seq_len(ncol(resp)), patt)) {
    ok <- okm[, cols[[1L]]]
    got <- pattern_residuals(
      resp[ok, cols, drop = FALSE],
      ph[ok, , drop = FALSE],
      formula
    )
    if (is.null(got)) {
      dead <- c(dead, colnames(resp)[cols])
      next
    }
    out[ok, cols] <- got
  }
  if (length(dead)) {
    cli::cli_warn(
      c(
        "{length(dead)} clock{?s} had too few complete samples to fit:
         {.val {capped_vals(dead)}}.",
        "i" = "{cli::qty(dead)}{?That column is/Those columns are} all
               {.code NA}."
      ),
      call = NULL
    )
  }
  out
}

say_fill_batch <- function(x, rhs_vars) {
  per_clock <- x[["coverage"]][["per_clock"]]
  n_batch <- n_batches(x)
  if (n_batch < 2L || MC_BATCH %in% rhs_vars) {
    return(invisible(NULL))
  }
  filled <- vapply(
    per_clock,
    function(recs) {
      any(vapply(
        recs,
        function(r) {
          !is.null(r) && as.integer(r[["score_imputed_partial"]]) > 0L
        },
        logical(1L)
      ))
    },
    logical(1L)
  )
  if (!any(filled)) {
    return(invisible(NULL))
  }
  cli::cli_inform(c(
    "!" = "The returned value has {n_batch} batch{?es}, and
           {.fn calc_clocks} filled some absent CpGs with a mean taken inside
           each batch.",
    "i" = "A fill of that kind can shift the scores of one batch against
           another.",
    "i" = "Add {.field {MC_BATCH}} to {.arg formula} so the model accounts for
           the batch."
  ))
  invisible(NULL)
}

#' Age Acceleration Or Difference
#'
#' Computes age acceleration or age difference for every clock in `x`.
#'
#' @inheritParams mc-params
#' @param formula A one-sided formula. The model fit against each clock's
#'   score. Default is `NULL`, which uses `~ Age`.
#' @param type One of "accel" or "diff". The quantity to compute for each
#'   clock. Default is `"accel"`.
#' @param data A data.frame. Extra sample metadata, joined to the `pheno` in `x`
#'   by sample id. Default is `NULL`.
#'
#' @inheritSection mc-params Clocks that use all the samples
#'
#' @details
#' The default `type = "accel"` calculates age acceleration.
#' It regresses each clock in `x` on `Age` and returns the residuals.
#' `type = "diff"` calculates the raw difference between each clock and
#' `Age`, and fits no model unless `formula` is given.
#'
#' `formula` replaces the default model completely. It does not add to it, so
#' `~ Plate` regresses each clock on the plate alone, and not on age.
#'
#' `data` carries the covariates the calculation needs. Pass a data.frame of
#' `ID` and `Plate` to `data`, with `formula = ~ Age + Plate`.
#' It may add a column, and it may repeat a column that scoring already used.
#' `calc_accel()` stops when a repeated column disagrees with the value
#' scoring used.
#'
#' @returns A data.frame. In long form, one row for each sample and clock,
#'   with the fitted value, an `accel_id` column that names the model, and,
#'   when `x` holds more than one batch, `mc_batch_id`. In wide form, one row
#'   for each sample, with one column for each clock.
#'
#' @seealso
#' [score_associations()] for how each clock tracks age against a reference.
#'
#' @examples
#' clocks <- c("Horvath1", "Hannum")
#' sim <- sim_DNAm(clocks, n = 20, Age = TRUE, Female = TRUE)
#' res <- calc_clocks(sim[["DNAm"]], clocks)
#'
#' # accel with no formula regresses each clock's score on ~ Age
#' head(calc_accel(res, data = sim[["pheno"]]))
#'
#' # a formula replaces the default model, so name every term you want
#' pheno <- sim[["pheno"]]
#' pheno[["Plate"]] <- sample(c("P1", "P2"), nrow(pheno), replace = TRUE)
#' head(calc_accel(res, formula = ~ Age + Plate, data = pheno))
#'
#' # diff with no formula is the raw difference from age, with no model fit
#' head(calc_accel(res, type = "diff", data = sim[["pheno"]]))
#'
#' # diff with a formula residualizes the difference
#' head(calc_accel(res, type = "diff", formula = ~ Age, data = sim[["pheno"]]))
#'
#' @export
calc_accel <- function(
  x,
  formula = NULL,
  type = c("accel", "diff"),
  data = NULL,
  long = TRUE
) {
  check_mc_result(x)
  checkmate::assert_data_frame(data, min.rows = 1, null.ok = TRUE)
  checkmate::assert_flag(long)
  type <- match.arg(type)
  formula <- accel_formula(formula, type)
  x <- finalized(x)

  pheno_id <- x[["provenance"]][["pheno_id"]]
  sample_id <- x[["provenance"]][["sample_id"]]
  rhs_vars <- if (is.null(formula)) character(0) else all.vars(formula)
  vars <- unique(c(if (type == "diff") "Age", rhs_vars))
  say_fill_batch(x, rhs_vars)

  pheno <- merge_accel_data(x[["pheno"]], data, pheno_id)
  pheno[[MC_BATCH]] <- x[["provenance"]][[MC_BATCH]][
    id_index(pheno[[pheno_id]], sample_id, "calc_accel batch")
  ]
  need <- setdiff(vars, names(pheno))
  if (length(need)) {
    cli::cli_abort(
      c(
        "The {.field pheno} of {.arg x} has no {cli::qty(need)}column{?s}
         {.val {capped_vals(need)}}.",
        "i" = "Add {cli::qty(need)}{?it/them} to {.arg data}.",
        "i" = "{.arg data} joins on the sample id column of {.arg x}."
      ),
      call = NULL
    )
  }
  check_pheno(pheno, ID = pheno_id, extra_columns = vars, sample_id = sample_id)
  # by id, never by row order
  ph <- pheno[
    id_index(sample_id, pheno[[pheno_id]], "calc_accel pheno"),
    ,
    drop = FALSE
  ]

  resp <- x[["scores"]]
  if (type == "diff") {
    resp <- resp - ph[["Age"]]
  }
  # diff with no formula is the difference itself -- no fit, no drop
  out <- if (is.null(formula)) {
    resp
  } else {
    residualize(resp, ph, rhs_vars, formula)
  }

  shape_scores(
    out,
    pheno_id,
    "accel",
    x[["provenance"]][[MC_BATCH]],
    long,
    accel_label(formula, type)
  )
}
