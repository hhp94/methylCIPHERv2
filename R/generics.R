# result-record methods for class "mc_result"

# compact view of a scored result: dims, clock ids, a corner of the scores
#' @export
print.mc_result <- function(x, n = 6L, k = 4L, ...) {
  s <- x[["scores"]]
  cat(sprintf("<mc_result> %d sample(s) x %d clock(s)\n", nrow(s), ncol(s)))
  ids <- colnames(s)
  shown <- utils::head(ids, 8L)
  cat(sprintf(
    "clocks: %s%s\n\n",
    paste(shown, collapse = ", "),
    if (length(ids) > length(shown)) sprintf(", ... (%d total)", length(ids)) else ""
  ))
  nr <- min(n, nrow(s))
  nc <- min(k, ncol(s))
  print(round(s[seq_len(nr), seq_len(nc), drop = FALSE], 3L))
  if (nr < nrow(s) || nc < ncol(s)) {
    cat(sprintf("... %d more row(s), %d more col(s)\n", nrow(s) - nr, ncol(s) - nc))
  }
  cat("\nExtract scores with $scores; per-clock QC with clocks_coverage().\n")
  invisible(x)
}

# scores as a plain n x k double matrix
#' @export
as.matrix.mc_result <- function(x, ...) {
  x[["scores"]]
}

# scores as a data.frame: the id column plus one column per clock. No pheno --
# that stays on $pheno / augment(), so it cannot leak into an analysis by accident.
#' @export
as.data.frame.mc_result <- function(x, ...) {
  s <- x[["scores"]]
  id <- x[["provenance"]][["pheno_id"]]
  out <- data.frame(
    .id = rownames(s),
    s,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  names(out)[1L] <- id
  out
}

# analysis-ready table: scores + id + the aligned covariates the record carries
# (Age/Female/...), plus an optional user data.frame joined by the id column.
# With adjust=, additionally append per-clock residual columns (age acceleration
# when adjust = "Age") -- see augment_residualize().
#
#   augment(res)                          # scores + the run's covariates
#   augment(res, data = pheno)            # + user covariates joined by id
#   augment(res, adjust = c("Age","Female"))  # + <clock>_resid columns
#' @export
augment <- function(x, data = NULL, adjust = NULL) {
  check_mc_result(x)
  id <- x[["provenance"]][["pheno_id"]]
  out <- as.data.frame(x)

  pheno <- x[["pheno"]]
  if (!is.null(pheno) && ncol(pheno) > 1L) {
    out <- merge(out, pheno, by = id, all.x = TRUE, sort = FALSE)
  }
  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (!id %in% names(data)) {
      cli::cli_abort(
        c(
          "{.arg data} must have the id column {.field {id}} to join on.",
          "i" = "It is joined to the scores by sample id."
        ),
        call = NULL
      )
    }
    # duplicate ids would fan out the join and multiply score rows -- refuse
    if (anyDuplicated(data[[id]])) {
      cli::cli_abort(
        c(
          "{.arg data} has duplicate ids in {.field {id}}.",
          "i" = "Each sample id must appear once, or the scores get multiplied."
        ),
        call = NULL
      )
    }
    # a column already on the table would be silently suffixed .x/.y -- flag it
    clash <- setdiff(intersect(names(out), names(data)), id)
    if (length(clash)) {
      cli::cli_warn(
        "{.arg data} {cli::qty(clash)}column{?s} {.field {clash}} {cli::qty(clash)}{?is/are}
         already in the table; both kept, suffixed {.val .x}/{.val .y}.",
        call = NULL
      )
    }
    out <- merge(out, data, by = id, all.x = TRUE, sort = FALSE)
  }

  if (!is.null(adjust)) {
    out <- augment_residualize(out, x[["provenance"]][["clocks"]], adjust)
  }
  out
}

# append one <clock>_resid column per clock = residual of lm(clock ~ adjust)
# fit within this table. adjust on "Age" gives age acceleration. The residual is
# cohort-dependent (fit over the supplied samples), so it shifts if the sample
# set changes.
augment_residualize <- function(out, clocks, adjust) {
  checkmate::assert_character(adjust, min.len = 1L, any.missing = FALSE)
  missing <- setdiff(adjust, names(out))
  if (length(missing)) {
    cli::cli_abort(
      c(
        "{.arg adjust} {cli::qty(missing)}column{?s} {.field {missing}}
         {cli::qty(missing)}{?is/are} not in the table.",
        "i" = "The record only carries covariates a clock required -- pass the rest
               via {.arg data}."
      ),
      call = NULL
    )
  }
  clocks <- intersect(clocks, names(out))
  form_rhs <- paste(sprintf("`%s`", adjust), collapse = " + ")
  for (cl in clocks) {
    form <- stats::as.formula(sprintf("`%s` ~ %s", cl, form_rhs))
    out[[paste0(cl, "_resid")]] <- tryCatch(
      as.numeric(stats::residuals(
        stats::lm(form, data = out, na.action = stats::na.exclude)
      )),
      error = function(e) rep(NA_real_, nrow(out))
    )
  }
  out
}

# rbind is refused -- re-run calc_clocks() on the combined DNAm instead
#' @export
rbind.mc_result <- function(...) {
  cli::cli_abort(
    c(
      "{.cls mc_result} records cannot be {.fn rbind}-ed.",
      "i" = "Re-run {.fn calc_clocks} on the combined DNAm -- batch-dependent
             clocks must see all samples at once."
    ),
    call = NULL
  )
}
