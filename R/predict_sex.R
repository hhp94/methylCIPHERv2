# karyotype call from the two DNAmSex_Wang sex PCs. rules come from the catalog.

# the one group that answers this question
SEX_GROUP <- "DNAmSex_Wang"

# companion columns to the declared karyotype output_column
RECORDED_SEX <- "recorded_sex"
SEX_MISMATCH <- "sex_mismatch"

# non-aneuploid calls the rule table emits, and binary Female labels. aneuploid calls are shown, not flagged.
BINARY_CALLS <- c(female = "Female", male = "Male")

# declared karyotype_call block for SEX_GROUP.
karyotype_spec <- function() {
  kc <- group_entry(SEX_GROUP)[["routing"]][["karyotype_call"]]
  if (is.null(kc)) {
    catalog_bug("Group %s declares no karyotype_call.", SEX_GROUP)
  }
  kc
}

# operand keys a rule set names, in first-appearance order
karyotype_keys <- function(kc) {
  keys <- unique(unlist(lapply(kc[["rules"]], names), use.names = FALSE))
  setdiff(keys, as.character(kc[["output_column"]]))
}

# operand key -> clock id. checked against declared order.
karyotype_inputs <- function(kc) {
  ids <- as.character(unlist(kc[["inputs"]]))
  keys <- karyotype_keys(kc)
  if (length(keys) != length(ids)) {
    catalog_bug(
      "%s: karyotype_call declares %d input(s) for %d rule operand(s).",
      SEX_GROUP,
      length(ids),
      length(keys)
    )
  }
  bad <- !endsWith(tolower(ids), tolower(keys))
  if (any(bad)) {
    catalog_bug(
      "%s: karyotype operand '%s' does not name input %s.",
      SEX_GROUP,
      keys[bad][[1L]],
      ids[bad][[1L]]
    )
  }
  stats::setNames(ids, keys)
}

# declared condition ("<0", ">0") -> predicate. threshold is read, not assumed.
karyotype_predicate <- function(cond) {
  cond <- as.character(cond)
  op <- sub("^([<>]=?).*$", "\\1", cond)
  rhs <- suppressWarnings(as.numeric(sub("^[<>]=?", "", cond)))
  if (!(op %in% c("<", "<=", ">", ">=")) || is.na(rhs)) {
    catalog_bug("%s: cannot read karyotype condition '%s'.", SEX_GROUP, cond)
  }
  function(x) match.fun(op)(x, rhs)
}

# operand scores plus declared rules -> one call per sample.
apply_karyotype <- function(scores, kc) {
  out_col <- as.character(kc[["output_column"]])
  n <- length(scores[[1L]])
  call <- rep(as.character(kc[["default"]]), n)

  for (rule in kc[["rules"]]) {
    hit <- rep(TRUE, n)
    for (key in setdiff(names(rule), out_col)) {
      if (is.null(scores[[key]])) {
        catalog_bug(
          "%s: karyotype rule names unknown operand '%s'.",
          SEX_GROUP,
          key
        )
      }
      hit <- hit & karyotype_predicate(rule[[key]])(scores[[key]])
    }
    # an NA score matches nothing.
    hit[is.na(hit)] <- FALSE
    call[hit] <- as.character(rule[[out_col]])
  }

  # unscorable samples get no call (not a default label).
  call[Reduce(`|`, lapply(scores, is.na))] <- NA_character_
  call
}

# every call the rule table can emit, its default included
karyotype_calls <- function(kc) {
  out_col <- as.character(kc[["output_column"]])
  unique(c(
    as.character(kc[["default"]]),
    vapply(kc[["rules"]], function(r) as.character(r[[out_col]]), character(1L))
  ))
}

# map pheno Female (1/0) onto the rule table labels. refuse non-0/1 values.
recorded_from_female <- function(female) {
  checkmate::assert_integerish(
    female,
    lower = 0,
    upper = 1,
    any.missing = TRUE,
    null.ok = FALSE,
    .var.name = "pheno$Female"
  )
  # na carries through: an unrecorded sex is not a disagreement
  ifelse(
    as.integer(female) == 1L,
    BINARY_CALLS[["female"]],
    BINARY_CALLS[["male"]]
  )
}

# left join recorded sex onto calls by id, never by row order.
attach_recorded <- function(out, pheno, pheno_id, pred, kc) {
  if (is.null(pheno) || !"Female" %in% names(pheno)) {
    return(out)
  }
  missing_labels <- setdiff(BINARY_CALLS, karyotype_calls(kc))
  if (length(missing_labels)) {
    catalog_bug(
      "%s: karyotype_call emits no '%s' call to compare a recorded sex against.",
      SEX_GROUP,
      missing_labels[[1L]]
    )
  }

  idx <- match(out[[pheno_id]], as.character(pheno[[pheno_id]]))
  if (anyNA(idx)) {
    stop(
      "predict_sex: a scored sample has no pheno row. This is a package bug --
       please report it.",
      call. = FALSE
    )
  }

  recorded <- recorded_from_female(pheno[["Female"]][idx])
  out[[RECORDED_SEX]] <- recorded
  # an aneuploid or unscored call is never a disagreement
  out[[SEX_MISMATCH]] <- !is.na(recorded) &
    pred %in% BINARY_CALLS &
    pred != recorded
  out
}

say_mismatch <- function(out) {
  n <- sum(out[[SEX_MISMATCH]])
  if (!n) {
    return(invisible(NULL))
  }
  cli::cli_inform(c(
    "!" = "{n} sample{?s} {cli::qty(n)}{?has/have} a predicted sex that does
           not match the {.field Female} column in {.arg pheno}.",
    "i" = "The {.field {SEX_MISMATCH}} column marks
           {cli::qty(n)}{?that sample/those samples}.",
    "i" = "A mismatch can come from the recorded sex or from the array data.
           Check both sources before you correct either one."
  ))
  invisible(NULL)
}

#' Predicted Sex Karyotype
#'
#' Predicts sex and identifies sex chromosome aneuploidy.
#'
#' @inheritParams mc-params
#' @param ... Passed to [calc_clocks()].
#'
#' @references
#' Wang Y, Hannon E, Grant OA, Gorrie-Stone TJ, Kumari M, Mill J, Zhai X,
#' McDonald-Maier KD, Schalkwyk LC (2021). DNA methylation-based sex
#' classifier to predict sex and identify sex chromosome aneuploidy.
#' *BMC Genomics*, 22(1), 484. \doi{10.1186/s12864-021-07675-2}
#'
#' @details
#' This is a re-implementation of the sex prediction algorithm of the
#' wateRmelon package.
#'
#' The returned data.frame has one row for each sample, with the two
#' `DNAmSex_Wang` scores and a `predicted_sex` column. `predicted_sex` is
#' one of `"Male"`, `"Female"`, `"47,XXY"`, or `"45,XO"`. A sample missing
#' either score gets `NA`, not a default call.
#'
#' When `pheno` has a `Female` column, coded `0` or `1`, the result also
#' carries `recorded_sex` and `sex_mismatch`. `sex_mismatch` is `TRUE` only
#' where `predicted_sex` disagrees with a binary `recorded_sex`. A
#' `"47,XXY"` or `"45,XO"` call is never flagged, because a binary `Female`
#' column cannot record it.
#'
#' @returns A data.frame. One row for each sample, with the two
#'   `DNAmSex_Wang` scores, `predicted_sex`, and, when `pheno` has a
#'   `Female` column, `recorded_sex` and `sex_mismatch`.
#'
#' @examples
#' ids <- c("DNAmSex_Wang_ChrX", "DNAmSex_Wang_ChrY")
#' sim <- sim_DNAm(ids, n = 6, Female = TRUE)
#' predict_sex(sim[["DNAm"]], sim[["pheno"]])
#'
#' @export
predict_sex <- function(DNAm, pheno = NULL, ...) {
  kc <- karyotype_spec()
  map <- karyotype_inputs(kc)

  # both members are scored together -- neither is interpretable alone
  res <- calc_clocks(DNAm, unname(map), pheno = pheno, ...)
  out <- as.data.frame(res, long = FALSE)

  scores <- lapply(map, function(id) out[[id]])
  pred <- apply_karyotype(scores, kc)
  out[[as.character(kc[["output_column"]])]] <- pred

  # the Female covariate is not required. comparison reads the caller's pheno.
  out <- attach_recorded(
    out,
    pheno,
    res[["provenance"]][["pheno_id"]],
    pred,
    kc
  )
  if (SEX_MISMATCH %in% names(out)) {
    say_mismatch(out)
  }
  out
}
