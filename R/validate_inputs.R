# DNAm/pheno structure checks, run before any clock is resolved


check_DNAm <- function(DNAm) {
  checkmate::assert_matrix(
    DNAm,
    mode = "double",
    min.rows = 1,
    min.cols = 1
  )
  checkmate::assert_character(colnames(DNAm), unique = TRUE, null.ok = FALSE)
  # sample ids are mandatory
  if (is.null(rownames(DNAm))) {
    cli::cli_abort(
      c(
        "{.arg DNAm} needs sample ids as rownames.",
        "i" = "If the rows are anonymous, name them yourself:
               {.code rownames(DNAm) <- paste0(\"sample\", seq_len(nrow(DNAm)))}"
      ),
      call = NULL
    )
  }
  checkmate::assert_character(rownames(DNAm), unique = TRUE, null.ok = FALSE)
  # cg... ids should be columns
  if (ncol(DNAm) < 2e5 && !any(startsWith(colnames(DNAm), "cg"))) {
    cli::cli_warn(
      c(
        if (any(startsWith(rownames(DNAm), "cg"))) {
          "DNAm looks transposed -- CpG ids (cg...) are in the rows."
        } else {
          "No DNAm column names look like CpG ids (cg...)."
        },
        "i" = "{.fn calc_clocks} wants samples in rows and CpGs in columns.
               Try {.code t(DNAm)} if yours is the other way around."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# note when a clock scores against the whole matrix, not its panel
note_full_panel_clocks <- function(clock_ids) {
  full <- clock_ids[vapply(clock_ids, clock_needs_full_panel, logical(1))]
  if (!length(full)) {
    return(invisible(NULL))
  }
  cli::cli_inform(c(
    "i" = "{.val {full}} take{?s/} per-sample moments over all CpGs -- a large
           subset is usually enough."
  ))
  invisible(NULL)
}

# pheno structure checks
check_pheno <- function(
  pheno,
  ID = NULL,
  extra_columns = NULL,
  sample_id = NULL
) {
  if (is.null(pheno)) {
    return(invisible(NULL))
  }
  checkmate::assert_data_frame(pheno, min.rows = 1)
  checkmate::assert_string(ID, null.ok = FALSE)
  checkmate::assert_choice(ID, names(pheno))
  checkmate::assert_character(
    pheno[[ID]],
    any.missing = FALSE,
    unique = TRUE,
    null.ok = FALSE
  )
  if ("Female" %in% extra_columns) {
    checkmate::assert_integerish(
      pheno[["Female"]],
      lower = 0,
      upper = 1,
      null.ok = FALSE,
      any.missing = TRUE
    )
  }
  if ("Age" %in% extra_columns) {
    checkmate::assert_numeric(
      pheno[["Age"]],
      finite = TRUE,
      null.ok = FALSE,
      any.missing = TRUE
    )
  }
  warn_missing_covariates(pheno, ID, extra_columns, sample_id)
  invisible(NULL)
}

# warn on NA in required covariates
warn_missing_covariates <- function(
  pheno,
  ID,
  extra_columns,
  sample_id
) {
  cols <- intersect(extra_columns, names(pheno))
  if (!length(cols)) {
    return(invisible(NULL))
  }
  # only rows that survive the id-join
  rows <- if (is.null(sample_id)) {
    seq_len(nrow(pheno))
  } else {
    idx <- match(sample_id, pheno[[ID]])
    idx[!is.na(idx)]
  }

  n_na <- vapply(cols, function(cl) sum(is.na(pheno[[cl]][rows])), integer(1L))
  n_na <- n_na[n_na > 0L]
  if (!length(n_na)) {
    return(invisible(NULL))
  }

  cli::cli_warn(
    c(
      "Missing values in {length(n_na)} pheno covariate{?s}:",
      bullets(vapply(
        seq_along(n_na),
        function(i) {
          cli::format_inline(
            "{.field {names(n_na)[[i]]}}: {n_na[[i]]} sample{?s}"
          )
        },
        character(1L)
      )),
      "i" = "Clocks that need them score NA for those samples."
    ),
    call = NULL
  )
  invisible(NULL)
}

# align pheno by id-join, keep id column + required covariates
resolve_pheno <- function(DNAm, pheno, pheno_id, keep) {
  if (is.null(pheno)) {
    return(NULL)
  }
  sample_id <- rownames(DNAm)
  keep <- unique(c(pheno_id, keep))

  missing <- setdiff(sample_id, pheno[[pheno_id]])
  if (length(missing)) {
    cli::cli_abort(
      c(
        "pheno is missing {length(missing)} sample id{?s} from DNAm:",
        "x" = "{.val {utils::head(missing, 10L)}}"
      ),
      call = NULL
    )
  }
  # id column + required covariates only
  out <- pheno[match(sample_id, pheno[[pheno_id]]), keep, drop = FALSE]
  rownames(out) <- NULL
  out
}
