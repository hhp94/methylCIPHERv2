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
        "{.arg DNAm} needs sample ids as rownames so scores can be matched
         to samples.",
        "i" = "If the rows are anonymous, you can name them with:
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
        "i" = "{.fn calc_clocks} expects samples in rows and CpGs in columns.
               Try {.code t(DNAm)} if yours is the other way around."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

# note when a clock scores against the whole matrix, not its panel.
# returns the full-panel ids so callers need not re-sweep the sequence
note_full_panel_clocks <- function(clock_ids) {
  full <- clock_ids[vapply(clock_ids, clock_needs_full_panel, logical(1))]
  if (!length(full)) {
    return(invisible(full))
  }
  cli::cli_inform(
    c(
      "i" = "{.val {full}} score{cli::qty(full)}{?s/} against every column of
             {.arg DNAm}, not just {cli::qty(full)}{?its/their} own panel.",
      "i" = "Pass every CpG you measured -- a pre-subset {.arg DNAm} changes
             {cli::qty(full)}{?this/these} score{?s}."
    )
  )
  invisible(full)
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
  checkmate::assert_choice(ID, names(pheno))
  checkmate::assert_character(
    pheno[[ID]],
    any.missing = FALSE,
    unique = TRUE,
    null.ok = FALSE
  )
  # required covariates must exist -- the score branches read them unguarded
  miss <- setdiff(extra_columns, names(pheno))
  if (length(miss)) {
    cli::cli_abort(
      c(
        "{.arg pheno} is missing {length(miss)} column{?s} the requested
         clocks need: {.field {miss}}.",
        "i" = "Add {cli::qty(miss)}{?it/them} to {.arg pheno} and try again."
      ),
      call = NULL
    )
  }
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
      "i" = "Clocks that need them will score NA for those samples."
    ),
    call = NULL
  )
  invisible(NULL)
}

# align pheno by id-join. with none supplied it is the id column alone
resolve_pheno <- function(DNAm, pheno, pheno_id, keep) {
  sample_id <- rownames(DNAm)
  out <- if (is.null(pheno)) {
    stats::setNames(data.frame(sample_id, stringsAsFactors = FALSE), pheno_id)
  } else {
    missing <- setdiff(sample_id, pheno[[pheno_id]])
    if (length(missing)) {
      cli::cli_abort(
        c(
          "{.arg pheno} is missing {length(missing)} sample id{?s} that appear
           in DNAm:",
          "x" = "{.val {utils::head(missing, 10L)}}",
          "i" = "Every DNAm row needs a matching id in the pheno id column."
        ),
        call = NULL
      )
    }
    # id column + required covariates only
    pheno[
      match(sample_id, pheno[[pheno_id]]),
      unique(c(pheno_id, keep)),
      drop = FALSE
    ]
  }
  # keyed by the id column, never by row names
  rownames(out) <- NULL
  out
}
