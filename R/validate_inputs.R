# dnam/pheno structure checks, run before any clock is resolved

# probe-id prefixes across 450K / EPICv1 / EPICv2 / MSA.
PROBE_ID_PREFIXES <- c("cg", "ch", "rs", "nv")

# EPICv2/MSA probe id with address suffix (cg00002033_TC11). no match on EPICv1/450K.
PROBE_REPLICATE_SUFFIX <- "_[BT][CO][0-9]+$"

# orientation and replicate suffixes are whole-array properties. a bounded stride sample decides both.
id_sample <- function(x, n = 2000L) {
  if (!length(x)) {
    return(character(0))
  }
  if (length(x) <= n) {
    return(x)
  }
  x[unique(as.integer(seq(1L, length(x), length.out = n)))]
}

has_probe_ids <- function(x) {
  length(x) > 0L &&
    any(vapply(
      PROBE_ID_PREFIXES,
      function(p) any(startsWith(x, p)),
      logical(1)
    ))
}

check_DNAm <- function(DNAm) {
  # diagnose orientation and suffixes before the matrix refusal.
  d <- dim(DNAm)
  if (length(d) != 2L) {
    cli::cli_abort(
      "{.arg DNAm} must be two-dimensional, with samples in rows and CpGs
       in columns.",
      call = NULL
    )
  }
  cn <- id_sample(colnames(DNAm))
  rn <- id_sample(rownames(DNAm))
  cn_probes <- has_probe_ids(cn)
  rn_probes <- has_probe_ids(rn)

  # probe ids in the rows is decisive. more rows than columns is only suspicious.
  transposed <- rn_probes && !cn_probes
  if (length(cn) && !cn_probes && (transposed || d[[1L]] > d[[2L]])) {
    cli::cli_warn(
      c(
        if (transposed) {
          "{.arg DNAm} looks transposed. The probe ids are in the rows."
        } else {
          "No {.arg DNAm} column name looks like a probe id
           ({.val {PROBE_ID_PREFIXES}}), and the matrix has more rows than
           columns."
        },
        "i" = "{.fn calc_clocks} reads the samples from the rows and the CpGs
               from the columns.",
        "i" = "Use {.code t(DNAm)} if the matrix has the other orientation."
      ),
      call = NULL
    )
  }

  # EPICv2/MSA replicate probes. panels use the unsuffixed id.
  suffixed <- cn[grepl(PROBE_REPLICATE_SUFFIX, cn)]
  if (length(suffixed)) {
    cli::cli_warn(
      c(
        "{.arg DNAm} holds EPICv2 or MSA chip probes. Most clocks need those
         probes deduplicated first.",
        "i" = "An EPICv2 or MSA chip carries several probes per CpG and
               suffixes each one with its address, for example
               {.val {suffixed[[1L]]}}.",
        "i" = "A clock panel names the plain CpG id, so a suffixed column
               counts as absent.",
        "i" = "Strip the suffixes with
               {.code sub(\"_[BT][CO][0-9]+$\", \"\", colnames(DNAm))}.",
        "i" = "Stripping is the first half. Collapse the duplicate columns to
               one column per CpG as well.",
        "i" = "{.fn calc_clocks} cannot collapse them, because that step needs
               the array manifest.",
        "i" = "{.fn clock_cpgs} shows the ids a clock expects."
      ),
      call = NULL
    )
  }

  if (is.data.frame(DNAm)) {
    conv <- if (transposed) "t(as.matrix(DNAm))" else "as.matrix(DNAm)"
    cli::cli_abort(
      c(
        "{.arg DNAm} is a {.cls data.frame}. {.fn calc_clocks} needs a numeric
         {.cls matrix}.",
        "i" = "Convert the {.cls data.frame} with {.code {conv}}.",
        if (!transposed && !cn_probes) {
          c(
            "i" = "Most methylation tables hold the CpGs in the rows.",
            "i" = "Check the orientation before you convert."
          )
        }
      ),
      call = NULL
    )
  }

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
        "{.arg DNAm} needs sample ids as {.fn rownames}.",
        "i" = "The score rows carry those ids.",
        "i" = "Name unnamed rows with
               {.code rownames(DNAm) <- paste0(\"sample\", seq_len(nrow(DNAm)))}."
      ),
      call = NULL
    )
  }
  checkmate::assert_character(rownames(DNAm), unique = TRUE, null.ok = FALSE)
  invisible(NULL)
}

# say_* emits to the user. note_* records into the block's collector.
say_full_panel_clocks <- function(clock_ids) {
  full <- clock_ids[vapply(clock_ids, clock_needs_full_panel, logical(1))]
  if (!length(full)) {
    return(invisible(full))
  }
  cli::cli_inform(
    c(
      "i" = "{.val {full}} score{cli::qty(full)}{?s/} against every column of
             {.arg DNAm}, not just {cli::qty(full)}{?its/their} own panel.",
      "i" = "Pass every CpG you measured.",
      "i" = "A subset of {.arg DNAm} changes {cli::qty(full)}{?this/these}
             score{?s}."
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
  # front-door pheno structure (cli).
  if (!ID %in% names(pheno)) {
    cli::cli_abort(
      c(
        "The sample id column {.val {ID}} is not in {.arg pheno}.",
        "i" = "{.arg pheno_id} names the column that holds the sample ids.",
        "i" = "{.arg pheno} has {.field {capped_vals(names(pheno))}}."
      ),
      call = NULL
    )
  }
  # deparse names an internal arg. .var.name must be spelled out.
  checkmate::assert_character(
    pheno[[ID]],
    any.missing = FALSE,
    unique = TRUE,
    null.ok = FALSE,
    .var.name = paste0("pheno$", ID)
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
    # gated on extra_columns. units check fires only when age is consumed.
    warn_age_units(pheno, ID, sample_id)
  }
  warn_missing_covariates(pheno, ID, extra_columns, sample_id)
  invisible(NULL)
}

# age unit bounds at the edge of biological possibility.
AGE_MAX_YEARS <- 122
# pre-birth as a fraction of a year (-0.5, -1) is ok. gestational age in weeks is not.
AGE_MIN_YEARS <- -2

# rows of pheno that survive the id-join, in sample order.
joined_rows <- function(pheno, ID, sample_id) {
  if (is.null(sample_id)) {
    return(seq_len(nrow(pheno)))
  }
  idx <- match(sample_id, pheno[[ID]])
  idx[!is.na(idx)]
}

# per-row age units gate.
warn_age_units <- function(pheno, ID, sample_id) {
  rows <- joined_rows(pheno, ID, sample_id)
  age <- pheno[["Age"]][rows]
  ids <- as.character(pheno[[ID]][rows])
  # na is an unrecorded age, not a units error. Inf cannot reach here.
  ok <- !is.na(age)

  high <- ok & age > AGE_MAX_YEARS
  if (any(high)) {
    at <- which.max(replace(age, !high, -Inf))
    cli::cli_warn(
      c(
        "{sum(high)} of {length(age)} {.field Age} {cli::qty(sum(high))}value{?s}
         {?is/are} above {.val {AGE_MAX_YEARS}}.",
        "x" = "The largest is {.val {signif(age[[at]], 6)}}, for sample
               {.val {ids[[at]]}}.",
        "i" = "{.field Age} must be in years.",
        "i" = "{.val {AGE_MAX_YEARS}} is the verified human maximum, so these
               values are usually a units mistake.",
        "i" = "Convert with {.code Age / 12} for months, {.code Age / 52} for
               weeks, or {.code Age / 365.25} for days.",
        "i" = "{.fn calc_accel} reads this column, so correct it before you
               measure age acceleration."
      ),
      call = NULL
    )
  }

  low <- ok & age < AGE_MIN_YEARS
  if (any(low)) {
    at <- which.min(replace(age, !low, Inf))
    cli::cli_warn(
      c(
        "{sum(low)} of {length(age)} {.field Age} {cli::qty(sum(low))}value{?s}
         {?is/are} below {.val {AGE_MIN_YEARS}}.",
        "x" = "The smallest is {.val {signif(age[[at]], 6)}}, for sample
               {.val {ids[[at]]}}.",
        "i" = "A small negative age is normal. Some cohorts code pre-birth as
               a fraction of a year.",
        "i" = "These values are lower than that. Gestational age in weeks runs
               from {.val {-40}} to {.val {0}}.",
        "i" = "Convert weeks with {.code Age / 52}."
      ),
      call = NULL
    )
  }

  # the flagged ids, on the same axis as sample_id
  invisible(ids[high | low])
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
  rows <- joined_rows(pheno, ID, sample_id)

  n_na <- vapply(cols, function(cl) sum(is.na(pheno[[cl]][rows])), integer(1L))
  n_na <- n_na[n_na > 0L]
  if (!length(n_na)) {
    return(invisible(NULL))
  }

  cli::cli_warn(
    c(
      "{length(n_na)} {.arg pheno} covariate{?s} {?has/have} missing values.",
      capped_bullets(seq_along(n_na), function(i) {
        vapply(
          i,
          function(k) {
            cli::format_inline(
              "{.field {names(n_na)[[k]]}}: {n_na[[k]]} sample{?s}"
            )
          },
          character(1L)
        )
      }),
      "i" = "A sample with a missing covariate scores {.code NA}.",
      "i" = "Fill the column in {.arg pheno}, or drop those samples before
             you score."
    ),
    call = NULL
  )
  invisible(NULL)
}

# align pheno by id-join. none supplied -> id column alone.
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
           in {.arg DNAm}:",
          "x" = "{.val {capped_vals(missing)}}",
          "i" = "Every {.arg DNAm} row needs a matching id in the
                 {.arg pheno_id} column."
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
