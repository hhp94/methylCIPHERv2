# DNAm/pheno validation and clock-id resolution

check_DNAm <- function(DNAm) {
  checkmate::assert_matrix(
    DNAm,
    mode = "double",
    min.rows = 1,
    min.cols = 1
  )
  checkmate::assert_character(colnames(DNAm), unique = TRUE, null.ok = FALSE)
  # sample ids are mandatory -- the package never invents them
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
  invisible(TRUE)
}

# Zhang2019 uses full-matrix moments
resolve_DNAm_extra <- function(clock_ids) {
  if ("Zhang2019" %in% clock_ids) {
    cli::cli_inform(c(
      "i" = "Zhang2019 takes per-sample moments over all CpGs -- a large subset
             is usually enough."
    ))
  }
  invisible(TRUE)
}

# pheno structure checks
check_pheno <- function(
  pheno,
  ID = NULL,
  extra_columns = NULL,
  sample_id = NULL
) {
  if (is.null(pheno)) {
    return(invisible(TRUE))
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
  invisible(TRUE)
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
    return(invisible(character(0)))
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
    return(invisible(character(0)))
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
  invisible(names(n_na))
}

# id column + required covariates only (row names dropped)
narrow_pheno <- function(pheno, keep) {
  out <- pheno[, keep, drop = FALSE]
  rownames(out) <- NULL
  out
}

# align pheno to sample_id by id-join, then narrow to the id column plus the
# covariates this run actually needs
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
  pheno <- pheno[match(sample_id, pheno[[pheno_id]]), , drop = FALSE]
  narrow_pheno(pheno, keep)
}

# typo-suggestion pools: matched names, recommended token values
suggestion_pools <- function() {
  routed <- sex_routed_members()
  callable <- setdiff(mc_index[["clock_id"]], names(routed$alias))
  groups <- unique(mc_index[["group_id"]])
  list(
    groups = stats::setNames(groups, groups),
    clocks = c(stats::setNames(callable, callable), routed$alias)
  )
}

# nearest pool entries for a typo (case-insensitive, substring-friendly)
did_you_mean <- function(tok, pool, n = 5L) {
  d <- utils::adist(tok, names(pool), ignore.case = TRUE, partial = TRUE)[1L, ]
  utils::head(unique(unname(pool[order(d, nchar(names(pool)))])), n)
}

# nearest-match bullets for unmatched tokens
suggestion_bullets <- function(toks, pools = suggestion_pools(), n = 5L) {
  unlist(lapply(toks, function(tok) {
    hits <- lapply(pools, function(pool) did_you_mean(tok, pool, n))
    if (length(hits) == 1L) {
      h <- hits[[1L]]
      return(c(
        "*" = cli::format_inline(
          "{.val {tok}} -- did you mean {.or {.val {h}}}?"
        )
      ))
    }
    lines <- vapply(
      names(hits),
      function(label) {
        h <- hits[[label]]
        cli::format_inline("{label}: {.or {.val {h}}}")
      },
      character(1L),
      USE.NAMES = FALSE
    )
    c(
      "*" = cli::format_inline("{.val {tok}}"),
      stats::setNames(lines, rep(" ", length(lines)))
    )
  }))
}

# user tokens -> catalog clock_ids
# precedence: "all" > tag > group_id > clock_id
resolve_clocks <- function(clocks) {
  checkmate::assert_character(
    clocks,
    min.len = 1L,
    any.missing = FALSE,
    min.chars = 1L,
    .var.name = "clocks"
  )

  members <- split(mc_index[["clock_id"]], mc_index[["group_id"]])
  clock_ids <- mc_index[["clock_id"]]

  # sex-routed members are internal -- request the alias instead
  routed <- sex_routed_members()
  asked_routed <- intersect(clocks, names(routed$alias))
  if (length(asked_routed)) {
    cli::cli_abort(
      c(
        "Can't request {length(asked_routed)} sex-specific model{?s}
         directly:",
        bullets(vapply(
          asked_routed,
          function(tok) {
            cli::format_inline(
              "{.val {tok}} -- try {.val {routed$alias[[tok]]}} instead"
            )
          },
          character(1L)
        )),
        "i" = "Sex is chosen per sample from {.arg pheno}."
      ),
      call = NULL
    )
  }
  callable <- setdiff(clock_ids, names(routed$alias))

  resolve_member <- function(tok) {
    if (!is.null(members[[tok]])) {
      return(intersect(members[[tok]], callable))
    }
    if (tok %in% callable) {
      return(tok)
    }
    NULL
  }

  resolve_one <- function(tok) {
    if (tok == "all") {
      return(callable)
    }
    tag <- MC_TAGS[[tok]]
    if (!is.null(tag)) {
      hits <- lapply(tag, resolve_member)
      dead <- tag[vapply(hits, is.null, logical(1L))]
      if (length(dead)) {
        cli::cli_abort(
          c(
            "Keyword {.val {tok}} points at {cli::qty(dead)} missing
             input{?s}: {.val {dead}}.",
            "i" = "This is a package bug -- please report it."
          ),
          call = NULL
        )
      }
      return(unique(unlist(hits, use.names = FALSE)))
    }
    resolve_member(tok)
  }

  resolved <- lapply(clocks, resolve_one)
  bad <- clocks[vapply(resolved, is.null, logical(1L))]

  if (length(bad)) {
    bad <- unique(bad)
    cli::cli_abort(
      c(
        "{length(bad)} unknown input{?s} in {.arg clocks}: {.val {bad}}.",
        "i" = "Closest matches:",
        suggestion_bullets(bad),
        "i" = "See {.fn list_clocks} or {.fn list_tags}
               ({.val {names(MC_TAGS)}})."
      ),
      call = NULL
    )
  }

  out <- unlist(resolved, use.names = FALSE)
  out[!duplicated(out)]
}

# clock_inputs closure, deps before dependents
resolve_clocks_sequence <- function(clocks) {
  st <- new.env(parent = emptyenv())
  st$out <- character(length(mc_index[["clock_id"]]))
  st$n <- 0L
  st$seen <- new.env(parent = emptyenv())

  visit <- function(id, stack) {
    if (!is.null(st$seen[[id]])) {
      return(invisible())
    }
    if (id %in% stack) {
      cycle <- c(stack[match(id, stack):length(stack)], id)
      cli::cli_abort(
        "Dependency cycle among clocks: {paste(cycle, collapse = ' -> ')}",
        call = NULL
      )
    }
    for (dep in clock_depends_on(id)) {
      visit(dep, c(stack, id))
    }
    st$n <- st$n + 1L
    st$out[[st$n]] <- id
    st$seen[[id]] <- TRUE
    invisible()
  }

  for (id in clocks) {
    visit(id, character(0))
  }
  st$out[seq_len(st$n)]
}

# collapse identical CpG panels
dedup_panels <- function(panels) {
  uniq <- list()
  idx <- integer(length(panels))
  for (i in seq_along(panels)) {
    hit <- 0L
    for (j in seq_along(uniq)) {
      if (identical(panels[[i]], uniq[[j]])) {
        hit <- j
        break
      }
    }
    if (!hit) {
      uniq[[length(uniq) + 1L]] <- panels[[i]]
      hit <- length(uniq)
    }
    idx[[i]] <- hit
  }
  list(uniq = uniq, idx = idx)
}

# per-clock normalization decision, keyed by clock id. Data-independent, so it
# is resolved once before any DNAm is touched.
resolve_normalize <- function(normalize, clock_sequence) {
  schemes <- vapply(clock_sequence, clock_norm_scheme, character(1))
  names(schemes) <- clock_sequence
  # constitutive normalization is on by default; everything else is opt-in
  out <- stats::setNames(schemes %in% NORM_CONSTITUTIVE, clock_sequence)

  if (!is.null(normalize) && length(normalize)) {
    checkmate::assert_logical(normalize, any.missing = FALSE)
    nm <- names(normalize)

    if (is.null(nm)) {
      if (length(normalize) != 1L) {
        cli::cli_abort(
          c(
            "{.arg normalize} must be one {.code TRUE}/{.code FALSE} or a
             named logical vector, got {length(normalize)} unnamed values.",
            "i" = "Name them by clock id, e.g.
                   {.code normalize = c(Horvath1 = TRUE)}."
          ),
          call = NULL
        )
      }
      # a bare policy is a wish, not a claim about any one clock: it reaches
      # the clocks that can honor it and passes over the ones that cannot
      out[clock_sequence[
        schemes %in% setdiff(NORM_SCHEMES, NORM_CONSTITUTIVE)
      ]] <- normalize
    } else {
      unknown <- setdiff(nm, clock_sequence)
      if (length(unknown)) {
        cli::cli_abort(
          c(
            "{.arg normalize} names {cli::qty(unknown)} clock{?s}
             {.val {unknown}} that {cli::qty(unknown)}{?is/are} not being
             scored.",
            "i" = "Name only clocks reached by {.arg clocks}."
          ),
          call = NULL
        )
      }
      # a request the catalog cannot express is an error; declining a scheme
      # the clock never declared is merely redundant
      unusable <- nm[normalize & !(schemes[nm] %in% NORM_SCHEMES)]
      if (length(unusable)) {
        declared <- unique(unname(schemes[unusable]))
        cli::cli_abort(
          c(
            "Cannot normalize {.val {unusable}}: {cli::qty(unusable)}
             {?it declares/they declare} {.val {declared}}.",
            "i" = "Only {.val {NORM_SCHEMES}} are expressible as a declared
                   panel plus a vendored target."
          ),
          call = NULL
        )
      }
      fixed <- nm[!normalize & schemes[nm] %in% NORM_CONSTITUTIVE]
      if (length(fixed)) {
        declared <- unique(unname(schemes[fixed]))
        cli::cli_abort(
          c(
            "Cannot decline normalization for {.val {fixed}}.",
            "i" = "{cli::qty(fixed)}{?Its/Their} {.val {declared}}
                   normalization is part of the clock definition, not
                   preprocessing."
          ),
          call = NULL
        )
      }
      out[nm] <- normalize
    }
  }

  out
}

# scoring + norm panels for the compute sequence (load packs first)
clock_panels <- function(clock_sequence, packs = NULL, normalize = NULL) {
  if (is.null(normalize)) {
    normalize <- resolve_normalize(NULL, clock_sequence)
  }
  list(
    clock_id = clock_sequence,
    score = dedup_panels(lapply(
      clock_sequence,
      clock_scoring_cpgs,
      packs = packs
    )),
    norm = dedup_panels(lapply(
      clock_sequence,
      function(cid) clock_norm_cpgs(cid, normalize[[cid]])
    ))
  )
}

# union of scoring + norm CpGs
panels_union <- function(panels) {
  unique(unlist(c(panels$score$uniq, panels$norm$uniq), use.names = FALSE))
}

# per-clock present/absent CpG sets over usable_cols
resolve_cpgs <- function(usable_cols, panels) {
  usable <- unique(usable_cols)
  clock_sequence <- panels$clock_id

  # split each distinct panel once
  split_panels <- function(d) {
    lapply(d$uniq, function(p) {
      hit <- match(p, usable, 0L) > 0L
      list(needed = p, present = p[hit], absent = p[!hit])
    })
  }
  score_parts <- split_panels(panels$score)
  norm_parts <- split_panels(panels$norm)

  per_clock <- lapply(seq_along(clock_sequence), function(i) {
    s <- score_parts[[panels$score$idx[[i]]]]
    nm <- norm_parts[[panels$norm$idx[[i]]]]
    list(
      clock_id = clock_sequence[[i]],
      score_needed = s$needed,
      score_present = s$present,
      score_absent = s$absent,
      norm_needed = nm$needed,
      norm_present = nm$present,
      norm_absent = nm$absent,
      norm_scheme = clock_norm_scheme(clock_sequence[[i]]),
      # the one declared panel fact: does this clock count over a norm panel?
      normalizes = length(nm$needed) > 0L
    )
  })
  names(per_clock) <- clock_sequence

  present_needed_union <- unique(unlist(
    lapply(c(score_parts, norm_parts), function(x) x$present),
    use.names = FALSE
  ))

  # distinct-panel parts + per-clock index (shared panels counted once)
  panel_index <- list(
    score = list(parts = score_parts, idx = panels$score$idx),
    norm = list(parts = norm_parts, idx = panels$norm$idx)
  )

  list(
    per_clock = per_clock,
    present_needed_union = present_needed_union,
    panel_index = panel_index
  )
}

# pre-score scoring-panel coverage gate
WARN_COVERAGE_MARGIN <- 1.1

# cap long failure lists
coverage_bullets <- function(lines) {
  shown <- utils::head(lines, 10L)
  if (length(lines) > length(shown)) {
    shown <- c(shown, sprintf("... and %d more", length(lines) - length(shown)))
  }
  bullets(shown)
}

check_coverage <- function(cpg_list, threshold = 0.75) {
  checkmate::assert_number(threshold, lower = 0, upper = 1)
  warn_below <- min(1, threshold * WARN_COVERAGE_MARGIN)

  panel_line <- function(id, present, needed, label) {
    sprintf(
      "%s: %d/%d %s CpGs (%.1f%%)",
      id,
      length(present),
      length(needed),
      label,
      100 * length(present) / length(needed)
    )
  }

  classify <- function(x) {
    if (!length(x$score_needed)) {
      return(list(level = "", line = NA_character_))
    }
    ratio <- length(x$score_present) / length(x$score_needed)
    level <- if (ratio == 0 || ratio < threshold) {
      "stop"
    } else if (ratio < warn_below) {
      "warn"
    } else {
      ""
    }
    list(
      level = level,
      line = panel_line(x$clock_id, x$score_present, x$score_needed, "scoring")
    )
  }

  graded <- lapply(cpg_list$per_clock, classify)
  levels <- vapply(graded, function(g) g$level, character(1L))
  lines_for <- function(lvl) {
    vapply(graded[levels == lvl], function(g) g$line, character(1L))
  }

  fail <- lines_for("stop")
  if (length(fail)) {
    cli::cli_abort(
      c(
        "{length(fail)} clock{?s} {?doesn't/don't} have enough CpGs to score
         ({.arg min_clocks_coverage} = {format(threshold)}):",
        coverage_bullets(fail),
        "i" = "Drop {cli::qty(fail)}{?it/them} from {.arg clocks}, or lower
               {.arg min_clocks_coverage}."
      ),
      call = NULL
    )
  }

  marginal <- lines_for("warn")
  if (length(marginal)) {
    cli::cli_warn(
      c(
        "{length(marginal)} clock{?s} just clear{?s/} {.arg min_clocks_coverage}
         = {format(threshold)}:",
        coverage_bullets(marginal),
        "i" = "Scores still run, but more of the panel is imputed."
      ),
      call = NULL
    )
  }

  # thin QN backgrounds warn only
  thin <- vapply(
    cpg_list$per_clock,
    function(x) {
      if (
        !length(x$norm_needed) ||
          length(x$norm_present) / length(x$norm_needed) >= threshold
      ) {
        return(NA_character_)
      }
      panel_line(x$clock_id, x$norm_present, x$norm_needed, "normalization")
    },
    character(1L)
  )
  thin <- thin[!is.na(thin)]
  if (length(thin)) {
    # the two schemes treat an absent background CpG differently
    thin_schemes <- unique(vapply(names(thin), clock_norm_scheme, character(1)))
    fate <- if (all(thin_schemes == "bmiq")) {
      "Absent background CpGs are dropped from the calibration fit."
    } else if (any(thin_schemes == "bmiq")) {
      "Absent background CpGs are dropped from a BMIQ fit, and filled from
       the reference mean for quantile normalization."
    } else {
      "Missing background CpGs are filled from the reference mean."
    }
    cli::cli_warn(
      c(
        "{length(thin)} clock{?s} {?has/have} a thin normalization background
         (under {.arg min_clocks_coverage} = {format(threshold)}):",
        coverage_bullets(thin),
        "i" = fate,
        "i" = "{.fn clocks_coverage} reports the panel counts per clock."
      ),
      call = NULL
    )
  }

  invisible(unique(c(names(levels)[levels != ""], names(thin))))
}

# per-sample observed fraction of the row-gate panel (norm if normalizes, else score)
row_coverage <- function(cov, score_miss, norm_miss) {
  if (is.null(cov)) {
    return(NULL)
  }
  qn <- isTRUE(cov[["normalizes"]])
  needed <- if (qn) cov[["norm_needed"]] else cov[["score_needed"]]
  present <- if (qn) cov[["norm_present"]] else cov[["score_present"]]
  miss <- if (qn) norm_miss else score_miss
  if (is.null(miss) || !length(needed) || needed == 0L) {
    return(NULL)
  }
  panel_ratio(present, miss, needed)
}

# per-sample coverage gate (warn only) over the hoisted coverage structure
check_row_coverage <- function(coverage, threshold = 0.75) {
  checkmate::assert_number(threshold, lower = 0, upper = 1)

  line_for <- function(id) {
    rc <- row_coverage(
      coverage$per_clock[[id]],
      coverage$sample_miss$score[[id]],
      coverage$sample_miss$norm[[id]]
    )
    if (is.null(rc)) {
      return(NA_character_)
    }
    cov <- rc[["cov"]]
    low <- !is.na(cov) & cov < threshold
    if (!any(low)) {
      return(NA_character_)
    }
    sprintf(
      "%s: %d of %d sample(s), worst %.1f%% of %d CpGs",
      id,
      sum(low),
      sum(!is.na(cov)),
      100 * min(cov[low]),
      rc[["needed"]]
    )
  }

  lines <- vapply(names(coverage$per_clock), line_for, character(1L))
  lines <- lines[!is.na(lines)]
  if (length(lines)) {
    cli::cli_warn(
      c(
        "{length(lines)} clock{?s} scored some samples under
         {.arg min_samples_coverage} = {format(threshold)}:",
        coverage_bullets(lines),
        "i" = "Those sample scores lean on imputed CpGs."
      ),
      call = NULL
    )
  }

  invisible(names(lines))
}
