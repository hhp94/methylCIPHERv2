# typo-suggestion pools: matched names, recommended token values
suggestion_pools <- function() {
  routed <- sex_routed_members()
  callable <- setdiff(mc_index[["clock_id"]], names(routed[["alias"]]))
  groups <- unique(mc_index[["group_id"]])
  list(
    groups = stats::setNames(groups, groups),
    clocks = c(stats::setNames(callable, callable), routed[["alias"]])
  )
}

# nearest pool entries for a typo (case-insensitive, substring-friendly)
did_you_mean <- function(tok, pool, n = 5L) {
  d <- utils::adist(tok, names(pool), ignore.case = TRUE, partial = TRUE)[1L, ]
  utils::head(unique(unname(pool[order(d, nchar(names(pool)))])), n)
}

# nearest-match bullets for unmatched tokens. cap is on token count.
suggestion_bullets <- function(toks, pools = suggestion_pools(), n = 5L) {
  cli_escape(unlist(lapply(capped_vals(toks), function(tok) {
    hits <- lapply(pools, function(pool) did_you_mean(tok, pool, n))
    if (length(hits) == 1L) {
      h <- hits[[1L]]
      return(c(
        "*" = cli::format_inline(
          "{.val {tok}}. Did you mean {.or {.val {h}}}?"
        )
      ))
    }
    # noun phrases, so a pool name cannot read as an argument name
    what <- c(groups = "group ids", clocks = "clock ids")
    lines <- vapply(
      names(hits),
      function(label) {
        h <- hits[[label]]
        cli::format_inline("{what[[label]]}: {.or {.val {h}}}")
      },
      character(1L),
      USE.NAMES = FALSE
    )
    c(
      "*" = cli::format_inline("{.val {tok}}"),
      stats::setNames(lines, rep(" ", length(lines)))
    )
  })))
}

# user tokens -> catalog clock_ids (all > tag > group_id > clock_id)
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

  # sex-routed members are internal -- request the alias
  routed <- sex_routed_members()
  asked_routed <- intersect(clocks, names(routed[["alias"]]))
  if (length(asked_routed)) {
    cli::cli_abort(
      c(
        "{length(asked_routed)} sex-specific model{?s} cannot be requested by
         name:",
        capped_bullets(asked_routed, function(toks) {
          vapply(
            toks,
            function(tok) {
              cli::format_inline(
                "{.val {tok}}. Request {.val {routed[['alias']][[tok]]}}
                 instead."
              )
            },
            character(1L)
          )
        }),
        "i" = "Request the family alias.",
        "i" = "The alias reads the sex of each sample from {.arg pheno}."
      ),
      call = NULL
    )
  }
  callable <- setdiff(clock_ids, names(routed[["alias"]]))

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
            "The keyword {.val {tok}} names {cli::qty(dead)} clock{?s} that
             the catalog does not contain: {.val {dead}}.",
            "i" = "This is a bug in {.pkg methylCIPHERv2}. Report it."
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
        "{length(bad)} name{?s} in {.arg clocks} {cli::qty(bad)}{?is/are} not
         a clock, a group or a keyword: {.val {capped_vals(bad)}}.",
        "i" = "Closest matches:",
        suggestion_bullets(bad),
        "i" = "See {.fn list_clocks} or {.fn list_clock_tags}
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
  st[["out"]] <- character(length(mc_index[["clock_id"]]))
  st[["n"]] <- 0L
  st[["seen"]] <- new.env(parent = emptyenv())

  visit <- function(id, stack) {
    if (!is.null(st[["seen"]][[id]])) {
      return(invisible())
    }
    if (id %in% stack) {
      cycle <- c(stack[match(id, stack):length(stack)], id)
      stop(
        sprintf(
          "Dependency cycle among clocks: %s",
          paste(cycle, collapse = " -> ")
        ),
        call. = FALSE
      )
    }
    for (dep in clock_depends_on(id)) {
      visit(dep, c(stack, id))
    }
    st[["n"]] <- st[["n"]] + 1L
    st[["out"]][[st[["n"]]]] <- id
    st[["seen"]][[id]] <- TRUE
    invisible()
  }

  for (id in clocks) {
    visit(id, character(0))
  }
  st[["out"]][seq_len(st[["n"]])]
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

# per-clock normalization decision (data-independent)
resolve_normalize <- function(normalize, clock_sequence) {
  schemes <- vapply(clock_sequence, clock_norm_scheme, character(1))
  names(schemes) <- clock_sequence
  # constitutive normalization is on by default, everything else is opt-in
  out <- stats::setNames(schemes %in% NORM_CONSTITUTIVE, clock_sequence)

  if (!is.null(normalize) && length(normalize)) {
    checkmate::assert_logical(normalize, any.missing = FALSE)
    nm <- names(normalize)

    if (is.null(nm)) {
      if (length(normalize) != 1L) {
        cli::cli_abort(
          c(
            "{.arg normalize} must be one {.code TRUE} or {.code FALSE}, or a
             named logical vector. It has {length(normalize)} unnamed values.",
            "i" = "Name them by clock id, for example
                   {.code normalize = c(Horvath1 = TRUE)}."
          ),
          call = NULL
        )
      }
      # bare policy applies where the clock can honor it
      out[clock_sequence[
        schemes %in% setdiff(NORM_SCHEMES, NORM_CONSTITUTIVE)
      ]] <- normalize
    } else {
      unknown <- setdiff(nm, clock_sequence)
      if (length(unknown)) {
        cli::cli_abort(
          c(
            "{.arg normalize} names {length(unknown)} clock{?s} that {?is/are}
             not being scored: {.val {capped_vals(unknown)}}.",
            "i" = "Name only a clock that {.arg clocks} reaches."
          ),
          call = NULL
        )
      }
      # unknown scheme is an error, declining an undeclared one is redundant
      unusable <- nm[normalize & !(schemes[nm] %in% NORM_SCHEMES)]
      if (length(unusable)) {
        declared <- unique(unname(schemes[unusable]))
        cli::cli_abort(
          c(
            "Cannot normalize {.val {capped_vals(unusable)}}.",
            "i" = "{cli::qty(declared)}The declared scheme{?s} {?is/are}
                   {.val {declared}}. Only {.val {NORM_SCHEMES}} are
                   expressible as a declared panel plus a vendored target."
          ),
          call = NULL
        )
      }
      fixed <- nm[!normalize & schemes[nm] %in% NORM_CONSTITUTIVE]
      if (length(fixed)) {
        declared <- unique(unname(schemes[fixed]))
        cli::cli_abort(
          c(
            "Cannot decline normalization for
             {.val {capped_vals(fixed)}}.",
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

# union of the named panel roles (both by default)
panels_union <- function(panels, roles = c("score", "norm")) {
  unique(unlist(
    lapply(roles, function(r) panels[[r]][["uniq"]]),
    use.names = FALSE
  ))
}

# per-clock present/absent CpG sets over usable_cols
resolve_cpgs <- function(usable_cols, panels) {
  usable <- unique(usable_cols)
  clock_sequence <- panels[["clock_id"]]

  # split each distinct panel once (factor keeps empty norm panels)
  split_panels <- function(d) {
    panels <- d[["uniq"]]
    grp <- factor(
      rep(seq_along(panels), lengths(panels)),
      levels = seq_along(panels)
    )
    hits <- split(
      match(unlist(panels, use.names = FALSE), usable, 0L) > 0L,
      grp
    )
    lapply(seq_along(panels), function(i) {
      p <- panels[[i]]
      hit <- hits[[i]]
      list(needed = p, present = p[hit], absent = p[!hit])
    })
  }
  score_parts <- split_panels(panels[["score"]])
  norm_parts <- split_panels(panels[["norm"]])

  per_clock <- lapply(seq_along(clock_sequence), function(i) {
    s <- score_parts[[panels[["score"]][["idx"]][[i]]]]
    nm <- norm_parts[[panels[["norm"]][["idx"]][[i]]]]
    list(
      clock_id = clock_sequence[[i]],
      score_needed = s[["needed"]],
      score_present = s[["present"]],
      score_absent = s[["absent"]],
      norm_needed = nm[["needed"]],
      norm_present = nm[["present"]],
      norm_absent = nm[["absent"]],
      # does this clock count over a norm panel?
      normalizes = length(nm[["needed"]]) > 0L
    )
  })
  names(per_clock) <- clock_sequence

  # distinct-panel parts + per-clock index
  panel_index <- list(
    score = list(parts = score_parts, idx = panels[["score"]][["idx"]]),
    norm = list(parts = norm_parts, idx = panels[["norm"]][["idx"]])
  )

  list(per_clock = per_clock, panel_index = panel_index)
}
