# run cohort parity tests
test_parity <- function(filter = "fixtures-parity", ...) {
  withr::with_envvar(
    c(MC_PARITY = "1"),
    devtools::test(filter = filter, ...)
  )
}

# score an mc_sim without unpacking it: sim_DNAm("Hannum", n = 4) |> sim_score()
sim_score <- function(sim, clocks = sim[["clocks"]], ...) {
  calc_clocks(sim[["DNAm"]], clocks, pheno = sim[["pheno"]], ...)
}

# accepted extensions, longest first so `.csv.gz` never reads as `.gz`
WRITE_SIM_EXTS <- c(".csv.gz", ".h5")

# write an n x p u(0,1) cohort for chunked-front-end benchmarks (format from path).
write_sim_DNAm <- function(
  n,
  p,
  path = file.path(tempdir(), "sim_DNAm.csv.gz"),
  transpose = FALSE,
  chunk_dims = "auto",
  gzip_level = 4L
) {
  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_count(p, positive = TRUE)
  checkmate::assert_string(path)
  checkmate::assert_flag(transpose)

  ext <- WRITE_SIM_EXTS[endsWith(tolower(path), WRITE_SIM_EXTS)]
  if (length(ext) != 1L) {
    stop(
      sprintf(
        "`path` must end in one of %s, not \"%s\".",
        paste(WRITE_SIM_EXTS, collapse = ", "),
        path
      ),
      call. = FALSE
    )
  }
  stem <- substr(path, 1L, nchar(path) - nchar(ext))
  pheno_path <- paste0(stem, "_pheno.csv")

  cpg_id <- sprintf("cg%08d", seq_len(n))
  sample_id <- paste0("sample", seq_len(p))

  # canonical cpgs x samples. doubles keep `n * p` off integer overflow.
  mat <- matrix(
    stats::runif(as.numeric(n) * as.numeric(p)),
    nrow = n,
    dimnames = list(cpg_id, sample_id)
  )
  if (transpose) {
    mat <- t(mat)
  }

  # half 1, half 0 (odd `p` leaves the extra sample male), shuffled in one draw
  n_female <- floor(p / 2)
  female <- sample(rep(c(1, 0), c(n_female, p - n_female)))
  pheno <- data.frame(
    ID = sample_id,
    Age = stats::rnorm(p, mean = 30, sd = 5),
    Female = female
  )

  id_col <- if (transpose) "sample_id" else "cpg_id"
  switch(
    ext,
    ".csv.gz" = {
      require_dev_ns("data.table")
      # as.data.table(keep.rownames=) is one copy.
      data.table::fwrite(
        data.table::as.data.table(mat, keep.rownames = id_col),
        path
      )
    },
    ".h5" = {
      require_dev_ns("hdf5r")
      f <- h5_member(hdf5r::H5File, "new")(path, mode = "w")
      on.exit(h5_member(f, "close_all")(), add = TRUE)
      create <- h5_member(f, "create_dataset")
      create("betas", mat, chunk_dims = chunk_dims, gzip_level = gzip_level)
      create(id_col, rownames(mat))
      create(if (transpose) "cpg_id" else "sample_id", colnames(mat))
    }
  )
  utils::write.csv(pheno, pheno_path, row.names = FALSE)

  c(DNAm = path, pheno = pheno_path)
}

# hdf5r overloads `[[`. take the r6 binding off the environment.
h5_member <- function(obj, name) {
  get(name, envir = obj)
}

require_dev_ns <- function(pkg, who = "write_sim_DNAm()") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("%s needs the %s package.", who, pkg),
      call. = FALSE
    )
  }
}

# --- roxygen doc lint ------------------------------------------------------
# enforces @param form from dev/WRITING.md. inherited params checked at the donor only.

# the type vocabulary. every @param opens with one of these.
DOC_TYPES <- c(
  "A numeric matrix.",
  "A character vector.",
  "A numeric vector.",
  "A string.",
  "A number between 0 and 1.",
  "A boolean.",
  "A single whole number.",
  "A data.frame.",
  "A list.",
  "A named list.",
  "An `mc_result` object.",
  "An `mc_citation` object.",
  "An `mc_sim` object.",
  "A one-sided formula.",
  "A named logical vector.",
  "One of ",
  # a catch-all S3 method argument, and a method that only ever throws
  "Any object.",
  "Nothing.",
  # the three fixed sentences for `...`
  "Passed to ",
  "Not used.",
  "Two or more "
)

# @returns names a type somewhere in its text (not held to the opening-fragment form).
DOC_TYPE_STEMS <- sub("\\.$", "", DOC_TYPES)

DOC_DEFAULT <- "Default is "

# returns a data.frame of violations, zero rows when clean.
lint_roxygen <- function(path = ".") {
  require_dev_ns("roxygen2", "lint_roxygen()")
  blocks <- roxygen2::parse_package(path, env = NULL)
  found <- do.call(rbind, lapply(blocks, lint_doc_block))
  if (is.null(found)) {
    found <- data.frame(
      topic = character(0),
      param = character(0),
      problem = character(0)
    )
  }
  rownames(found) <- NULL
  found
}

lint_doc_block <- function(block) {
  tags <- roxygen2::block_get_tags(block, "param")
  returns <- c(
    roxygen2::block_get_tags(block, "returns"),
    roxygen2::block_get_tags(block, "return")
  )
  if (!length(tags) && !length(returns)) {
    return(NULL)
  }
  topic <- doc_block_name(block)
  defaults <- doc_block_defaults(block)
  rows <- lapply(tags, function(tag) {
    param <- tag[["val"]][["name"]]
    problems <- doc_param_problems(
      param,
      doc_squish(tag[["val"]][["description"]]),
      defaults
    )
    if (!length(problems)) {
      return(NULL)
    }
    data.frame(topic = topic, param = param, problem = problems)
  })
  do.call(rbind, c(rows, list(lint_doc_returns(topic, returns))))
}

# @returns must name a type somewhere in its text.
lint_doc_returns <- function(topic, returns) {
  if (!length(returns)) {
    return(NULL)
  }
  text <- doc_squish(returns[[1L]][["val"]])
  if (doc_has_type(text)) {
    return(NULL)
  }
  data.frame(
    topic = topic,
    param = "@returns",
    problem = "does not name a type from DOC_TYPES"
  )
}

doc_has_type <- function(text) {
  any(vapply(
    DOC_TYPE_STEMS,
    function(stem) grepl(stem, text, fixed = TRUE),
    logical(1)
  ))
}

doc_param_problems <- function(param, text, defaults) {
  out <- character(0)
  if (!any(startsWith(text, DOC_TYPES))) {
    out <- c(out, "does not open with a type from DOC_TYPES")
  }
  # a block with no formals (the donor) states its own defaults
  if (!param %in% names(defaults)) {
    return(out)
  }
  says <- grepl(DOC_DEFAULT, text, fixed = TRUE)
  if (defaults[[param]] && !says) {
    out <- c(out, "has a default but does not state it")
  }
  if (!defaults[[param]] && says) {
    out <- c(out, "states a default it does not have")
  }
  out
}

# the assigned name, or the @name tag for a NULL block like mc-params
doc_block_name <- function(block) {
  named <- roxygen2::block_get_tag_value(block, "name")
  if (!is.null(named)) {
    return(named)
  }
  call <- block[["call"]]
  ok <- is.call(call) &&
    length(call) >= 3L &&
    is.name(call[[1L]]) &&
    as.character(call[[1L]]) %in% c("<-", "=")
  if (ok) as.character(call[[2L]]) else NA_character_
}

# one flag per formal: does it carry a default? empty for a NULL block.
doc_block_defaults <- function(block) {
  call <- block[["call"]]
  ok <- is.call(call) &&
    length(call) >= 3L &&
    is.call(call[[3L]]) &&
    identical(call[[3L]][[1L]], quote(`function`))
  if (!ok) {
    return(logical(0))
  }
  args <- as.list(call[[3L]][[2L]])
  args <- args[names(args) != "..."]
  vapply(
    args,
    function(v) !(is.symbol(v) && !nzchar(as.character(v))),
    logical(1)
  )
}

# roxygen wraps descriptions, so collapse the newlines before matching
doc_squish <- function(x) {
  gsub("\\s+", " ", trimws(x))
}

# --- @seealso cross-reference lint -----------------------------------------
# closed @seealso groups from dev/WRITING.md. every link must be two-way. reads man/.

# every \link target under a node. the target is in Rd_option for \link[=x]{y}.
rd_link_targets <- function(node) {
  if (identical(attr(node, "Rd_tag"), "\\link")) {
    opt <- attr(node, "Rd_option")
    tgt <- if (is.null(opt)) {
      paste(unlist(node), collapse = "")
    } else {
      paste(unlist(opt), collapse = "")
    }
    return(sub("^=", "", tgt))
  }
  if (is.list(node)) {
    return(unlist(lapply(node, rd_link_targets), use.names = FALSE))
  }
  character(0)
}

# tag -> the text of every node carrying it, for one parsed Rd
rd_nodes <- function(db, tag) {
  db[vapply(db, function(e) identical(attr(e, "Rd_tag"), tag), logical(1))]
}

# returns a data.frame of violations, zero rows when clean.
lint_seealso <- function(path = ".") {
  rds <- list.files(
    file.path(path, "man"),
    pattern = "\\.Rd$",
    full.names = TRUE
  )
  alias <- list()
  links <- list()
  for (rd in rds) {
    db <- tools::parse_Rd(rd)
    nm <- trimws(paste(unlist(rd_nodes(db, "\\name")), collapse = ""))
    alias[[nm]] <- trimws(vapply(
      rd_nodes(db, "\\alias"),
      function(e) paste(unlist(e), collapse = ""),
      character(1)
    ))
    see <- rd_nodes(db, "\\seealso")
    if (length(see)) {
      links[[nm]] <- unique(rd_link_targets(see[[1L]]))
    }
  }
  # a link target is an alias, so resolve it to the topic that owns it
  owner <- function(a) {
    a <- sub("\\(\\)$", "", a)
    hit <- names(alias)[vapply(alias, function(x) a %in% x, logical(1))]
    if (length(hit)) hit[[1L]] else NA_character_
  }

  found <- list()
  for (nm in names(links)) {
    for (t in links[[nm]]) {
      o <- owner(t)
      if (is.na(o)) {
        # an R CMD check WARNING, and the only one this catches without check
        found[[length(found) + 1L]] <- data.frame(
          topic = nm,
          target = t,
          problem = "link has no topic"
        )
        next
      }
      back <- links[[o]]
      two_way <- !is.null(back) &&
        any(vapply(back, function(b) identical(owner(b), nm), logical(1)))
      if (!two_way) {
        found[[length(found) + 1L]] <- data.frame(
          topic = nm,
          target = o,
          problem = "link is one-way"
        )
      }
    }
  }
  out <- do.call(rbind, found)
  if (is.null(out)) {
    out <- data.frame(
      topic = character(0),
      target = character(0),
      problem = character(0)
    )
  }
  rownames(out) <- NULL
  out
}

# render every man page to one text file (Rd2txt console form).
dump_roxygen <- function(
  path = ".",
  out = file.path(tempdir(), "roxygen-dump.txt")
) {
  rds <- sort(list.files(
    file.path(path, "man"),
    pattern = "\\.Rd$",
    full.names = TRUE
  ))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  for (rd in rds) {
    writeLines(c(strrep("=", 72), basename(rd), strrep("=", 72)), con)
    tools::Rd2txt(rd, out = con, options = list(underline_titles = FALSE))
    writeLines("", con)
  }
  out
}
