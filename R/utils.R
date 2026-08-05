# shared cross-cutting helpers

# cap lists before they reach cli. one cap for every message.
MC_MSG_CAP <- 10L

# escape braces in cli bullets so data cannot become a template.
cli_escape <- function(x) {
  out <- gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
  stats::setNames(out, names(x))
}

# name a character vector as cli "*" bullets
bullets <- function(x) {
  stats::setNames(cli_escape(x), rep("*", length(x)))
}

# cli "*" bullets. cap first, then format.
capped_bullets <- function(x, fmt = identity, n = MC_MSG_CAP) {
  bullets(fmt(utils::head(x, n)))
}

# per-element {.val} markup, for capped_bullets(x, val_lines)
val_lines <- function(x) {
  vapply(x, function(v) cli::format_inline("{.val {v}}"), character(1L))
}

# capped values for the inline "{.val {capped_vals(x)}}" form
capped_vals <- function(x, n = MC_MSG_CAP) {
  utils::head(x, n)
}

# comma-joined head, for a plain stop()/warning() that is not a cli template
capped <- function(x, n = MC_MSG_CAP) {
  paste(utils::head(x, n), collapse = ", ")
}

# index that reindexes an id-keyed right side into `key`'s order. `what` names
# the call site for the defect message.
# unmatched: "stop" a package bug, "drop" keep the matches, "na" caller's own.
id_index <- function(key, id, what, unmatched = c("stop", "drop", "na")) {
  unmatched <- match.arg(unmatched)
  bug <- function(fmt, ...) {
    stop(
      sprintf(
        paste0(
          "id_index(%s): ",
          fmt,
          " This is a package bug -- please report it."
        ),
        what,
        ...
      ),
      call. = FALSE
    )
  }

  dup <- anyDuplicated(id)
  if (dup) {
    bug("the join key repeats id %s.", id[[dup]])
  }
  idx <- match(key, id)
  miss <- is.na(idx)
  if (!any(miss)) {
    return(idx)
  }
  switch(
    unmatched,
    stop = bug("%d id(s) have no row: %s.", sum(miss), capped(key[miss])),
    drop = idx[!miss],
    na = idx
  )
}

# polynomial eval, lowest degree first (horner-style, 1-row safe)
poly_eval <- function(x, coef) {
  out <- rep(0, length(x))
  for (k in seq_along(coef)) {
    out <- out + coef[[k]] * x^(k - 1L)
  }
  out
}
