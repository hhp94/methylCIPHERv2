# browsable catalog menu for the clocks= argument

# token a user should pass to get this clock (alias for routed members)
request_token <- function(clock_id, alias) {
  ifelse(clock_id %in% names(alias), alias[clock_id], clock_id)
}

#' @export
list_clocks <- function(group = NULL, tag = NULL, pattern = NULL) {
  checkmate::assert_character(group, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_subset(tag, names(MC_TAGS), empty.ok = TRUE)
  checkmate::assert_string(pattern, null.ok = TRUE)

  routed <- sex_routed_members()
  idx <- mc_index
  callable <- !idx[["clock_id"]] %in% names(routed[["alias"]])

  # how many callable clocks a group token expands to
  group_size <- table(idx[["group_id"]][callable])

  out <- data.frame(
    clock_id = idx[["clock_id"]],
    group_id = idx[["group_id"]],
    request_as = unname(request_token(idx[["clock_id"]], routed[["alias"]])),
    callable = callable,
    group_size = as.integer(group_size[idx[["group_id"]]]),
    covariates = vapply(
      idx[["covariates_required"]],
      paste,
      character(1L),
      collapse = ", "
    ),
    external = idx[["external_group"]],
    batch_dependent = idx[["batch_dependent"]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out[["group_size"]][is.na(out[["group_size"]])] <- 0L

  # tags follow request_as so routed members inherit their alias's tags
  tag_ids <- lapply(names(MC_TAGS), resolve_clocks)
  names(tag_ids) <- names(MC_TAGS)
  out[["tags"]] <- vapply(
    out[["request_as"]],
    function(id) {
      paste(
        names(MC_TAGS)[vapply(tag_ids, function(x) id %in% x, logical(1L))],
        collapse = ", "
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )

  if (!is.null(group)) {
    unknown <- setdiff(group, out[["group_id"]])
    if (length(unknown)) {
      pool <- suggestion_pools()[["group"]]
      cli::cli_abort(
        c(
          "Unknown group{cli::qty(length(unknown))}{?s}: {.val {unknown}}.",
          bullets(vapply(
            unknown,
            function(tok) {
              cli::format_inline(
                "did you mean {.or {.val {did_you_mean(tok, pool)}}}?"
              )
            },
            character(1L)
          )),
          "i" = "Run {.fn list_clocks} with no arguments to see every group."
        ),
        call = NULL
      )
    }
    out <- out[out[["group_id"]] %in% group, , drop = FALSE]
  }

  if (!is.null(tag)) {
    keep <- unique(unlist(lapply(tag, resolve_clocks), use.names = FALSE))
    out <- out[out[["request_as"]] %in% keep, , drop = FALSE]
  }

  if (!is.null(pattern)) {
    hit <- grepl(pattern, out[["clock_id"]], ignore.case = TRUE) |
      grepl(pattern, out[["group_id"]], ignore.case = TRUE)
    out <- out[hit, , drop = FALSE]
  }

  out <- out[order(out[["group_id"]], out[["clock_id"]]), , drop = FALSE]
  row.names(out) <- NULL
  out
}
