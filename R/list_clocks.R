# browsable catalog menu for the clocks= argument

# default menu columns.
LIST_CLOCKS_DEFAULT_COLS <- c(
  "clock_id",
  "group_id",
  "request_as",
  "covariates",
  "external",
  "tags"
)

# token a user should pass to get this clock (alias for routed members)
request_token <- function(clock_id, alias) {
  ifelse(clock_id %in% names(alias), alias[clock_id], clock_id)
}

#' Epigenetic Clock Catalog
#'
#' Lists the clocks in the catalog, with the group, tags, and token to
#' request each one.
#'
#' @param group A character vector. Keeps only the clocks in these groups.
#'   Default is `NULL`, which keeps every group.
#' @param tag A character vector. Keeps only the clocks that carry one of
#'   these tags. Default is `NULL`, which applies no tag filter.
#' @param pattern A string. A regular expression matched against the clock id
#'   and the group id. Default is `NULL`, which applies no pattern filter.
#' @inheritParams mc-params
#'
#' @details
#' Valid values for `tag` are the names of [list_clock_tags()].
#'
#' `request_as` names the token to pass to `clocks` in [calc_clocks()]. It
#' differs from `clock_id` for a clock that only another clock can request.
#' `covariates` names the [calc_clocks()] `pheno` columns a clock needs, and
#' `external` is `TRUE` for a clock whose weights are a download.
#'
#' `all_columns = TRUE` adds four more columns.
#'
#' - `callable` is `FALSE` for a clock that only another clock can request.
#' - `group_size` counts the callable clocks a group token expands to.
#' - `batch_dependent` is `TRUE` for a clock whose score depends on the other
#'   samples scored with it.
#' - `normalize` names the background normalization a clock gets, and is
#'   empty for a clock that gets none. `"bmiq"` is optional, and you turn it
#'   on with the `normalize` argument of [calc_clocks()]. `"quantile"` is
#'   part of the clock's definition and is always applied.
#'
#' @returns A data.frame. One row for each clock in the catalog, including a
#'   clock that scores only as part of another clock.
#'
#' @seealso
#' - [list_clock_tags()] for the tags a `tag` value accepts.
#' - [clock_cpgs()] for the CpGs a set of clocks needs.
#' - [list_mc_assets()] for the assets an external clock needs.
#'
#' @examples
#' list_clocks(pattern = "^Horvath")
#' nrow(list_clocks(tag = "mortality"))
#' list_clocks(group = "Dunedin", all_columns = TRUE)
#'
#' @export
list_clocks <- function(
  group = NULL,
  tag = NULL,
  pattern = NULL,
  all_columns = FALSE
) {
  checkmate::assert_character(group, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_subset(tag, names(MC_TAGS), empty.ok = TRUE)
  checkmate::assert_string(pattern, null.ok = TRUE)
  checkmate::assert_flag(all_columns)

  routed <- sex_routed_members()
  idx <- mc_index
  callable <- !idx[["clock_id"]] %in% names(routed[["alias"]])

  # callable clocks a group token expands to
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
    # the scheme actually applied, so a declared but inexpressible one reads ""
    normalize = vapply(
      idx[["clock_id"]],
      function(id) {
        scheme <- clock_norm_scheme(id)
        if (scheme %in% NORM_SCHEMES) scheme else ""
      },
      character(1L),
      USE.NAMES = FALSE
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out[["group_size"]][is.na(out[["group_size"]])] <- 0L

  # tags follow request_as (routed members inherit the alias)
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
      pool <- suggestion_pools()[["groups"]]
      cli::cli_abort(
        c(
          "{length(unknown)} name{?s} in {.arg group} {cli::qty(unknown)}{?is/are}
           not a group: {.val {capped_vals(unknown)}}.",
          capped_bullets(unknown, function(toks) {
            vapply(
              toks,
              function(tok) {
                cli::format_inline(
                  "{.val {tok}}. Did you mean
                   {.or {.val {did_you_mean(tok, pool)}}}?"
                )
              },
              character(1L)
            )
          }),
          "i" = "Call {.fn list_clocks} with no arguments to see every group."
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
  if (all_columns) {
    return(out)
  }
  out[, LIST_CLOCKS_DEFAULT_COLS, drop = FALSE]
}
