# Shared @param definitions. Every topic that takes one of these writes
# `@inheritParams mc-params` instead of copying the text.
#
# A param belongs here only if its default is the same at every call site, or
# a topic overrides it locally. `groups` is the overriding case: the donor
# carries `Default is "all"` for the three asset verbs, and `load_mc_assets()`,
# which has no default, writes its own @param and wins.
# `x` here is an mc_result, so a topic whose `x` is something else
# (`print.mc_sim`, `cite_clocks.character`) must not inherit from this file.

#' Shared Parameters
#'
#' Parameter definitions reused across the package.
#'
#' @param DNAm A numeric matrix. The methylation beta values, with samples in
#'   the rows and CpGs in the columns.
#' @param x An `mc_result` object. The value returned by [calc_clocks()].
#' @param clocks A character vector. The clocks to score, named by clock id,
#'   group id, or tag.
#' @param pheno A data.frame. The sample metadata, with one row for each
#'   sample. Default is `NULL`.
#' @param normalize A named logical vector. Turns background normalization on
#'   for the clocks that support it. Default is `NULL`, which leaves the
#'   optional schemes off.
#' @param ext_data A string. The path to the directory that holds the clock
#'   assets. Default is `NULL`, which uses the assets directory.
#' @param ask A boolean. Asks for consent before a download or a delete.
#'   Default is `TRUE`. Pass `FALSE` to download or delete without asking, in
#'   a non-interactive session.
#' @param all_columns A boolean. Returns every column, including the ones the
#'   frame leaves out by default. Default is `FALSE`.
#' @param groups A character vector. The asset groups to act on. One or more
#'   of `"PCBrainAge"`, `"PCClocks"`, `"SystemsAge"` and `"Zhang2019"`, or
#'   `"all"` for every group. Repeated values are ignored, and an empty vector
#'   selects nothing. Default is `"all"`.
#'
#' @section The assets directory:
#' Four clock groups keep their weights in downloadable assets, outside the
#' package. `ext_data` says where to read them from, and accepts three forms.
#'
#' - `NULL` reads from the assets directory, and downloads any asset that is
#'   missing. Use [set_mc_assets_dir()] to choose that directory.
#' - A path reads only that directory, and never downloads. A missing asset
#'   is an error.
#' - Assets already in memory from [load_mc_assets()] are used directly.
#'
#' @name mc-params
#' @keywords internal
NULL
