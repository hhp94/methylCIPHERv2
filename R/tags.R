# keyword macros for clocks= (expand to groups/clocks)
# TODO: move membership upstream into methylCIPHER-meta
MC_TAGS <- list(
  gestational = c("Bohlin", "Knight", "Mayne", "LeePlacentalAge"),
  mitotic = c("EpiTOC", "EpiTOC2", "MiAge", "RepliTali"),
  mortality = c("GrimAge", "ZhangMortality")
)

# keyword registry: tag -> group/clock tokens it expands to.
# "all" is a token too, but not a tag.
#' Clock Tag Registry
#'
#' Lists the keyword tags that expand to a group of clocks, for use as a
#' `clocks` or `tag` value.
#'
#' @returns A named list. Each name is a tag, and each element is a character
#'   vector of the group or clock tokens it expands to.
#'
#' @seealso
#' - [list_clocks()] for the clocks a tag expands to.
#' - [clock_cpgs()] for the CpGs a set of clocks needs.
#' - [list_mc_assets()] for the assets an external clock needs.
#'
#' @examples
#' list_clock_tags()
#'
#' @export
list_clock_tags <- function() {
  MC_TAGS
}
