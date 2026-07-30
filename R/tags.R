# keyword macros for clocks= (expand to groups/clocks)
# TODO: move membership upstream into methylCIPHER-meta
MC_TAGS <- list(
  gestational = c("Bohlin", "Knight", "Mayne", "LeePlacentalAge"),
  mitotic = c("EpiTOC", "EpiTOC2", "MiAge", "RepliTali"),
  mortality = c("GrimAge", "ZhangMortality")
)

# the keyword registry: tag -> the group/clock tokens it expands to.
# "all" is a token too, but not a tag -- it is every callable clock
#' @export
list_clock_tags <- function() {
  MC_TAGS
}
