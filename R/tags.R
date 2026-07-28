# keyword macros for clocks= (expand to groups/clocks)
# TODO: move membership upstream into methylCIPHER-meta
MC_TAGS <- list(
  gestational = c("Bohlin", "Knight", "Mayne", "LeePlacentalAge"),
  mitotic = c("EpiTOC", "EpiTOC2", "MiAge", "RepliTali"),
  mortality = c("GrimAge", "ZhangMortality")
)

# print the keyword registry (returns it invisibly)
#' @export
list_tags <- function() {
  cli::cli_text("Keywords you can pass to {.arg clocks}:")
  for (tag in names(MC_TAGS)) {
    ids <- resolve_clocks(tag)
    cli::cli_bullets(c(
      "*" = cli::format_inline(
        "{.strong {tag}} ({length(ids)} clock{?s}): {.val {MC_TAGS[[tag]]}}"
      )
    ))
  }
  cli::cli_bullets(c(
    "i" = "{.code \"all\"} scores every callable clock.
           {.fn list_clocks} lists them one by one."
  ))
  invisible(MC_TAGS)
}
