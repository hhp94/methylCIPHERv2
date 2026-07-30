# every CpG a clocks= request needs, over the same panels the scorer resolves

#' @export
clock_cpgs <- function(
  clocks,
  normalize = NULL,
  ext_data = NULL,
  ask = TRUE
) {
  # the sequence, not the request: a composite reads its dependencies' panels
  clock_sequence <- resolve_clocks_sequence(resolve_clocks(clocks))
  normalize <- resolve_normalize(normalize, clock_sequence)
  packs <- load_mc_assets(pack_groups_needed(clock_sequence), ext_data, ask)
  clock_panels_union(clock_sequence, packs, normalize)
}

# scoring panels plus, where a clock normalizes, its background panel
clock_panels_union <- function(clock_ids, packs, normalize) {
  panels <- clock_panels(clock_ids, packs, normalize)
  score <- panels[["score"]]
  # an empty scoring panel is fine only for a sex-routed alias (owns no panel)
  unresolved <- clock_ids[vapply(
    seq_along(clock_ids),
    function(i) {
      !length(score[["uniq"]][[score[["idx"]][[i]]]]) &&
        !length(clock_depends_on(clock_ids[[i]]))
    },
    logical(1)
  )]

  if (length(unresolved)) {
    cli::cli_abort(
      c(
        "Couldn't resolve any scoring CpGs for {.val {unresolved}}.",
        "i" = "If these are external clocks, try loading their packs first
               with {.fn load_mc_assets}."
      ),
      call = NULL
    )
  }

  cpgs <- panels_union(panels)
  cpgs[nzchar(cpgs) & !is.na(cpgs)]
}
