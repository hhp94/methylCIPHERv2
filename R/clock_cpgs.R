# every CpG a clocks= request needs measured -- panels plus declared moment refs

#' CpGs Required To Score Clocks
#'
#' Lists the CpGs needed to score a set of clocks, including background CpGs
#' for normalization.
#'
#' @inheritParams mc-params
#'
#' @inheritSection mc-params The assets directory
#'
#' @details
#' A clock built from other clocks also needs their CpGs. Turning
#' `normalize` on for a clock adds its background panel to the returned set.
#'
#' @returns A character vector. The CpGs needed to score `clocks`, with
#'   duplicates removed.
#'
#' @seealso
#' - [list_clocks()] for the clocks a `clocks` value accepts.
#' - [list_clock_tags()] for the tags a `clocks` value accepts.
#' - [list_mc_assets()] for the assets an external clock needs.
#'
#' @examples
#' cpgs <- clock_cpgs(c("Horvath1", "Hannum"))
#' length(cpgs)
#'
#' # normalizing Horvath1 adds its background panel to the union
#' norm_cpgs <- clock_cpgs(c("Horvath1", "Hannum"), normalize = c(Horvath1 = TRUE))
#' length(norm_cpgs)
#'
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
  sequence_cpgs(clock_sequence, packs, normalize)
}

# panels plus declared moment refs for a resolved sequence (one union).
sequence_cpgs <- function(clock_sequence, packs = NULL, normalize = NULL) {
  cpgs <- union(
    clock_panels_union(clock_sequence, packs, normalize),
    unlist(resolve_moment_domains(clock_sequence), use.names = FALSE)
  )
  cpgs[nzchar(cpgs) & !is.na(cpgs)]
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
        "{length(unresolved)} clock{?s} {cli::qty(unresolved)}{?has/have} no
         scoring CpGs: {.val {capped_vals(unresolved)}}.",
        "i" = "An external clock keeps its CpGs in an asset.",
        "i" = "Call {.fn load_mc_assets} to load the assets, or
               {.fn list_mc_assets} to see which assets are on disk."
      ),
      call = NULL
    )
  }

  # the caller filters: one blank/NA screen over the whole answer, not two
  panels_union(panels)
}
