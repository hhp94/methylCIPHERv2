# n x length(cpgs) U(0,1) beta matrix over the ambient RNG (unseeded)
#' @export
random_betas <- function(cpgs, n = 10L) {
  matrix(
    stats::runif(n * length(cpgs)),
    nrow = n,
    dimnames = list(paste0("sample", seq_len(n)), cpgs)
  )
}

#' @export
sim_DNAm <- function(
  clocks,
  n = 10,
  Age = FALSE,
  Female = FALSE,
  remove = 0,
  normalize = NULL,
  ext_data = NULL,
  ask = TRUE
) {
  checkmate::assert_flag(Age)
  checkmate::assert_flag(Female)
  checkmate::assert_int(remove, lower = 0)

  clock_sequence <- resolve_clocks_sequence(resolve_clocks(clocks))
  normalize <- resolve_normalize(normalize, clock_sequence)
  packs <- load_mc_assets(pack_groups_needed(clock_sequence), ext_data, ask)
  cpgs <- clock_cpgs(clock_sequence, packs, normalize)
  if (remove > 0) {
    n_drop <- min(remove, length(cpgs))
    cpgs <- cpgs[-sample.int(length(cpgs), n_drop)]
  }
  ID <- paste0("sample", seq_len(n))
  DNAm <- random_betas(cpgs, n = n)
  pheno <- data.frame(ID = ID)
  if (Age) {
    pheno[["Age"]] <- stats::rnorm(n, mean = 45, sd = 5)
  }
  if (Female) {
    pheno[["Female"]] <- numeric(n)
    pheno[["Female"]][sample.int(n, floor(n / 2))] <- 1
  }
  out <- list(
    DNAm = DNAm,
    pheno = pheno
  )
  class(out) <- c("mc_sim", "list")
  out
}

# print mc_sim (DNAm + pheno preview)
#' @export
print.mc_sim <- function(x, n = 6, p = 6, ...) {
  DNAm <- x[["DNAm"]]
  pheno <- x[["pheno"]]
  nr <- nrow(DNAm)
  nc <- ncol(DNAm)
  ni <- min(n, nr)
  pi <- min(p, nc)

  cat(sprintf("<mc_sim> %d sample(s) x %d CpG(s)\n\n", nr, nc))
  cat(sprintf("DNAm [showing %d x %d]:\n", ni, pi))
  print(DNAm[seq_len(ni), seq_len(pi), drop = FALSE])
  if (ni < nr || pi < nc) {
    cat(sprintf("... %d more row(s), %d more col(s)\n", nr - ni, nc - pi))
  }

  cat(sprintf(
    "\npheno [showing %d of %d row(s)]:\n",
    min(n, nrow(pheno)),
    nrow(pheno)
  ))
  print(utils::head(pheno, n))

  invisible(x)
}

# CpG union a request needs -- the same panels the scorer resolves
clock_cpgs <- function(clock_ids, packs, normalize) {
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
        "No scoring CpGs for {.val {unresolved}}.",
        "i" = "External packs may need to be loaded first."
      ),
      call = NULL
    )
  }

  cpgs <- panels_union(panels)
  cpgs[nzchar(cpgs) & !is.na(cpgs)]
}
