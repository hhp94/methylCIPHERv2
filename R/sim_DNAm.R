# n x length(cpgs) U(0,1) beta matrix over the ambient RNG (unseeded)
random_betas <- function(cpgs, n = 10L) {
  matrix(
    stats::runif(n * length(cpgs)),
    nrow = n,
    dimnames = list(paste0("sample", seq_len(n)), cpgs)
  )
}

#' Simulated Methylation Data
#'
#' Builds a random beta matrix and a matching `pheno` data.frame for a set of
#' clocks.
#'
#' @inheritParams mc-params
#' @param n A single whole number. The number of samples to simulate. Default
#'   is `10`.
#' @param Age A boolean. Adds an `Age` column to `pheno`, drawn from a
#'   normal distribution. Default is `FALSE`.
#' @param Female A boolean. Adds a `Female` column to `pheno`, with
#'   about half the samples set to `1`. Default is `FALSE`.
#' @param remove A single whole number. The number of CpGs to drop at random
#'   from the simulated panel. Default is `0`.
#' @param suffix A string. Appended to every sample id, so two simulated
#'   matrices stay distinct. Default is `NULL`, which leaves the ids as
#'   given.
#'
#' @inheritSection mc-params The assets directory
#'
#' @returns An `mc_sim` object. It holds the simulated `DNAm` matrix, the
#'   matching `pheno` data.frame, the `clocks` argument as given, and the
#'   `suffix`, which is `NULL` when no suffix was set.
#'
#' @examples
#' sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 10, Age = TRUE, Female = TRUE)
#' dim(sim[["DNAm"]])
#' head(sim[["pheno"]])
#'
#' @export
sim_DNAm <- function(
  clocks,
  n = 10,
  Age = FALSE,
  Female = FALSE,
  remove = 0,
  normalize = NULL,
  ext_data = NULL,
  ask = TRUE,
  suffix = NULL
) {
  checkmate::assert_flag(Age)
  checkmate::assert_flag(Female)
  checkmate::assert_int(remove, lower = 0)
  # opt-in sample-id suffix.
  if (!is.null(suffix)) {
    checkmate::assert_string(suffix, min.chars = 1L)
  }

  cpgs <- clock_cpgs(clocks, normalize, ext_data, ask)
  if (remove > 0) {
    n_drop <- min(remove, length(cpgs))
    cpgs <- cpgs[-sample.int(length(cpgs), n_drop)]
  }
  # suffixed ids make two simulated blocks disjoint for rbind's first gate
  ID <- paste0("sample", seq_len(n))
  if (!is.null(suffix)) {
    ID <- paste0(ID, "_", suffix)
  }
  DNAm <- random_betas(cpgs, n = n)
  # one id source for both.
  rownames(DNAm) <- ID
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
    pheno = pheno,
    # the clocks argument as given, so a dev helper can score without restating it
    clocks = clocks,
    # null unless the ids were suffixed -- reported, never fed back in
    suffix = suffix
  )
  class(out) <- c("mc_sim", "list")
  out
}

# dnam then pheno, in the shared printer grammar (R/print.R)
#' Print Method For An mc_sim Object
#'
#' Prints a compact summary of an `mc_sim` object, with a preview of `DNAm`
#' and `pheno`.
#'
#' @param x An `mc_sim` object. The value returned by [sim_DNAm()].
#' @param n A single whole number. The number of sample rows to preview from
#'   `DNAm` and `pheno`. Default is `6`.
#' @param p A single whole number. The number of CpG columns to preview from
#'   `DNAm`. Default is `6`.
#' @param ... Not used.
#'
#' @returns An `mc_sim` object. Returns `x`, invisibly, after printing it.
#'
#' @examples
#' sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 10)
#' print(sim)
#'
#' @export
print.mc_sim <- function(x, n = 6, p = 6, ...) {
  DNAm <- x[["DNAm"]]
  pheno <- x[["pheno"]]

  cat(
    fmt_header("mc_sim", nrow(DNAm), "sample", ncol(DNAm), "CpG"),
    "\n",
    sep = ""
  )
  print_block("DNAm", DNAm, min(n, nrow(DNAm)), min(p, ncol(DNAm)), "CpG")
  print_block(
    "pheno",
    pheno,
    min(n, nrow(pheno)),
    ncol(pheno),
    "column",
    cut_cols = FALSE
  )
  # an unsuffixed sim says nothing here, so the line only appears when asked for
  if (!is.null(x[["suffix"]])) {
    cat("\n", fmt_section("suffix", x[["suffix"]]), "\n", sep = "")
  }

  invisible(x)
}
