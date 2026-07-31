# n x length(cpgs) U(0,1) beta matrix over the ambient RNG (unseeded)
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
  ask = TRUE,
  suffix = NULL
) {
  checkmate::assert_flag(Age)
  checkmate::assert_flag(Female)
  checkmate::assert_int(remove, lower = 0)
  # opt-in sample-id suffix. a default would silently rename every existing call
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
  # one id source for both, rather than two expressions that happen to agree
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
    # NULL unless the ids were suffixed -- reported, never fed back in
    suffix = suffix
  )
  class(out) <- c("mc_sim", "list")
  out
}

# DNAm then pheno, in the shared printer grammar (R/print.R)
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
