# run cohort parity tests
test_parity <- function(filter = "fixtures-parity", ...) {
  withr::with_envvar(
    c(MC_PARITY = "1"),
    devtools::test(filter = filter, ...)
  )
}

# accepted extensions, longest first so `.csv.gz` never reads as `.gz`
WRITE_SIM_EXTS <- c(".csv.gz", ".pq", ".duckdb")

# Write a sim cohort to disk as a chunked-front-end fixture. Canonical shape is
# cpgs x samples with the id column `cpg_id`; `transpose` writes samples x cpgs
# and names it `sample_id`. Format comes from the extension. pheno always lands
# beside it as <stem>_pheno.csv.
#
# Not an S3 generic: this file is build-ignored, so an `S3method()` entry would
# point at a function the installed package does not have -- and an
# *unregistered* method is unreachable from UseMethod even under load_all().
write_sim_DNAm <- function(
  x,
  path = file.path(tempdir(), "sim_DNAm.csv.gz"),
  transpose = FALSE,
  ...
) {
  checkmate::assert_class(x, "mc_sim")
  checkmate::assert_string(path)
  checkmate::assert_flag(transpose)

  ext <- WRITE_SIM_EXTS[endsWith(tolower(path), WRITE_SIM_EXTS)]
  if (length(ext) != 1L) {
    stop(
      sprintf(
        "`path` must end in one of %s, not \"%s\".",
        paste(WRITE_SIM_EXTS, collapse = ", "),
        path
      ),
      call. = FALSE
    )
  }
  stem <- substr(path, 1L, nchar(path) - nchar(ext))
  pheno_path <- paste0(stem, "_pheno.csv")

  # DNAm is stored samples x cpgs, so the canonical write is the transposed one
  mat <- if (transpose) x[["DNAm"]] else t(x[["DNAm"]])
  out <- cbind(rownames(mat), as.data.frame(mat, check.names = FALSE))
  names(out)[[1]] <- if (transpose) "sample_id" else "cpg_id"

  switch(
    ext,
    ".csv.gz" = {
      require_dev_ns("data.table")
      data.table::fwrite(out, path)
    },
    ".pq" = {
      require_dev_ns("arrow")
      arrow::write_parquet(out, path)
    },
    ".duckdb" = {
      require_dev_ns("duckdb")
      require_dev_ns("DBI")
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path)
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
      DBI::dbWriteTable(con, "betas", out, overwrite = TRUE)
    }
  )
  utils::write.csv(x[["pheno"]], pheno_path, row.names = FALSE)

  c(DNAm = path, pheno = pheno_path)
}

require_dev_ns <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("write_sim_DNAm() needs the %s package for this format.", pkg),
      call. = FALSE
    )
  }
}
