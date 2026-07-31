# run cohort parity tests
test_parity <- function(filter = "fixtures-parity", ...) {
  withr::with_envvar(
    c(MC_PARITY = "1"),
    devtools::test(filter = filter, ...)
  )
}

# accepted extensions, longest first so `.csv.gz` never reads as `.gz`
WRITE_SIM_EXTS <- c(".csv.gz", ".h5")

# write an n x p u(0,1) cohort for chunked-front-end benchmarks. format from the path extension. not an s3 method (this file is build-ignored).
write_sim_DNAm <- function(
  n,
  p,
  path = file.path(tempdir(), "sim_DNAm.csv.gz"),
  transpose = FALSE,
  chunk_dims = "auto",
  gzip_level = 4L
) {
  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_count(p, positive = TRUE)
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

  cpg_id <- sprintf("cg%08d", seq_len(n))
  sample_id <- paste0("sample", seq_len(p))

  # canonical cpgs x samples; the doubles keep `n * p` off integer overflow
  mat <- matrix(
    stats::runif(as.numeric(n) * as.numeric(p)),
    nrow = n,
    dimnames = list(cpg_id, sample_id)
  )
  if (transpose) {
    mat <- t(mat)
  }

  # half 1, half 0 (odd `p` leaves the extra sample male), shuffled in one draw
  n_female <- floor(p / 2)
  female <- sample(rep(c(1, 0), c(n_female, p - n_female)))
  pheno <- data.frame(
    ID = sample_id,
    Age = stats::rnorm(p, mean = 30, sd = 5),
    Female = female
  )

  id_col <- if (transpose) "sample_id" else "cpg_id"
  switch(
    ext,
    ".csv.gz" = {
      require_dev_ns("data.table")
      # as.data.table(keep.rownames=) is one copy; cbind(as.data.frame()) is two
      data.table::fwrite(
        data.table::as.data.table(mat, keep.rownames = id_col),
        path
      )
    },
    ".h5" = {
      require_dev_ns("hdf5r")
      f <- h5_member(hdf5r::H5File, "new")(path, mode = "w")
      on.exit(h5_member(f, "close_all")(), add = TRUE)
      create <- h5_member(f, "create_dataset")
      create("betas", mat, chunk_dims = chunk_dims, gzip_level = gzip_level)
      create(id_col, rownames(mat))
      create(if (transpose) "cpg_id" else "sample_id", colnames(mat))
    }
  )
  utils::write.csv(pheno, pheno_path, row.names = FALSE)

  c(DNAm = path, pheno = pheno_path)
}

# hdf5r overloads `[[`, and `$` is banned in r/. take the r6 binding off the environment
h5_member <- function(obj, name) {
  get(name, envir = obj)
}

require_dev_ns <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("write_sim_DNAm() needs the %s package for this format.", pkg),
      call. = FALSE
    )
  }
}
