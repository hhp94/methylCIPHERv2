# partial NA -> cohort mean, fully absent -> vendor ref

# NA scan over needed columns
scan_missing_cpgs <- function(DNAm, needed_cpgs) {
  cols <- colnames(DNAm)
  present_needed <- intersect(needed_cpgs, cols)
  out <- list(
    has_na = FALSE,
    usable_cols = present_needed,
    partial_na_cols = character(0),
    all_na_cols = character(0)
  )
  if (!anyNA(DNAm)) {
    return(out)
  }

  nr <- nrow(DNAm)

  row_miss <- slideimp::mat_miss(DNAm, col = FALSE)
  dead <- rownames(DNAm)[row_miss == ncol(DNAm)]
  if (length(dead)) {
    cli::cli_abort(
      c(
        "{length(dead)} sample{?s} {cli::qty(dead)}{?has/have} no observed CpGs
         (all NA): {.val {utils::head(dead, 10L)}}.",
        "i" = "Remove or fix {cli::qty(dead)}{?it/them} before scoring."
      ),
      call = NULL
    )
  }

  sub <- DNAm[, present_needed, drop = FALSE]
  col_miss <- slideimp::mat_miss(sub, col = TRUE)
  all_na <- present_needed[col_miss == nr]
  partial <- present_needed[col_miss > 0 & col_miss < nr]

  out$has_na <- TRUE
  out$all_na_cols <- all_na
  out$partial_na_cols <- partial
  out$usable_cols <- setdiff(present_needed, all_na)
  out
}

# cohort column means for partial-NA fill
build_partial_cache <- function(DNAm, cache_cpgs) {
  if (!length(cache_cpgs)) {
    return(NULL)
  }
  sub <- DNAm[, cache_cpgs, drop = FALSE]
  slideimp::mean_imp_col(sub, cores = 1L)
}
