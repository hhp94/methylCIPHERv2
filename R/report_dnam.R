# DNAm-input QC for report(): format check, array guess, beta-range sanity,
# missingness, and per-clock probe coverage over the raw matrix. Computations
# only -- they return a structured list that report_render.R draws.

# structural gate: abort on input report() cannot assess, collect softer issues
report_check_dnam <- function(DNAm) {
  if (!is.matrix(DNAm) || !is.numeric(DNAm)) {
    cli::cli_abort(
      c(
        "{.arg DNAm} must be a numeric matrix (samples x CpGs).",
        "i" = "report() cannot assess an object of class {.cls {class(DNAm)}}."
      ),
      call = NULL
    )
  }
  if (is.null(rownames(DNAm))) {
    cli::cli_abort(
      c(
        "{.arg DNAm} has no sample ids in its rownames.",
        "i" = "Name the rows (samples) before running report()."
      ),
      call = NULL
    )
  }
  if (is.null(colnames(DNAm))) {
    cli::cli_abort(
      c(
        "{.arg DNAm} has no CpG ids in its colnames.",
        "i" = "report() needs CpG ids (cg...) as column names."
      ),
      call = NULL
    )
  }

  issues <- character(0)
  cols_cg <- any(startsWith(colnames(DNAm), "cg"))
  rows_cg <- any(startsWith(rownames(DNAm), "cg"))
  transposed <- !cols_cg && rows_cg
  if (transposed) {
    issues <- c(
      issues,
      "CpG ids (cg...) are in the rows -- the matrix looks transposed
       (want samples x CpGs; try t(DNAm))."
    )
  } else if (!cols_cg) {
    issues <- c(issues, "No column names look like CpG ids (cg...).")
  }
  if (anyDuplicated(colnames(DNAm))) {
    issues <- c(issues, "Duplicate CpG ids among the columns.")
  }
  if (anyDuplicated(rownames(DNAm))) {
    issues <- c(issues, "Duplicate sample ids among the rows.")
  }

  list(
    ok = length(issues) == 0L,
    issues = issues,
    n_samples = nrow(DNAm),
    n_probes = ncol(DNAm),
    transposed = transposed
  )
}

# approximate cg-probe counts per Illumina methylation array (manifest-version
# dependent -- used only as a nearest-size guess, never as a hard identity)
MC_ARRAY_SIZES <- c(
  `27K` = 27578L,
  `450K` = 485512L,
  EPICv1 = 865859L,
  EPICv2 = 937690L
)

# best-guess array from cg-probe count, with an EPICv2 replicate-suffix override
detect_array <- function(DNAm) {
  cols <- colnames(DNAm)
  n_cg <- sum(startsWith(cols, "cg"))
  n <- if (n_cg > 0L) n_cg else ncol(DNAm)
  # EPICv2 relabels replicate probes as cg#######_<two letters><two digits>
  n_v2 <- sum(grepl("_[A-Za-z]{2}[0-9]{2}$", cols))

  d <- abs(n - MC_ARRAY_SIZES)
  nearest <- names(MC_ARRAY_SIZES)[which.min(d)]
  pct <- n / MC_ARRAY_SIZES[[nearest]]
  confident <- abs(1 - pct) < 0.05

  guess <- if (n_v2 >= 10L) {
    "EPICv2"
  } else if (confident) {
    nearest
  } else {
    "subset/unknown"
  }
  note <- if (identical(guess, "subset/unknown")) {
    sprintf(
      "%d probes is %.0f%% of the closest full array (%s) -- likely a filtered
       or targeted subset, not a whole array.",
      n,
      100 * pct,
      nearest
    )
  } else if (n_v2 >= 10L && !identical(nearest, "EPICv2")) {
    sprintf("EPICv2 replicate-probe suffixes detected (%d probes).", n_v2)
  } else {
    sprintf("%d cg probes match %s (%.0f%% of its full size).", n, nearest, 100 * pct)
  }

  list(
    guess = guess,
    n_probes = ncol(DNAm),
    n_cg = n_cg,
    nearest = nearest,
    pct_of_nearest = as.numeric(pct),
    epicv2_suffixes = n_v2,
    note = note
  )
}

# value-range and scale sanity over the whole matrix
check_beta_range <- function(DNAm) {
  rng <- suppressWarnings(range(DNAm, na.rm = TRUE))
  n_below <- sum(DNAm < 0, na.rm = TRUE)
  n_above <- sum(DNAm > 1, na.rm = TRUE)
  n_nan <- sum(is.nan(DNAm))
  n_inf <- sum(is.infinite(DNAm))

  # offending probes via per-column ranges (memory-light on large arrays)
  cr <- matrixStats::colRanges(DNAm, na.rm = TRUE)
  oor <- which(cr[, 1L] < 0 | cr[, 2L] > 1)
  oor_probes <- colnames(DNAm)[oor]

  scale_note <- if (!is.finite(rng[1L]) || !is.finite(rng[2L])) {
    "All values are NA or non-finite."
  } else if (rng[1L] >= 0 && rng[2L] <= 1) {
    "Values lie within [0, 1] -- consistent with beta values."
  } else if (rng[1L] < 0 && rng[2L] > 4) {
    "Negative values with a wide range -- these look like M-values, not betas."
  } else if (rng[1L] >= 0 && rng[2L] > 50) {
    "Values run well above 1 -- these look like percentages; divide by 100."
  } else {
    "Some values fall outside [0, 1] -- check for data errors or the wrong scale."
  }

  list(
    min = rng[1L],
    max = rng[2L],
    n_below0 = as.integer(n_below),
    n_above1 = as.integer(n_above),
    n_nan = as.integer(n_nan),
    n_inf = as.integer(n_inf),
    n_out_of_range = length(oor),
    out_of_range_probes = oor_probes,
    scale_note = scale_note
  )
}

# external clock ids whose group is in `groups` (i.e. cannot be assessed)
external_ids_in_groups <- function(ids, groups) {
  if (!length(groups)) {
    return(character(0))
  }
  ids[vapply(
    ids,
    function(id) clock_is_external(id) && clock_group_id(id) %in% groups,
    logical(1L)
  )]
}

# load external packs only when they can be resolved without a surprise download:
# honor an explicit `assets`, else use whatever is already cached and skip the rest
report_load_packs <- function(clock_sequence, assets, ask) {
  needed <- pack_groups_needed(clock_sequence)
  if (!length(needed)) {
    return(list(packs = NULL, skipped = character(0)))
  }
  if (!is.null(assets)) {
    return(list(packs = load_mc_assets(needed, assets, ask), skipped = character(0)))
  }
  cached <- names(mc_cached_files(needed))
  skipped <- setdiff(needed, cached)
  packs <- if (length(cached)) load_mc_assets(cached, NULL, ask) else NULL
  list(packs = packs, skipped = skipped)
}

# evenly-spaced column subsample for the O(samples^2) / row-scan stats (bounds
# cost on full arrays without a seed; deterministic, unlike random sampling)
rp_subsample_cols <- function(sub, cap = 20000L) {
  if (ncol(sub) <= cap) {
    return(sub)
  }
  sub[, round(seq(1L, ncol(sub), length.out = cap)), drop = FALSE]
}

# MAD-outlier flag (robust); FALSE for everything when spread is degenerate
rp_mad_outlier <- function(x, k = 4) {
  m <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(FALSE, length(x)))
  }
  is.finite(x) & abs(x - m) > k * s
}

# suspiciously similar sample pairs (near-identical -> swap / replicate)
detect_duplicates <- function(cols, ids, threshold = 0.999, max_n = 2000L) {
  n <- nrow(cols)
  empty <- data.frame(
    a = character(0),
    b = character(0),
    cor = numeric(0),
    stringsAsFactors = FALSE
  )
  if (n < 2L || n > max_n) {
    return(empty)
  }
  cm <- suppressWarnings(stats::cor(t(cols), use = "pairwise.complete.obs"))
  pr <- which(upper.tri(cm) & cm >= threshold, arr.ind = TRUE)
  if (!nrow(pr)) {
    return(empty)
  }
  data.frame(
    a = ids[pr[, 1L]],
    b = ids[pr[, 2L]],
    cor = cm[cbind(pr[, 1L], pr[, 2L])],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# per-sample missingness, distribution shape, outliers, and duplicate pairs
report_samples <- function(DNAm) {
  ids <- rownames(DNAm)

  # data-quality QC is a property of the whole array, not the clock subset:
  # clock CpGs are age-selected and less sharply bimodal, so scoping the
  # distribution to them mis-flags a genuinely bimodal dataset.
  row_miss <- as.integer(slideimp::mat_miss(DNAm, col = FALSE))
  na_frac <- row_miss / max(1L, ncol(DNAm))

  cols <- rp_subsample_cols(DNAm)
  mean_beta <- matrixStats::rowMeans2(cols, na.rm = TRUE)
  median_beta <- matrixStats::rowMedians(cols, na.rm = TRUE)
  # betas are bimodal (mass near 0 and 1); a high middle fraction is a flat /
  # unimodal sample -- a bad or wrongly-scaled distribution
  mid_frac <- matrixStats::rowMeans2(cols > 0.2 & cols < 0.8, na.rm = TRUE)

  # bimodality is judged from the density *shape* -- a real trough between two
  # extreme modes -- not from the middle-band fraction, which real methylation
  # (lots of intermediate CpGs) drives to ~0.5 even when clearly bimodal.
  dens <- sample_densities(cols)
  bimodal_ratio <- rep(NA_real_, length(ids))
  bimodal_ratio[dens$idx] <- dens$ratio
  cohort_ratio <- stats::median(dens$ratio, na.rm = TRUE)
  cohort_bimodal_ok <- is.finite(cohort_ratio) && cohort_ratio < BIMODAL_RATIO_MAX

  # per-sample flags are all *relative* to the cohort (a sample unlike its
  # peers). A whole cohort that is not bimodal is the separate one-off
  # `cohort_bimodal_ok` signal -- flagging every sample would be noise.
  hi_na <- na_frac > 0.1 | rp_mad_outlier(na_frac)
  mean_out <- rp_mad_outlier(mean_beta)
  # only call a single sample non-bimodal when the cohort otherwise is (odd one
  # out); if the whole cohort is flat, the cohort flag already says so
  distro_out <- if (cohort_bimodal_ok) {
    is.finite(bimodal_ratio) & bimodal_ratio >= BIMODAL_RATIO_MAX
  } else {
    rep(FALSE, length(ids))
  }

  flags <- vapply(
    seq_along(ids),
    function(i) {
      paste(c(
        if (hi_na[i]) "high-NA",
        if (mean_out[i]) "mean-outlier",
        if (distro_out[i]) "distro-outlier"
      ), collapse = ",")
    },
    character(1L)
  )

  table <- data.frame(
    sample_id = ids,
    na_frac = na_frac,
    mean_beta = mean_beta,
    median_beta = median_beta,
    mid_frac = mid_frac,
    flags = flags,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    table = table,
    n_flagged = sum(nzchar(flags)),
    cohort_bimodal_ratio = cohort_ratio,
    cohort_bimodal_ok = cohort_bimodal_ok,
    duplicates = detect_duplicates(cols, ids),
    densities = list(x = dens$x, y = dens$y, flagged = distro_out[dens$idx])
  )
}

# a healthy beta density is bimodal: the middle trough sits well below the two
# extreme modes. This ratio (center density / smaller flanking mode) is ~0.1 for
# real methylation and ~1 for flat/uniform/mis-scaled data. Flag at >= this.
BIMODAL_RATIO_MAX <- 0.6

# center-to-mode ratio of one density curve on `grid`; NA if unusable
bimodality_ratio <- function(y, grid) {
  if (all(!is.finite(y))) {
    return(NA_real_)
  }
  low <- suppressWarnings(max(y[grid <= 0.35], na.rm = TRUE))
  high <- suppressWarnings(max(y[grid >= 0.65], na.rm = TRUE))
  center <- mean(y[grid >= 0.4 & grid <= 0.6], na.rm = TRUE)
  peaks <- min(low, high)
  if (!is.finite(peaks) || peaks <= 0 || !is.finite(center)) {
    return(NA_real_)
  }
  center / peaks
}

# compact per-sample beta density curves on a shared [0,1] grid, capped in count
# so the stored object stays small, plus each curve's bimodality ratio. `idx` is
# which samples (rows of `cols`) were used.
sample_densities <- function(cols, grid_n = 128L, max_curves = 120L) {
  grid <- seq(0, 1, length.out = grid_n)
  n <- nrow(cols)
  idx <- if (n > max_curves) round(seq(1L, n, length.out = max_curves)) else seq_len(n)
  y <- vapply(
    idx,
    function(i) {
      x <- cols[i, ]
      x <- x[is.finite(x)]
      if (length(x) < 10L) {
        return(rep(NA_real_, grid_n))
      }
      stats::density(x, from = 0, to = 1, n = grid_n)$y
    },
    numeric(grid_n)
  )
  ratio <- apply(y, 2L, bimodality_ratio, grid = grid)
  list(x = grid, y = y, idx = idx, ratio = ratio)
}

# one row per assessed clock: needed / present / absent-from-matrix / all-NA / coverage
build_coverage_table <- function(cpg_list, cols, mna, seq_ids, output_ids) {
  per <- cpg_list$per_clock
  rows <- lapply(seq_ids, function(id) {
    needed <- per[[id]][["score_needed"]]
    in_matrix <- intersect(needed, cols)
    absent <- setdiff(needed, cols)
    all_na <- intersect(in_matrix, mna$all_na_cols)
    usable <- setdiff(in_matrix, mna$all_na_cols)
    data.frame(
      clock_id = id,
      group_id = clock_group_id(id),
      external = clock_is_external(id),
      role = if (id %in% output_ids) "returned" else "dependency",
      needed = length(needed),
      present = length(usable),
      absent = length(absent),
      all_na = length(all_na),
      coverage = if (length(needed)) length(usable) / length(needed) else NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  # clocks with no scoring panel of their own (sex-routed aliases) carry NA needed
  out <- out[out$needed > 0L, , drop = FALSE]
  out[order(out$coverage, out$clock_id), , drop = FALSE]
}

# assemble the DNAm-side report list
report_dnam <- function(DNAm, clocks = "all", assets = NULL, ask = TRUE) {
  DNAm <- coerce_dnam(DNAm)
  fmt <- report_check_dnam(DNAm)
  arr <- detect_array(DNAm)
  beta <- check_beta_range(DNAm)

  clock_ids <- resolve_clocks(clocks)
  seq_ids <- resolve_clocks_sequence(clock_ids)
  output_ids <- drop_routed_members(c(clock_ids, setdiff(seq_ids, clock_ids)))

  pk <- report_load_packs(seq_ids, assets, ask)
  skip_ids <- external_ids_in_groups(seq_ids, pk$skipped)
  assess_ids <- setdiff(seq_ids, skip_ids)

  panels <- clock_panels(assess_ids, pk$packs)
  needed_union <- panels_union(panels)
  mna <- scan_missing_cpgs(DNAm, needed_union)
  cpg_list <- resolve_cpgs(mna$usable_cols, panels)
  cov_tbl <- build_coverage_table(
    cpg_list,
    colnames(DNAm),
    mna,
    assess_ids,
    output_ids
  )

  samples <- report_samples(DNAm)

  missing <- list(
    total_needed_union = length(needed_union),
    n_absent_from_matrix = length(setdiff(needed_union, colnames(DNAm))),
    n_partial_probes = length(mna$partial_na_cols),
    n_full_na_probes = length(mna$all_na_cols),
    partial_probes = mna$partial_na_cols,
    full_na_probes = mna$all_na_cols
  )

  list(
    format = fmt,
    array = arr,
    beta = beta,
    missing = missing,
    samples = samples,
    coverage = cov_tbl,
    skipped_external = pk$skipped
  )
}
