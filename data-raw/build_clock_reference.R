# build per-clock age-association reference (maintainer-run).
# input: TranslAGE AgeSexAssociations table. output: data-raw/clock_reference.csv.

suppressMessages(devtools::load_all(quiet = TRUE))

# local maintainer input (identifiable, per-dataset) -- NOT committed
SRC <- Sys.getenv(
  "MC_AGESEX_CSV",
  "/Users/dsb52/Documents/GitHub/TranslAGE-workflows/WORKFLOWS/1_Workflows/00_AgeSexAssociations/results/age_sex_regression.csv"
)
# shipped with the package (installed, read at runtime via system.file)
OUT <- "inst/extdata/clock_reference.csv"
MIN_N <- 20L # drop tiny datasets: their correlations are noise

# --- meta-analysis helpers -------------------------------------------------

# random-effects pool (DerSimonian-Laird). returns pooled estimate and 95% prediction interval.
re_meta <- function(y, se) {
  ok <- is.finite(y) & is.finite(se) & se > 0
  y <- y[ok]
  se <- se[ok]
  k <- length(y)
  if (k < 1L) {
    return(c(pooled = NA, lo = NA, hi = NA, tau = NA, k = 0))
  }
  if (k < 3L) {
    return(c(pooled = mean(y), lo = NA, hi = NA, tau = NA, k = k))
  }
  v <- se^2
  w <- 1 / v
  mu_fixed <- sum(w * y) / sum(w)
  Q <- sum(w * (y - mu_fixed)^2)
  c_term <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (k - 1)) / c_term)
  wr <- 1 / (v + tau2)
  mu <- sum(wr * y) / sum(wr)
  se_mu <- sqrt(1 / sum(wr))
  tcrit <- stats::qt(0.975, df = k - 2)
  half <- tcrit * sqrt(tau2 + se_mu^2)
  c(pooled = mu, lo = mu - half, hi = mu + half, tau = sqrt(tau2), k = k)
}

# correlation meta-analysis on the Fisher-z scale, back-transformed
re_meta_cor <- function(r, n) {
  ok <- is.finite(r) & is.finite(n) & n > 3 & abs(r) < 1
  m <- re_meta(atanh(r[ok]), 1 / sqrt(n[ok] - 3))
  c(
    pooled = tanh(m[["pooled"]]),
    lo = tanh(m[["lo"]]),
    hi = tanh(m[["hi"]]),
    k = m[["k"]]
  )
}

rnd <- function(x, d) unname(round(as.numeric(x), d))

# --- load, filter to catalog clocks ---------------------------------------

d <- utils::read.csv(SRC, stringsAsFactors = FALSE)
catalog <- resolve_clocks("all")
d <- d[d$clock %in% catalog & is.finite(d$n) & d$n >= MIN_N, , drop = FALSE]

matched <- sort(unique(d$clock))
cat(sprintf(
  "matched %d of %d catalog clocks; %d have no reference data\n",
  length(matched),
  length(catalog),
  length(setdiff(catalog, matched))
))

# --- per-clock summary ------------------------------------------------------

rows <- lapply(matched, function(cl) {
  s <- d[d$clock == cl, , drop = FALSE]
  ar <- re_meta_cor(s$pearson_all, s$n)
  asl <- re_meta(s$age_coef, s$age_se)
  sx <- re_meta(s$sex_coef, s$sex_se)
  ax <- re_meta(s$ia_age_x_sex_coef, s$ia_age_x_sex_se)
  data.frame(
    clock = cl,
    k_datasets = as.integer(ar[["k"]]),
    total_n = sum(s$n, na.rm = TRUE),
    # expected age correlation and its plausible (prediction) interval
    age_r = rnd(ar[["pooled"]], 3),
    age_r_lo = rnd(ar[["lo"]], 3),
    age_r_hi = rnd(ar[["hi"]], 3),
    # expected age slope (units of the clock per year) -- less range-dependent
    age_slope = rnd(asl[["pooled"]], 4),
    age_slope_lo = rnd(asl[["lo"]], 4),
    age_slope_hi = rnd(asl[["hi"]], 4),
    # expected sex effect (clock units, male vs female) + how often it matters
    sex_coef = rnd(sx[["pooled"]], 4),
    sex_coef_lo = rnd(sx[["lo"]], 4),
    sex_coef_hi = rnd(sx[["hi"]], 4),
    sex_frac_sig = rnd(mean(s$sex_p < 0.05, na.rm = TRUE), 2),
    # expected age x sex interaction (sex-specific rate of aging) + frequency
    age_x_sex = rnd(ax[["pooled"]], 5),
    age_x_sex_lo = rnd(ax[["lo"]], 5),
    age_x_sex_hi = rnd(ax[["hi"]], 5),
    age_x_sex_frac_sig = rnd(mean(s$ia_age_x_sex_p < 0.05, na.rm = TRUE), 2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
})
ref <- do.call(rbind, rows)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(ref, OUT, row.names = FALSE)
cat(sprintf("wrote %s (%d clocks x %d cols)\n", OUT, nrow(ref), ncol(ref)))
