
<!-- README.md is generated from README.Rmd. Please edit that file -->

# methylCIPHERv2

<!-- badges: start -->

<!-- badges: end -->

`{methylCIPHERv2}` is a lightweight epigenetic clocks R package. The
main function, `calc_clocks()`, has 3 important inputs (see
`calc_clocks()` docs for more details):

- `DNAm`: A DNAm beta-matrix where samples are in the rows and CpGs are
  in the columns
- `clocks`: A list of clocks (e.g., Horvath1), clock "group" (e.g.,
  GrimAge), or tags ("mitotic").
- `pheno`: A data.frame that contains Age and Female (1/0), and an ID
  column that links to the sample names of the DNAm beta matrix.

## Installation

You can install the development version of `{methylCIPHERv2}` with:

``` r
pak::pkg_install("hhp94/methylCIPHERv2")
```

## Workflow

``` r
library(methylCIPHERv2)
```

### Browse the clocks

`list_clocks()` returns a data.frame showing all the supported clocks.

``` r
head(list_clocks(), n = 3)
#>    clock_id     group_id request_as covariates external        tags
#> 1 Bohlin251       Bohlin  Bohlin251               FALSE gestational
#> 2  Bohlin96       Bohlin   Bohlin96               FALSE gestational
#> 3  AdaptAge CausalityAge   AdaptAge               FALSE
```

A clock can be selected by `clock_id`, `group_id`, or `tag`.
`list_clock_tags()` returns every tag a `clocks` argument accepts.
`clock_cpgs()` returns the full CpG panel behind any set of clocks.

``` r
grim_age_cpgs <- clock_cpgs("GrimAge")
head(grim_age_cpgs, n = 3)
#> [1] "cg13947317" "cg21272576" "cg03522107"
length(grim_age_cpgs)
#> [1] 1030
```

### Download the external weights

A few clocks require downloadable assets. `list_mc_assets()` reports the
required download status.

``` r
list_mc_assets()
#>     group_id n_clocks n_cpgs   size downloaded superseded superseded_size
#> 1 PCBrainAge        1 357852  6.68M      FALSE          0               0
#> 2   PCClocks       14  78464  8.88M      FALSE          0               0
#> 3 SystemsAge       13 125175 22.46M      FALSE          0               0
#> 4  Zhang2019        1 319607  5.02M      FALSE          0               0
```

For example, `"SystemsAge"` needs its asset for `calc_clocks()`.
`download_mc_assets()` downloads these assets onto your computer.

``` r
# pass "SystemsAge" to download that group alone
download_mc_assets(groups = "all", ask = TRUE) # downloading requires confirmation
```

``` r
list_mc_assets()
#>     group_id n_clocks n_cpgs   size downloaded superseded superseded_size
#> 1 PCBrainAge        1 357852  6.68M       TRUE          0               0
#> 2   PCClocks       14  78464  8.88M       TRUE          0               0
#> 3 SystemsAge       13 125175 22.46M       TRUE          0               0
#> 4  Zhang2019        1 319607  5.02M       TRUE          0               0
```

The function `clear_mc_assets()` deletes the downloaded assets.
`vignette("assets")` covers the assets directory and the advanced path
settings.

### Calculate Epigenetic Clocks

We simulate some data that carries the two kinds of missing data: CpGs
that are absent altogether, and CpGs that are present but missing in
some samples.

``` r
set.seed(1)
sim <- sim_DNAm(
  c("Horvath1", "Bohlin96"),
  n = 6,
  Age = TRUE,
  remove = 10 # drop CpGs entirely, to simulate completely missing CpGs
)

# simulate partially missing CpGs
sim[["DNAm"]][sample.int(length(sim[["DNAm"]]), 20)] <- NA

# count partial NAs
sum(is.na(sim[["DNAm"]]))
#> [1] 20
```

The two kinds of missingness are handled differently. A CpG that is
absent altogether is either filled from a reference value that ships
with the clock, or dropped. A CpG that is present but missing in some
samples is column mean imputed (i.e., CpG-wise).

Some clocks need a `Female` covariate. If not readily available,
`predict_sex()` scores the two `DNAmSex_Wang` ChrX and ChrY clocks and
returns a predicted sex for every sample. When `pheno` provides
`Female`, `predict_sex()` also returns `recorded_sex` and
`sex_mismatch`, and reports how many samples disagree.

``` r
# The betas are random here, so the mismatches are expected.
sex_sim <- sim_DNAm(
  c("DNAmSex_Wang_ChrX", "DNAmSex_Wang_ChrY"),
  n = 6,
  Female = TRUE
)

predict_sex(sex_sim[["DNAm"]], sex_sim[["pheno"]])
#>        ID DNAmSex_Wang_ChrX DNAmSex_Wang_ChrY predicted_sex recorded_sex
#> 1 sample1          37.76046         -4.768303        Female       Female
#> 2 sample2          35.86986         -4.877966        Female       Female
```

Clock calculation is performed by the `calc_clocks()` function.

``` r
res <- calc_clocks(
  DNAm = sim[["DNAm"]],
  clocks = c("Horvath1", "Bohlin96"),
  pheno = sim[["pheno"]]
)
```

`calc_clocks()` can also run in batches and be combined afterwards.

``` r
batch_1 <- sim_DNAm(c("Knight"), n = 6, suffix = "1")
batch_2 <- sim_DNAm(c("Knight"), n = 6, suffix = "2")

r1 <- calc_clocks(batch_1[["DNAm"]], "Knight", pheno = batch_1[["pheno"]])
r2 <- calc_clocks(batch_2[["DNAm"]], "Knight", pheno = batch_2[["pheno"]])

rbind(r1, r2)
```

`as.data.frame()` or `as.matrix()` finalizes the scores for further
analysis. `as.data.frame()` returns long or wide.

``` r
# long, the default
head(as.data.frame(res), n = 3)
#>        ID clock_id       score
#> 1 sample1 Horvath1   0.1544864
#> 2 sample2 Horvath1  72.9475426
#> 3 sample3 Horvath1 107.7384507
# wide
head(as.data.frame(res, long = FALSE), n = 3)
#>        ID    Horvath1 Bohlin96
#> 1 sample1   0.1544864 187.9799
#> 2 sample2  72.9475426 239.7353
#> 3 sample3 107.7384507 227.8278
```

### Check the coverage

Many factors affect the validity of a calculated clock. One of the most
important is missing data, both for each clock and for each sample.

`clocks_coverage()` reports missing data by `clock_id`. Its last column,
`missing_cpgs`, names every CpG that was missing (dropped here for
width).

``` r
cov <- clocks_coverage(res, all_columns = TRUE)
h1 <- cov[cov[["clock_id"]] == "Horvath1", ]
cov[setdiff(names(cov), grep("missing_cpgs|norm|role", names(cov), value = TRUE))]
#>   clock_id group_id policy score_needed score_present score_used
#> 1 Horvath1 Horvath1   omit          353           344        344
#> 2 Bohlin96   Bohlin   omit           96            95         95
#>   score_imputed_partial score_imputed_full score_dropped
#> 1                    17                  0             9
#> 2                     3                  0             1
```

`Horvath1` uses the `omit` policy, so the 9 absent CpGs were dropped
rather than filled, and the 20 missing values landed on 17 distinct CpGs
that were imputed.

`samples_coverage()` covers missingness per sample per clock.

``` r
head(samples_coverage(res), n = 3)
#>        id clock_id panel n_observed n_needed  coverage
#> 1 sample1 Horvath1 score        343      353 0.9716714
#> 2 sample2 Horvath1 score        341      353 0.9660057
#> 3 sample3 Horvath1 score        342      353 0.9688385
```

### Age acceleration

`calc_accel()` regresses each score on chronological age and returns the
residual (`type = "accel"`) or calculates the raw difference between a
clock and chronological age (`type = "diff"`). The value returned by
`calc_clocks()` keeps only the id column and the covariates required for
the calculation, so chronological age is supplied through `data` here.

``` r
# long, the default
head(calc_accel(res, data = sim[["pheno"]]), n = 3)
#>        ID clock_id  accel_id      accel
#> 1 sample1 Horvath1 Age_accel -75.472920
#> 2 sample2 Horvath1 Age_accel   9.107829
#> 3 sample3 Horvath1 Age_accel  26.007752
# can also return wide
head(calc_accel(res, data = sim[["pheno"]], long = FALSE), n = 3)
#>        ID Horvath1_Age_accel Bohlin96_Age_accel
#> 1 sample1         -75.472920         -37.499651
#> 2 sample2           9.107829           8.015624
#> 3 sample3          26.007752           5.579168
```

### Check the age associations

`score_associations()` compares each clock's observed correlation with
age against a reference range that ships with the package. Age is
supplied through `age` here, for the same reason it goes through `data`
above.

``` r
# Random data for illustration purposes only
score_associations(res, age = sim[["pheno"]][["Age"]])
#>   clock_id n obs_age_r exp_age_r exp_lo exp_hi outside wrong_sign
#> 1 Horvath1 6    -0.124     0.827  0.226  0.972    TRUE       TRUE
```

### Cite the clocks

`cite_clocks()` returns the reference for every clock in a result.

``` r
cite_clocks(res)
```

Or pass the desired `clock_id` directly.

``` r
cite_clocks("Horvath1")
```
