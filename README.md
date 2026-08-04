
<!-- README.md is generated from README.Rmd. Please edit that file -->

# methylCIPHERv2

<!-- badges: start -->

<!-- badges: end -->

`{methylCIPHERv2}` scores CpG-based DNA methylation ages, commonly
called epigenetic clocks. It reads a beta matrix and returns one score
per sample per clock, together with the coverage counts that say how
much of each clock's CpG panel was present in the data.

`calc_clocks()` is the single entry point. Every other function either
describes the available clocks before the calculation, or inspects the
results after it.

## Installation

`{methylCIPHERv2}` is not on CRAN yet.

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
head(list_clocks())
#>    clock_id     group_id request_as covariates external        tags
#> 1 Bohlin251       Bohlin  Bohlin251               FALSE gestational
#> 2  Bohlin96       Bohlin   Bohlin96               FALSE gestational
#> 3  AdaptAge CausalityAge   AdaptAge               FALSE            
#> 4   CausAge CausalityAge    CausAge               FALSE            
#> 5    DamAge CausalityAge     DamAge               FALSE            
#> 6 CellDRIFT    CellDRIFT  CellDRIFT               FALSE
```

A clock can be selected by `clock_id`, `group_id`, or `tag`.
`list_clock_tags()` returns every tag a `clocks` argument accepts.

``` r
list_clock_tags()
#> $gestational
#> [1] "Bohlin"          "Knight"          "Mayne"           "LeePlacentalAge"
#> 
#> $mitotic
#> [1] "EpiTOC"    "EpiTOC2"   "MiAge"     "RepliTali"
#> 
#> $mortality
#> [1] "GrimAge"        "ZhangMortality"
```

`clock_cpgs()` returns the full CpG panel behind any set of clocks.

``` r
grim_age_cpgs <- clock_cpgs("GrimAge")
head(grim_age_cpgs)
#> [1] "cg13947317" "cg21272576" "cg03522107" "cg03259703" "cg05461666"
#> [6] "cg15860924"
length(grim_age_cpgs)
#> [1] 1030
```

### Download the external weights

Most clocks ship inside the package. A few require downloadable assets.
`list_mc_assets()` reports each group, its size, and whether it is
already on disk.

``` r
list_mc_assets()
#>     group_id n_clocks n_cpgs   size downloaded superseded superseded_size
#> 1 PCBrainAge        1 357852  6.68M      FALSE          0               0
#> 2   PCClocks       14  78464  8.88M      FALSE          0               0
#> 3 SystemsAge       13 125175 22.46M      FALSE          0               0
#> 4  Zhang2019        1 319607  5.02M      FALSE          0               0
```

For example, `"SystemsAge"` needs its asset. `download_mc_assets()`
downloads these assets.

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

The function `clear_mc_assets()` deletes the downloaded assets, also
with confirmation. `vignette("assets")` covers the assets directory and
the advanced path settings.

### Calculate Epigenetic Clocks

Start from simulated data, carrying the two kinds of missing data: CpGs
that are absent altogether, and CpGs that are present but missing in
some samples.

``` r
set.seed(1)
sim <- sim_DNAm(
  "Horvath1",
  n = 6,
  Age = TRUE,
  remove = 10 # drop CpGs entirely, to simulate completely missing CpGs
)

# simulate partially missing CpGs
sim[["DNAm"]][sample.int(length(sim[["DNAm"]]), 20)] <- NA

sim
#> <mc_sim> 6 sample(s) x 343 CpG(s)
#> 
#> $DNAm [6 of 6 row(s), 6 of 343 CpG(s)]
#>         cg00075967 cg00374717 cg00864867 cg00945507 cg01027739 cg01353448
#> sample1  0.7774452 0.38611409  0.5995658  0.1079436  0.5530363 0.69273156
#> sample2  0.9347052 0.01339033  0.4935413  0.7237109  0.5297196 0.47761962
#> sample3  0.2121425 0.38238796  0.1862176  0.4112744  0.7893562 0.86120948
#> sample4  0.6516738 0.86969085  0.8273733  0.8209463  0.0233312 0.43809711
#> sample5  0.1255551 0.34034900  0.6684667  0.6470602  0.4772301 0.24479728
#> sample6  0.2672207 0.48208012  0.7942399  0.7829328  0.7323137 0.07067905
#> ... 337 more CpG(s)
#> 
#> $pheno [6 of 6 row(s), 2 column(s)]
#>        ID      Age
#> 1 sample1 44.98672
#> 2 sample2 40.81368
#> 3 sample3 45.59427
#> 4 sample4 42.49496
#> 5 sample5 50.12672
#> 6 sample6 46.08904

sum(is.na(sim[["DNAm"]]))
#> [1] 20
```

The two kinds of missingness are handled differently. A CpG that is
absent altogether is either filled from a reference value that ships
with the clock, or dropped. The `policy` column of `clocks_coverage()`
says which one a clock uses. A CpG that is present but missing in some
samples is imputed from the samples that do have it.

Some clocks need a `Female` covariate, which the `covariates` column of
`list_clocks()` marks. `predict_sex()` scores the two `DNAmSex_Wang`
clocks and returns a predicted sex for every sample. When `pheno` does
carry `Female`, `predict_sex()` also returns `recorded_sex` and
`sex_mismatch`, and reports how many samples disagree.

``` r
sex_sim <- sim_DNAm(
  c("DNAmSex_Wang_ChrX", "DNAmSex_Wang_ChrY"),
  n = 6,
  Female = TRUE
)

predict_sex(sex_sim[["DNAm"]], sex_sim[["pheno"]])
#> ! 3 samples have a predicted sex that does not match the Female column in
#>   `pheno`.
#> ℹ The sex_mismatch column marks those samples.
#> ℹ A mismatch can come from the recorded sex or from the array data. Check both
#>   sources before you correct either one.
#>        ID DNAmSex_Wang_ChrX DNAmSex_Wang_ChrY predicted_sex recorded_sex
#> 1 sample1          37.76046         -4.768303        Female       Female
#> 2 sample2          35.86986         -4.877966        Female       Female
#> 3 sample3          38.98598         -5.952733        Female         Male
#> 4 sample4          36.64618         -4.247936        Female         Male
#> 5 sample5          37.02444         -5.314503        Female       Female
#> 6 sample6          37.35877         -6.296180        Female         Male
#>   sex_mismatch
#> 1        FALSE
#> 2        FALSE
#> 3         TRUE
#> 4         TRUE
#> 5        FALSE
#> 6         TRUE
```

The betas are random here, so every predicted sex is arbitrary and the
mismatches are expected.

Clock calculation is performed by the `calc_clocks()` function.

``` r
res <- calc_clocks(
  sim[["DNAm"]],
  "Horvath1",
  pheno = sim[["pheno"]]
)
res
#> <mc_result> 6 sample(s) x 1 clock(s)
#> 
#> $scores [6 of 6 row(s), 1 of 1 clock(s)]
#>         Horvath1
#> sample1 67.52456
#> sample2 30.69383
#> sample3  2.85157
#> sample4 27.32398
#> sample5 76.43164
#> sample6 77.37217
#> 
#> $pheno [6 of 6 row(s), 1 column(s)]
#>        ID
#> 1 sample1
#> 2 sample2
#> 3 sample3
#> 4 sample4
#> 5 sample5
#> 6 sample6
```

`calc_clocks()` can also run in batches, which are combined afterwards.

``` r
batch_1 <- sim_DNAm(c("Knight"), n = 6, suffix = "1")
batch_2 <- sim_DNAm(c("Knight"), n = 6, suffix = "2")

r1 <- calc_clocks(batch_1[["DNAm"]], "Knight", pheno = batch_1[["pheno"]])
r2 <- calc_clocks(batch_2[["DNAm"]], "Knight", pheno = batch_2[["pheno"]])

rbind(r1, r2)
#> <mc_result> 12 sample(s) x 1 clock(s)
#> 
#> $scores [6 of 12 row(s), 1 of 1 clock(s)]
#>               Knight
#> sample1_1 -12.297759
#> sample2_1  -6.127835
#> sample3_1   2.200692
#> sample4_1  -6.358762
#> sample5_1  -1.096765
#> sample6_1   5.706422
#> ... 6 more row(s)
#> 
#> $pheno [6 of 12 row(s), 1 column(s)]
#>          ID
#> 1 sample1_1
#> 2 sample2_1
#> 3 sample3_1
#> 4 sample4_1
#> 5 sample5_1
#> 6 sample6_1
#> ... 6 more row(s)
#> 
#> $provenance [2 batch(es)]
#> 349572599947870a, 38f26e39f4c6a41e
```

`as.data.frame()` or `as.matrix()` finalize the scores for further
analysis. `as.data.frame()` returns long or wide.

``` r
# long, the default
head(as.data.frame(res))
#>        ID clock_id    score
#> 1 sample1 Horvath1 67.52456
#> 2 sample2 Horvath1 30.69383
#> 3 sample3 Horvath1  2.85157
#> 4 sample4 Horvath1 27.32398
#> 5 sample5 Horvath1 76.43164
#> 6 sample6 Horvath1 77.37217
# wide
head(as.data.frame(res, long = FALSE))
#>        ID Horvath1
#> 1 sample1 67.52456
#> 2 sample2 30.69383
#> 3 sample3  2.85157
#> 4 sample4 27.32398
#> 5 sample5 76.43164
#> 6 sample6 77.37217
```

### Check the coverage

Many factors affect the validity of a calculated clock. The most
important is missing data, both for each clock and for each sample.

`clocks_coverage()` reports missing data by `clock_id`. Its last column,
`missing_cpgs`, names every CpG that was missing (dropped here for
width).

``` r
cov <- clocks_coverage(res, all_columns = TRUE)
cov[setdiff(names(cov), grep("missing_cpgs|norm|role", names(cov), value = TRUE))]
#>   clock_id group_id policy score_needed score_present score_used
#> 1 Horvath1 Horvath1   omit          353           343        343
#>   score_imputed_partial score_imputed_full score_dropped
#> 1                    19                  0            10
```

`Horvath1` uses the `omit` policy, so the 10 absent CpGs were dropped
rather than filled, and the 20 missing values landed on 19 distinct CpGs
that were imputed.

`samples_coverage()` covers missingness per sample per clock.

``` r
head(samples_coverage(res))
#>        id clock_id panel n_observed n_needed  coverage
#> 1 sample1 Horvath1 score        340      353 0.9631728
#> 2 sample2 Horvath1 score        339      353 0.9603399
#> 3 sample3 Horvath1 score        342      353 0.9688385
#> 4 sample4 Horvath1 score        338      353 0.9575071
#> 5 sample5 Horvath1 score        340      353 0.9631728
#> 6 sample6 Horvath1 score        339      353 0.9603399
```

### Age acceleration

`calc_accel()` regresses each score on chronological age and returns the
residual. The value returned by `calc_clocks()` keeps only the id column
and the covariates the calculation required, so chronological age is
supplied through `data`, joined on the sample id.

``` r
head(calc_accel(res, data = sim[["pheno"]]))
#>        ID clock_id  accel_id      accel
#> 1 sample1 Horvath1 Age_accel  20.655930
#> 2 sample2 Horvath1 Age_accel   6.056435
#> 3 sample3 Horvath1 Age_accel -47.253674
#> 4 sample4 Horvath1 Age_accel  -6.270204
#> 5 sample5 Horvath1 Age_accel   2.180420
#> 6 sample6 Horvath1 Age_accel  24.631093
```

### Check the age associations

`calc_accel()` works sample by sample. `score_associations()` looks at
all the samples at once, and compares each clock's observed correlation
with age against a reference range that ships with the package. Age is
supplied through `age` here, for the same reason it goes through `data`
above.

``` r
score_associations(res, age = sim[["pheno"]][["Age"]])
#>   clock_id n obs_age_r exp_age_r exp_lo exp_hi outside wrong_sign
#> 1 Horvath1 6     0.552     0.827  0.226  0.972   FALSE      FALSE
```

The scores come from random data, so the correlation itself carries no
meaning, and six samples make the reference range wide. What matters in
real use are the two flags. `outside` marks a clock whose observed
correlation falls outside the reference range, and `wrong_sign` marks
one that runs in the opposite direction from a reference correlation
stronger than 0.3.

### Cite the clocks

`cite_clocks()` returns the reference for every clock in a result, so a
methods section can be assembled from the run itself.

``` r
cite_clocks(res)
#> <mc_citation> 1 clock(s) x 1 paper(s)
#> $bibtex [1 paper(s)]
#> @article{Horvath_2013_24138928,
#>   title = {{DNA} methylation age of human tissues and cell types},
#>   author = {Horvath, Steve},
#>   year = {2013},
#>   journal = {Genome biology},
#>   volume = {14},
#>   number = {10},
#>   pages = {R115},
#>   doi = {10.1186/gb-2013-14-10-r115},
#>   pmid = {24138928},
#>   url = {https://doi.org/10.1186/gb-2013-14-10-r115}
#> }
#> ℹ `as.data.frame(x)` returns the clock-to-paper table.
#> ℹ `writeLines(toBibtex(x), "refs.bib")` writes the bibtex to a file.
#> ℹ `citation("methylCIPHERv2")` cites the package itself.
```

Or pass the desired `clock_id` directly.

``` r
cite_clocks("Horvath1")
```
