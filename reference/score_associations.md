# Clock Age Associations

Compares each clock's observed correlation with age against a reference
range shipped with the package.

## Usage

``` r
score_associations(x, age = NULL)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- age:

  A numeric vector. The age for each sample, in the same order as the
  samples in `x`. Default is `NULL`, which uses the `Age` column of the
  `pheno` in `x`.

## Value

A data.frame. One row for each clock with a reference entry and enough
samples. Columns are the clock id, the sample count, the observed and
reference age correlations, the reference range, and the two flags
described above.

## Details

Each row compares a clock's observed score-age correlation with its
reference correlation and expected range. The `outside` column marks a
clock whose observed correlation falls outside that range. The
`wrong_sign` column marks a clock whose observed correlation has the
opposite sign from a reference correlation stronger than 0.3.

A clock needs at least 5 samples with a finite score and age, and
variation in both, to appear in the result. A clock with no entry in the
shipped reference table is left out.

Rows are ordered by the gap between the observed and the reference
correlation, most negative first.

## Clocks that use all the samples

Some clocks depend on information from all the samples, such as a
z-score. When `x` holds more than one batch, these clocks take their
value from every sample in `x`, and not from one batch alone. This is
the same calculation as
[`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md).

## See also

[`calc_accel()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_accel.md)
for the age acceleration of each sample.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20)
res <- calc_clocks(sim[["DNAm"]], clocks)
score_associations(res, age = runif(20, 20, 80))
#>   clock_id  n obs_age_r exp_age_r exp_lo exp_hi outside wrong_sign
#> 1 Horvath1 20     0.006     0.827  0.226  0.972    TRUE      FALSE
#> 2   Hannum 20     0.038     0.845  0.261  0.976    TRUE      FALSE
```
