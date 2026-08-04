# Age Acceleration Or Difference

Computes age acceleration or age difference for every clock in `x`.

## Usage

``` r
calc_accel(
  x,
  formula = NULL,
  type = c("accel", "diff"),
  data = NULL,
  long = TRUE
)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- formula:

  A one-sided formula. The model fit against each clock's score. Default
  is `NULL`, which uses `~ Age`.

- type:

  One of "accel" or "diff". The quantity to compute for each clock.
  Default is `"accel"`.

- data:

  A data.frame. Extra sample metadata, joined to the `pheno` in `x` by
  sample id. Default is `NULL`.

- long:

  A boolean. Returns one row for each sample and clock when `TRUE`, and
  one row for each sample, with one column for each clock, when `FALSE`.
  Default is `TRUE`.

## Value

A data.frame. In long form, one row for each sample and clock, with the
fitted value, an `accel_id` column that names the model, and, when `x`
holds more than one batch, `mc_batch_id`. In wide form, one row for each
sample, with one column for each clock.

## Details

This function recalculates any clock that depends on sample-wise
information, such as a z-score, from all the available samples when `x`
holds more than one batch. This is the same calculation as
[`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md).

The default `type = "accel"` calculates the well-known age acceleration.
It regresses each clock in `x` on `Age` and returns the residuals.
`type = "diff"` calculates the raw difference between each clock and
`Age`, and fits no model unless `formula` is given.

`formula` replaces the default model completely. It does not add to it,
so `~ Plate` regresses each clock on the plate alone, and not on age.

`data` may carry the covariates the calculation needs, as in
`data.frame(ID, Plate)` passed to `data`, with
`formula = ~ Age + Plate`. It may add a column, and it may repeat a
column that scoring already used. `calc_accel()` stops when a repeated
column disagrees with the value scoring used.

## See also

[`score_associations()`](https://hhp94.github.io/methylCIPHERv2/reference/score_associations.md)
for how each clock tracks age against a reference.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20, Age = TRUE, Female = TRUE)
res <- calc_clocks(sim[["DNAm"]], clocks)

# accel with no formula regresses each clock's score on ~ Age
head(calc_accel(res, data = sim[["pheno"]]))
#>        ID clock_id  accel_id      accel
#> 1 sample1 Horvath1 Age_accel   9.112710
#> 2 sample2 Horvath1 Age_accel  -6.811005
#> 3 sample3 Horvath1 Age_accel  -1.613302
#> 4 sample4 Horvath1 Age_accel -11.763407
#> 5 sample5 Horvath1 Age_accel -26.564985
#> 6 sample6 Horvath1 Age_accel  64.947493

# a formula replaces the default model, so name every term you want
pheno <- sim[["pheno"]]
pheno[["Plate"]] <- sample(c("P1", "P2"), nrow(pheno), replace = TRUE)
head(calc_accel(res, formula = ~ Age + Plate, data = pheno))
#>        ID clock_id        accel_id      accel
#> 1 sample1 Horvath1 Age_Plate_accel   9.630788
#> 2 sample2 Horvath1 Age_Plate_accel  -7.553565
#> 3 sample3 Horvath1 Age_Plate_accel  -1.048109
#> 4 sample4 Horvath1 Age_Plate_accel -11.169701
#> 5 sample5 Horvath1 Age_Plate_accel -27.247784
#> 6 sample6 Horvath1 Age_Plate_accel  65.466570

# diff with no formula is the raw difference from age, with no model fit
head(calc_accel(res, type = "diff", data = sim[["pheno"]]))
#>        ID clock_id accel_id      accel
#> 1 sample1 Horvath1     diff  47.341750
#> 2 sample2 Horvath1     diff  32.228082
#> 3 sample3 Horvath1     diff  33.856083
#> 4 sample4 Horvath1     diff  22.035982
#> 5 sample5 Horvath1     diff   8.973753
#> 6 sample6 Horvath1     diff 103.118082

# diff with a formula residualizes the difference
head(calc_accel(res, type = "diff", formula = ~ Age, data = sim[["pheno"]]))
#>        ID clock_id accel_id      accel
#> 1 sample1 Horvath1 Age_diff   9.112710
#> 2 sample2 Horvath1 Age_diff  -6.811005
#> 3 sample3 Horvath1 Age_diff  -1.613302
#> 4 sample4 Horvath1 Age_diff -11.763407
#> 5 sample5 Horvath1 Age_diff -26.564985
#> 6 sample6 Horvath1 Age_diff  64.947493
```
