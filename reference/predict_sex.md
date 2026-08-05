# Predicted Sex Karyotype

Predicts sex and identifies sex chromosome aneuploidy.

## Usage

``` r
predict_sex(DNAm, pheno = NULL, ...)
```

## Arguments

- DNAm:

  A numeric matrix. The methylation beta values, with samples in the
  rows and CpGs in the columns.

- pheno:

  A data.frame. The sample metadata, with one row for each sample.
  Default is `NULL`.

- ...:

  Passed to
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

## Value

A data.frame. One row for each sample, with the `DNAmSex_Wang_ChrX` and
`DNAmSex_Wang_ChrY` scores, `predicted_sex`, and, when `pheno` has a
`Female` column, `recorded_sex` and `sex_mismatch`.

## Details

This is a re-implementation of the sex prediction algorithm of the
wateRmelon package.

`predicted_sex` is one of `"Male"`, `"Female"`, `"47,XXY"`, or
`"45,XO"`. A sample missing either score gets `NA`, not a default call.

When `pheno` has a `Female` column, coded `0` or `1`, the result also
carries `recorded_sex` and `sex_mismatch`. `sex_mismatch` is `TRUE` only
where `predicted_sex` disagrees with a binary `recorded_sex`. A
`"47,XXY"` or `"45,XO"` call is never flagged, because a binary `Female`
column cannot record it.

## References

Wang Y, Hannon E, Grant OA, Gorrie-Stone TJ, Kumari M, Mill J, Zhai X,
McDonald-Maier KD, Schalkwyk LC (2021). DNA methylation-based sex
classifier to predict sex and identify sex chromosome aneuploidy. *BMC
Genomics*, 22(1), 484.
[doi:10.1186/s12864-021-07675-2](https://doi.org/10.1186/s12864-021-07675-2)

## Examples

``` r
sim <- sim_DNAm("DNAmSex_Wang", n = 6, Female = TRUE)
predict_sex(sim[["DNAm"]], sim[["pheno"]])
#> ! 3 samples have a predicted sex that does not match the Female column in
#>   `pheno`.
#> ℹ The sex_mismatch column marks those samples.
#> ℹ A mismatch can come from the recorded sex or from the array data. Check both
#>   sources before you correct either one.
#>        ID DNAmSex_Wang_ChrX DNAmSex_Wang_ChrY predicted_sex recorded_sex
#> 1 sample1          38.45090         -4.196821        Female       Female
#> 2 sample2          39.03248         -5.242193        Female       Female
#> 3 sample3          36.47824         -5.032105        Female         Male
#> 4 sample4          38.08351         -6.235841        Female         Male
#> 5 sample5          37.04119         -4.226353        Female       Female
#> 6 sample6          38.61471         -6.533602        Female         Male
#>   sex_mismatch
#> 1        FALSE
#> 2        FALSE
#> 3         TRUE
#> 4         TRUE
#> 5        FALSE
#> 6         TRUE
```
