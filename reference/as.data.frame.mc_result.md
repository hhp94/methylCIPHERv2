# Data Frame Method For An mc_result Object

Converts the
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
output to a data.frame containing just the clocks, in long or wide
format.

## Usage

``` r
# S3 method for class 'mc_result'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., long = TRUE)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- row.names:

  A character vector. Not used by this method. Default is `NULL`.

- optional:

  A boolean. Not used by this method. Default is `FALSE`.

- ...:

  Not used.

- long:

  A boolean. Returns one row for each sample and clock when `TRUE`, and
  one row for each sample, with one column for each clock, when `FALSE`.
  Default is `TRUE`.

## Value

A data.frame. In long form, one row for each sample and clock, with the
score and, when `x` holds more than one batch, an `mc_batch_id` column.
In wide form, one row for each sample, with one column for each clock.

## Details

This function recalculates any clock that depends on sample-wise
information, such as a z-score, from all the available samples when `x`
holds more than one batch. This is the same calculation as
[`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md).

## See also

- [`as.matrix.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/as.matrix.mc_result.md)
  for the scores as a numeric matrix.

- [`rbind.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/rbind.mc_result.md)
  for two runs combined into one object.

- [`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md)
  for a cross-sample score recomputed after a bind.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20)
res <- calc_clocks(sim[["DNAm"]], clocks)

head(as.data.frame(res))
#>        ID clock_id      score
#> 1 sample1 Horvath1   4.354637
#> 2 sample2 Horvath1 171.777818
#> 3 sample3 Horvath1  90.340034
#> 4 sample4 Horvath1  67.993015
#> 5 sample5 Horvath1  27.835940
#> 6 sample6 Horvath1  69.351179
head(as.data.frame(res, long = FALSE))
#>        ID   Horvath1   Hannum
#> 1 sample1   4.354637 67.80994
#> 2 sample2 171.777818 51.69911
#> 3 sample3  90.340034 96.12068
#> 4 sample4  67.993015 46.29666
#> 5 sample5  27.835940 70.10949
#> 6 sample6  69.351179 53.91436
```
