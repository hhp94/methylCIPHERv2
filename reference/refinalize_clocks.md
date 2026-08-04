# Scores Recomputed From All Samples

Recalculates every clock that depends on sample-wise information, such
as a z-score, using the samples now in `x`.

## Usage

``` r
refinalize_clocks(x)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

## Value

An `mc_result` object. The same as `x`, with any score that is computed
from all its samples together computed again from every sample in `x`.

## Details

A clock of that kind is calculated once from every sample in a
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
call, not one sample at a time. After
[`rbind()`](https://rdrr.io/r/base/cbind.html) combines several such
values, each score still holds the value its own input calculated.
`refinalize_clocks()` calculates it again from every sample in `x`.

`refinalize_clocks()` changes nothing, and returns `x` unchanged, when
`x` holds no clock of that kind.

## See also

- [`rbind.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/rbind.mc_result.md)
  for two runs combined into one object.

- [`as.data.frame.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/as.data.frame.mc_result.md)
  for the scores as a data.frame.

- [`as.matrix.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/as.matrix.mc_result.md)
  for the scores as a numeric matrix.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim1 <- sim_DNAm(clocks, n = 10)
sim2 <- sim_DNAm(clocks, n = 10, suffix = "b")

res1 <- calc_clocks(sim1[["DNAm"]], clocks)
res2 <- calc_clocks(sim2[["DNAm"]], clocks)
combined <- rbind(res1, res2)

# a no-op here, because neither Horvath1 nor Hannum is scored from all
# samples together
refinalize_clocks(combined)
#> ℹ `refinalize_clocks()` changes only clocks that are scored from all the
#>   samples together. `x` has none.
```
