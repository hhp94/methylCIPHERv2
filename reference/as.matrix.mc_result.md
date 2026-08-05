# Matrix Method For An mc_result Object

Converts the
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
output to a matrix containing just the clocks.

## Usage

``` r
# S3 method for class 'mc_result'
as.matrix(x, ...)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- ...:

  Not used.

## Value

A numeric matrix. The scores, with samples in the rows and clocks in the
columns.

## Clocks that use all the samples

Some clocks depend on information from all the samples, such as a
z-score. When `x` holds more than one batch, these clocks take their
value from every sample in `x`, and not from one batch alone. This is
the same calculation as
[`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md).

## See also

- [`as.data.frame.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/as.data.frame.mc_result.md)
  for the scores as a data.frame.

- [`rbind.mc_result()`](https://hhp94.github.io/methylCIPHERv2/reference/rbind.mc_result.md)
  for two runs combined into one object.

- [`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md)
  for a cross-sample score recomputed after a bind.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 10)
res <- calc_clocks(sim[["DNAm"]], clocks)
as.matrix(res)
#>            Horvath1     Hannum
#> sample1  147.710772  61.949225
#> sample2    2.314275 115.146393
#> sample3   28.122199  42.505738
#> sample4   71.980463  75.028900
#> sample5    2.486043  46.237341
#> sample6   97.074834  30.736410
#> sample7  105.082994  19.254252
#> sample8   50.147048  87.501966
#> sample9  119.151136  66.955452
#> sample10 105.897593   2.489874
```
