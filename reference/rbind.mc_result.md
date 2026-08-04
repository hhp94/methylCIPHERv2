# Combined Batches Of Scores

Stacks two or more outputs from
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
runs into one object of multiple batches.

## Usage

``` r
# S3 method for class 'mc_result'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  Two or more `mc_result` objects.

- deparse.level:

  A single whole number. Not used by this method. Default is `1`.

## Value

An `mc_result` object. It holds the stacked scores, `pheno`, and
`coverage` of every input, under one `mc_batch_id` label for each.

## Details

Each input must use disjoint sample ids, the same scored clocks, the
same `pheno_id`, and the same `normalize` setting.
[`rbind()`](https://rdrr.io/r/base/cbind.html) stops when any of those
differ between inputs.

The combined value gets one `mc_batch_id` label for each input. A clock
that depends on sample-wise information, such as a z-score, keeps the
value each input calculated on its own samples. Call
[`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md)
to calculate it again from every sample in the combined value.

## See also

- [`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md)
  for a cross-sample score recomputed after a bind.

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
combined
#> <mc_result> 20 sample(s) x 2 clock(s)
#> 
#> $scores [6 of 20 row(s), 2 of 2 clock(s)]
#>          Horvath1   Hannum
#> sample1  59.68812 92.04382
#> sample2  49.21227 36.27282
#> sample3 155.59404 60.09078
#> sample4 174.61964 52.25514
#> sample5 104.41792 15.09647
#> sample6  70.62510 68.06703
#> ... 14 more row(s)
#> 
#> $pheno [6 of 20 row(s), 1 column(s)]
#>        ID
#> 1 sample1
#> 2 sample2
#> 3 sample3
#> 4 sample4
#> 5 sample5
#> 6 sample6
#> ... 14 more row(s)
#> 
#> $provenance [2 batch(es)]
#> a683cd7ddcd49f1f, c877d86e3e1f0aa8
```
