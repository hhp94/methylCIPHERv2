# Print Method For An mc_result Object

Prints the scores and the `pheno` table for an `mc_result` object.

## Usage

``` r
# S3 method for class 'mc_result'
print(x, n = 6, p = 6, ...)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- n:

  A single whole number. The number of sample rows to print. Default is
  `6`.

- p:

  A single whole number. The number of clock columns to print for the
  scores table. Default is `6`.

- ...:

  Not used.

## Value

An `mc_result` object. Returns `x`, invisibly, after printing it.

## Details

The output lists batch labels only when `x` spans more than one
`mc_batch_id`.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20)
res <- calc_clocks(sim[["DNAm"]], clocks)
print(res, n = 3, p = 2)
#> <mc_result> 20 sample(s) x 2 clock(s)
#> 
#> $scores [3 of 20 row(s), 2 of 2 clock(s)]
#>          Horvath1   Hannum
#> sample1 133.34041 74.86782
#> sample2 128.90027 93.86625
#> sample3  85.39902 36.82550
#> ... 17 more row(s)
#> 
#> $pheno [3 of 20 row(s), 1 column(s)]
#>        ID
#> 1 sample1
#> 2 sample2
#> 3 sample3
#> ... 17 more row(s)
```
