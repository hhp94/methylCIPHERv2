# Print Method For An mc_sim Object

Prints a compact summary of an `mc_sim` object, with a preview of `DNAm`
and `pheno`.

## Usage

``` r
# S3 method for class 'mc_sim'
print(x, n = 6, p = 6, ...)
```

## Arguments

- x:

  An `mc_sim` object. The value returned by
  [`sim_DNAm()`](https://hhp94.github.io/methylCIPHERv2/reference/sim_DNAm.md).

- n:

  A single whole number. The number of sample rows to preview from
  `DNAm` and `pheno`. Default is `6`.

- p:

  A single whole number. The number of CpG columns to preview from
  `DNAm`. Default is `6`.

- ...:

  Not used.

## Value

An `mc_sim` object. Returns `x`, invisibly, after printing it.

## Examples

``` r
sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 10)
print(sim)
#> <mc_sim> 10 sample(s) x 418 CpG(s)
#> 
#> $DNAm [6 of 10 row(s), 6 of 418 CpG(s)]
#>         cg00075967 cg00374717 cg00864867 cg00945507 cg01027739 cg01353448
#> sample1  0.5091074  0.1445521  0.2979470 0.05896919  0.8873238  0.9567291
#> sample2  0.3261320  0.5621892  0.7734151 0.02944346  0.6389631  0.8700143
#> sample3  0.4894682  0.0195675  0.3327870 0.34282234  0.9617274  0.5808928
#> sample4  0.3227571  0.3651035  0.8183872 0.32656887  0.5875493  0.1702297
#> sample5  0.6099244  0.3754307  0.3533365 0.41680224  0.8958470  0.5705086
#> sample6  0.1581487  0.5671309  0.5591794 0.28953103  0.7384070  0.1078505
#> ... 4 more row(s), 412 more CpG(s)
#> 
#> $pheno [6 of 10 row(s), 1 column(s)]
#>        ID
#> 1 sample1
#> 2 sample2
#> 3 sample3
#> 4 sample4
#> 5 sample5
#> 6 sample6
#> ... 4 more row(s)
```
