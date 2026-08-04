# Clock Coverage Counts

Reports the CpG counts behind each clock's score in `x`, one row for
each clock and batch.

## Usage

``` r
clocks_coverage(x, all_columns = FALSE)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- all_columns:

  A boolean. Returns every column, including the ones the frame leaves
  out by default. Default is `FALSE`.

## Value

A data.frame. One row for each clock and batch, with the CpG counts of
its scoring panel, and the columns above that apply to `x`.

## Details

Every row gives the clock's `group_id`, its `policy`, and the six
`score_*` counts for the CpGs in its own scoring panel.

A clock assembled only from other clocks' scores reads no CpGs of its
own, and gets a row of `NA` counts. Read the coverage of the clocks it
depends on instead.

Four more columns appear only where they say something about `x`.

- `role` appears when `x` holds a clock that scores as part of another
  clock.

- `normalizes`, and the five `norm_*` counts for the normalizing panel,
  appear when a clock in `x` normalizes. There is no `norm_used`.

- `missing_cpgs` lists the CpGs absent from a clock's panel, and appears
  when a CpG is absent.

- `mc_batch_id` appears when `x` holds more than one batch. A batch is
  the set of samples from one
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
  call, and [`rbind()`](https://rdrr.io/r/base/cbind.html) combines
  batches.

Pass `all_columns = TRUE` to keep `role`, `normalizes`, the `norm_*`
counts, and `missing_cpgs` in every frame. Use it where the code that
reads the frame names a column directly. `mc_batch_id` is the one
exception, and still appears only when `x` holds more than one batch.

## See also

[`samples_coverage()`](https://hhp94.github.io/methylCIPHERv2/reference/samples_coverage.md)
for the same panels counted for each sample.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20)
res <- calc_clocks(sim[["DNAm"]], clocks)

clocks_coverage(res)
#>   clock_id group_id policy score_needed score_present score_used
#> 1 Horvath1 Horvath1   omit          353           353        353
#> 2   Hannum   Hannum   omit           71            71         71
#>   score_imputed_partial score_imputed_full score_dropped
#> 1                     0                  0             0
#> 2                     0                  0             0
clocks_coverage(res, all_columns = TRUE)
#>   clock_id group_id     role policy normalizes score_needed score_present
#> 1 Horvath1 Horvath1 returned   omit      FALSE          353           353
#> 2   Hannum   Hannum returned   omit      FALSE           71            71
#>   score_used score_imputed_partial score_imputed_full score_dropped norm_needed
#> 1        353                     0                  0             0           0
#> 2         71                     0                  0             0           0
#>   norm_present norm_imputed_partial norm_imputed_full norm_dropped missing_cpgs
#> 1            0                    0                 0            0             
#> 2            0                    0                 0            0             
```
