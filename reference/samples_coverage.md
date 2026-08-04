# Sample Coverage Counts

Reports each sample's CpG coverage for every clock in `x`, one row for
each sample, clock, and panel.

## Usage

``` r
samples_coverage(x)
```

## Arguments

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

## Value

A data.frame. One row for each sample, clock, and panel, with
`n_observed`, `n_needed`, `coverage`, and, when `x` holds more than one
batch, `mc_batch_id`.

## Details

Only the clocks in the returned scores of `x` get a row. A clock that
scores as part of another clock gets none. A clock assembled only from
other clocks' scores gets none. A clock scored separately for each sex
has no row for a sample outside the sex it scored.

A clock that normalizes has a second row for each sample, under
`panel = "norm"`, for the panel used to normalize it.

`samples_coverage()` warns when a row's `coverage` is under the
strictest `min_samples_coverage` value used to score `x`. The
`mc_batch_id` column appears only when `x` holds more than one batch.

## See also

[`clocks_coverage()`](https://hhp94.github.io/methylCIPHERv2/reference/clocks_coverage.md)
for the same panels counted for each clock.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20)
res <- calc_clocks(sim[["DNAm"]], clocks)

head(samples_coverage(res))
#>        id clock_id panel n_observed n_needed coverage
#> 1 sample1 Horvath1 score        353      353        1
#> 2 sample2 Horvath1 score        353      353        1
#> 3 sample3 Horvath1 score        353      353        1
#> 4 sample4 Horvath1 score        353      353        1
#> 5 sample5 Horvath1 score        353      353        1
#> 6 sample6 Horvath1 score        353      353        1
```
