# Epigenetic Clock Scores

Scores CpG-based epigenetic clocks on a matrix of methylation beta
values.

## Usage

``` r
calc_clocks(
  DNAm,
  clocks,
  pheno = NULL,
  pheno_id = "ID",
  min_clocks_coverage = 0.75,
  min_samples_coverage = 0.75,
  normalize = NULL,
  ext_data = NULL,
  ask = TRUE
)
```

## Arguments

- DNAm:

  A numeric matrix. The methylation beta values, with samples in the
  rows and CpGs in the columns.

- clocks:

  A character vector. The clocks to score, named by clock id, group id,
  or tag. See
  [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md).

- pheno:

  A data.frame. The sample metadata, with one row for each sample.
  Default is `NULL`.

- pheno_id:

  A string. The name of the column in `pheno` that holds the sample ids.
  Default is `"ID"`.

- min_clocks_coverage:

  A number between 0 and 1. The smallest fraction of a clock's CpGs that
  must be present for that clock to score. Default is `0.75`.

- min_samples_coverage:

  A number between 0 and 1. The smallest fraction of a clock's CpGs that
  must be present for a sample to score without a warning. Default is
  `0.75`.

- normalize:

  A named logical vector. Turns background normalization on for the
  clocks that support it. Default is `NULL`, which leaves the optional
  schemes off.

- ext_data:

  A string. The path to the directory that holds the clock assets.
  Default is `NULL`, which uses the assets directory.

- ask:

  A boolean. Asks for confirmation before the assets directory changes.
  Default is `TRUE`. Pass `FALSE` to continue without asking, in a
  non-interactive session.

## Value

An `mc_result` object. It holds the scores, the narrowed `pheno`, the
coverage counts, and the provenance of the run.

## Details

[`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md)
and
[`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md)
show every value `clocks` accepts.

`normalize` turns on the schemes that a clock declares as optional. It
cannot turn off a scheme that is part of the clock. The `normalize`
column of `list_clocks(all_columns = TRUE)` gives the scheme each clock
uses.

The two coverage arguments differ. `calc_clocks()` stops when a clock
has too few CpGs present, so every clock in the returned scores passed
`min_clocks_coverage`. A clock just above that floor, and a clock whose
normalization panel falls below it, each raise a warning and still
score. A sample with too few CpGs present raises a warning and still
scores. Pass the returned value to
[`clocks_coverage()`](https://hhp94.github.io/methylCIPHERv2/reference/clocks_coverage.md)
or
[`samples_coverage()`](https://hhp94.github.io/methylCIPHERv2/reference/samples_coverage.md)
to see the counts.

`calc_clocks()` narrows `pheno` before it stores it. The returned value
keeps the id column and the covariates that the clocks need, and drops
the other columns.

## The assets directory

Four clock groups keep their weights in downloadable assets, outside the
package. `ext_data` says where to read them from, and accepts three
forms.

- `NULL` reads from the assets directory, and downloads any asset that
  is missing. Use
  [`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md)
  to choose that directory.

- A path reads only that directory, and never downloads. A missing asset
  is an error.

- Assets already in memory from
  [`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md)
  are used directly.

## Examples

``` r
clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 20)

res <- calc_clocks(sim[["DNAm"]], clocks)
res
#> <mc_result> 20 sample(s) x 2 clock(s)
#> 
#> $scores [6 of 20 row(s), 2 of 2 clock(s)]
#>          Horvath1   Hannum
#> sample1 121.57529 37.88488
#> sample2 131.05186 50.32948
#> sample3 102.20911 60.26443
#> sample4  46.72342 69.17409
#> sample5  73.39157 62.24323
#> sample6 127.55273 69.88226
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

# pheno is narrowed to the id column and the covariates the clocks need
pheno <- data.frame(ID = rownames(sim[["DNAm"]]), Age = runif(20, 20, 80))
res <- calc_clocks(sim[["DNAm"]], clocks, pheno = pheno)
head(res[["pheno"]])
#>        ID
#> 1 sample1
#> 2 sample2
#> 3 sample3
#> 4 sample4
#> 5 sample5
#> 6 sample6
```
