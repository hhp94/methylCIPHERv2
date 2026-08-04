# Epigenetic Clock Catalog

Lists the clocks in the catalog, with the group, tags, and token to
request each one.

## Usage

``` r
list_clocks(group = NULL, tag = NULL, pattern = NULL, all_columns = FALSE)
```

## Arguments

- group:

  A character vector. Keeps only the clocks in these groups. Default is
  `NULL`, which keeps every group.

- tag:

  A character vector. Keeps only the clocks that carry one of these
  tags. Default is `NULL`, which applies no tag filter.

- pattern:

  A string. A regular expression matched against the clock id and the
  group id. Default is `NULL`, which applies no pattern filter.

- all_columns:

  A boolean. Returns every column, including the ones the frame leaves
  out by default. Default is `FALSE`.

## Value

A data.frame. One row for each clock in the catalog, including a clock
that scores only as part of another clock.

## Details

Valid values for `tag` are the names of
[`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md).

`request_as` names the token to pass to `clocks` in
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).
It differs from `clock_id` for a clock that only another clock can
request. `covariates` names the
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
`pheno` columns a clock needs, and `external` is `TRUE` for a clock
whose weights are a download.

`all_columns = TRUE` adds four more columns.

- `callable` is `FALSE` for a clock that only another clock can request.

- `group_size` counts the callable clocks a group token expands to.

- `batch_dependent` is `TRUE` for a clock whose score depends on the
  other samples scored with it.

- `normalize` names the background normalization a clock gets, and is
  empty for a clock that gets none. `"bmiq"` is optional, and you turn
  it on with the `normalize` argument of
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).
  `"quantile"` is part of the clock's definition and is always applied.

## See also

- [`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md)
  for the tags a `tag` value accepts.

- [`clock_cpgs()`](https://hhp94.github.io/methylCIPHERv2/reference/clock_cpgs.md)
  for the CpGs a set of clocks needs.

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets an external clock needs.

## Examples

``` r
list_clocks(pattern = "^Horvath")
#>   clock_id group_id request_as covariates external tags
#> 1 Horvath1 Horvath1   Horvath1               FALSE     
#> 2 Horvath2 Horvath2   Horvath2               FALSE     
nrow(list_clocks(tag = "mortality"))
#> [1] 13
list_clocks(group = "Dunedin", all_columns = TRUE)
#>        clock_id group_id    request_as callable group_size covariates external
#> 1   DunedinPACE  Dunedin   DunedinPACE     TRUE          2               FALSE
#> 2 DunedinPoAm38  Dunedin DunedinPoAm38     TRUE          2               FALSE
#>   batch_dependent normalize tags
#> 1           FALSE  quantile     
#> 2           FALSE               
```
