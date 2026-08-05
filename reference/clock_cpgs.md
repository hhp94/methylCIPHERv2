# CpGs Required To Score Clocks

Lists the CpGs needed to score a set of clocks, including background
CpGs for normalization.

## Usage

``` r
clock_cpgs(clocks, normalize = NULL, ext_data = NULL, ask = TRUE)
```

## Arguments

- clocks:

  A character vector. The clocks to score, named by clock id, group id,
  or tag. See
  [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md).

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

A character vector. The CpGs needed to score `clocks`, with duplicates
removed.

## Details

A clock built from other clocks also needs their CpGs. A clock whose
normalization is part of its definition always adds its background
panel. For a clock whose normalization is optional, `normalize` adds the
background panel only when it is on.

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

## See also

- [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md)
  for the clocks a `clocks` value accepts.

- [`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md)
  for the tags a `clocks` value accepts.

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets an external clock needs.

## Examples

``` r
cpgs <- clock_cpgs(c("Horvath1", "Hannum"))
length(cpgs)
#> [1] 418

# normalizing Horvath1 adds its background panel to the union
norm_cpgs <- clock_cpgs(c("Horvath1", "Hannum"), normalize = c(Horvath1 = TRUE))
length(norm_cpgs)
#> [1] 21432
```
