# Clock Tag Registry

Lists the keyword tags that expand to a group of clocks, for use as a
`clocks` or `tag` value.

## Usage

``` r
list_clock_tags()
```

## Value

A named list. Each name is a tag, and each element is a character vector
of the group or clock tokens it expands to.

## See also

- [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md)
  for the clocks a tag expands to.

- [`clock_cpgs()`](https://hhp94.github.io/methylCIPHERv2/reference/clock_cpgs.md)
  for the CpGs a set of clocks needs.

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets an external clock needs.

## Examples

``` r
list_clock_tags()
#> $gestational
#> [1] "Bohlin"          "Knight"          "Mayne"           "LeePlacentalAge"
#> 
#> $mitotic
#> [1] "EpiTOC"    "EpiTOC2"   "MiAge"     "RepliTali"
#> 
#> $mortality
#> [1] "GrimAge"        "ZhangMortality"
#> 
```
