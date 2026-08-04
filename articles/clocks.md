# The Clock Catalog

``` r

library(methylCIPHERv2)
```

## What the columns mean

**Clock** is the identifier for a single scored clock, and the name to
pass to
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).
**Group** is the family it belongs to, and a group name can be requested
on its own to score every member.

**Covariates** lists the columns that `pheno` must carry for that clock.
An empty cell means the clock reads methylation values alone.

## Scoring from the catalog

Any value in the **Clock** column works as the `clocks` argument.

``` r

res <- calc_clocks(DNAm, clocks = "Horvath1", pheno = pheno)
```

A group name scores every member of that family.

``` r

res <- calc_clocks(DNAm, clocks = "GrimAge", pheno = pheno)
```

A tag scores every clock that carries it.
[`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md)
returns the tags with their members.

``` r

names(list_clock_tags())
#> [1] "gestational" "mitotic"     "mortality"
```

## Filtering the catalog in code

[`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md)
returns the table above as a data frame, so the catalog is available in
a script as well as on this page. It takes a tag directly.

``` r

mitotic <- list_clocks(tag = "mitotic")
mitotic[["clock_id"]]
#> [1] "EpiTOC"        "EpiTOC2"       "HypoClock"     "MiAge"        
#> [5] "RepliTali"     "RepliTaliNorm"
```

The result is an ordinary data frame, so the other columns filter with
base R.

``` r

catalog <- list_clocks()
bundled <- catalog[!catalog[["external"]], ]
nrow(bundled)
#> [1] 108
```
