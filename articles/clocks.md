# The Clock Catalog

## What the columns mean

**Clock** is the identifier for a single scored clock. **Group** is the
family it belongs to.

**Covariates** lists the columns that `pheno` must carry for that clock.
An empty cell means the clock reads methylation values alone.

**Normalize** names the background normalization scheme a clock uses,
and an empty cell means it uses none. A `quantile` cell is part of the
clock definition and always runs. A `bmiq` cell is optional, and stays
off until the `normalize` argument of
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
turns it on.

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
returns this catalog as a data frame, so it is available in a script as
well as on this page. It takes a tag directly.

``` r

mitotic <- list_clocks(tag = "mitotic")
mitotic[["clock_id"]]
#> [1] "EpiTOC"        "EpiTOC2"       "HypoClock"     "MiAge"        
#> [5] "RepliTali"     "RepliTaliNorm"
```

The result is an ordinary data frame, so you can filter the other
columns with base R.

``` r

all_clocks <- list_clocks()
bundled <- all_clocks[!all_clocks[["external"]], ]
nrow(bundled)
#> [1] 108
```
