# Simulated Methylation Data

Builds a random beta matrix and a matching `pheno` data.frame for a set
of clocks.

## Usage

``` r
sim_DNAm(
  clocks,
  n = 10,
  Age = FALSE,
  Female = FALSE,
  remove = 0,
  normalize = NULL,
  ext_data = NULL,
  ask = TRUE,
  suffix = NULL
)
```

## Arguments

- clocks:

  A character vector. The clocks to score, named by clock id, group id,
  or tag. See
  [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md).

- n:

  A single whole number. The number of samples to simulate. Default is
  `10`.

- Age:

  A boolean. Adds an `Age` column to `pheno`, drawn from a normal
  distribution. Default is `FALSE`.

- Female:

  A boolean. Adds a `Female` column to `pheno`, with about half the
  samples set to `1`. Default is `FALSE`.

- remove:

  A single whole number. The number of CpGs to drop at random from the
  simulated panel. Default is `0`.

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

- suffix:

  A string. Appended to every sample id, so two simulated matrices stay
  distinct. Default is `NULL`, which leaves the ids as given.

## Value

An `mc_sim` object. It holds the simulated `DNAm` matrix, the matching
`pheno` data.frame, the `clocks` argument as given, and the `suffix`,
which is `NULL` when no suffix was set.

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
sim <- sim_DNAm(c("Horvath1", "Hannum"), n = 10, Age = TRUE, Female = TRUE)
dim(sim[["DNAm"]])
#> [1]  10 418
head(sim[["pheno"]])
#>        ID      Age Female
#> 1 sample1 44.13629      0
#> 2 sample2 47.61322      0
#> 3 sample3 44.30379      1
#> 4 sample4 52.24111      0
#> 5 sample5 40.72297      0
#> 6 sample6 58.15295      1
```
