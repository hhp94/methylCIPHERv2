# Managing Clock Assets

Most clocks are included with the package and need no setup. Some clock
groups are too large, so their weights live in downloadable files (i.e.,
assets).

This vignette shows how
[methylCIPHERv2](https://github.com/hhp94/methylCIPHERv2) resolves the
assets directory, and how to change it.

``` r

library(methylCIPHERv2)
```

## Discover

[`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
shows every clock group that needs a downloaded asset.

``` r

list_mc_assets()
#>     group_id n_clocks n_cpgs   size downloaded superseded superseded_size
#> 1 PCBrainAge        1 357852  6.68M      FALSE          0               0
#> 2   PCClocks       14  78464  8.88M      FALSE          0               0
#> 3 SystemsAge       13 125175 22.46M      FALSE          0               0
#> 4  Zhang2019        1 319607  5.02M      FALSE          0               0
```

`downloaded` says whether the current asset is on disk. `superseded`
counts older assets for that group, which a later version of the package
replaced. Superseded assets waste space and nothing reads them.

## Pathing

[`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md)
returns the path in effect.

``` r

get_mc_assets_dir()
#> "C:/Users/<you>/AppData/Local/R/cache/R/methylCIPHERv2"
```

To override that path for this session, run
[`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md),
which returns the previous setting. The same function restores the
default. This is useful for an HPC workflow, where the path can point at
a folder of assets you uploaded ahead of time.

``` r

old <- set_mc_assets_dir("/your/override")

# restore the previous setting
set_mc_assets_dir(old)

# or clear the override and fall back to the default
set_mc_assets_dir(NULL)
```

## Durable override

The default path is a cache directory, which a cleanup tool or the
operating system can wipe (more likely on macOS and Linux). For a
durable path, set the `MC_ASSETS_DIR` environment variable in a user or
project `.Renviron` file.

``` r

# open the user .Renviron file
file.edit("~/.Renviron")
```

Add this line to it, then restart R.

    MC_ASSETS_DIR=/your/override

## Downloading

[`download_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/download_mc_assets.md)
writes the assets to disk, and always asks for confirmation. Use
`ask = FALSE` in a script, a scheduled job, or a container build.
Consider the size of the download before you do.

``` r

download_mc_assets("SystemsAge", ask = TRUE) # the default, prompts first

# calc_clocks() routes through the same mechanism. ask = FALSE
# authorizes a download when one is needed.
calc_clocks(DNAm, clocks = "SystemsAge", ask = FALSE)
```

Pass `"all"` for every group, or a vector of group names for some of
them.

``` r

download_mc_assets(c("SystemsAge", "PCClocks"))
```

## Loading an asset into memory

Assets can be pre-loaded into memory once with
[`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md),
then re-used for the whole session. Pass the loaded object to the
`ext_data` argument of
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).
This reduces the overhead of a repeated scoring workflow.

``` r

assets <- load_mc_assets("SystemsAge")
res <- calc_clocks(DNAm, "SystemsAge", ext_data = assets)
```

## A directory that never downloads

The `ext_data` argument also takes a path. A path restricts the run to
the assets already in that directory. Nothing is downloaded, and a
missing asset is an error rather than a prompt.

``` r

res <- calc_clocks(DNAm, "SystemsAge", ext_data = "/data/methylation-assets")
```

Use this form on a machine with no network access, or in a pipeline that
must not reach the internet. Leave `ext_data` at its default of `NULL`
to allow a download.

## Deleting the assets

[`clear_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/clear_mc_assets.md)
empties the directory. It lists the files and asks for confirmation, and
it takes `ask` like the two verbs above.

``` r

clear_mc_assets()
```

The delete covers the current assets and the superseded ones together.
Pass a group name to delete one group.

``` r

clear_mc_assets("SystemsAge")
```

A file the package did not write is never touched.
