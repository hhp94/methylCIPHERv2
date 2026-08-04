# Assets Directory Override

Sets or clears the directory this session uses for clock assets.

## Usage

``` r
set_mc_assets_dir(path = NULL)
```

## Arguments

- path:

  A string. The new assets directory. Default is `NULL`, which clears
  the override.

## Value

A string. The override that was in place, returned invisibly, or `NULL`
when no override was set.

## Details

`set_mc_assets_dir()` creates the directory if it does not exist yet. It
stops if the directory cannot be created or is not writable. Clearing
the override falls back to the `mc.assets_dir` option, the
`MC_ASSETS_DIR` environment variable, or the default cache directory, in
that order.

The return value restores the previous state exactly, in the manner of
[`setwd()`](https://rdrr.io/r/base/getwd.html). It is the override that
was in place, and it is `NULL` when no override was set. Passing `NULL`
back clears the override, rather than pinning it to the directory that
happened to be in effect.

For the directory in effect, which is always a path, call
[`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md).

## See also

- [`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md)
  for the directory in effect.

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets in the directory.

- [`download_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/download_mc_assets.md)
  to write the assets to disk.

- [`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md)
  to read the assets into memory.

- [`clear_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/clear_mc_assets.md)
  to delete the assets from disk.

## Examples

``` r
if (FALSE) { # interactive()
old <- set_mc_assets_dir(tempdir())
get_mc_assets_dir()
set_mc_assets_dir(old)
}
```
