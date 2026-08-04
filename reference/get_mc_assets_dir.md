# Active Assets Directory

Returns the directory that this session uses for clock assets.

## Usage

``` r
get_mc_assets_dir()
```

## Value

A string. The absolute path to the assets directory in effect for this
session.

## Details

The directory comes from the `mc.assets_dir` option, the `MC_ASSETS_DIR`
environment variable, or the default cache directory, in that order. Set
an override with
[`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md).

## See also

- [`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md)
  to choose another directory.

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
get_mc_assets_dir()
}
```
