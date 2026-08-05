# Loaded Clock Assets

Reads the assets for one or more clock groups into memory.

## Usage

``` r
load_mc_assets(groups, ext_data = NULL, ask = TRUE)
```

## Arguments

- groups:

  A character vector. The asset groups to load. One or more of
  `"PCBrainAge"`, `"PCClocks"`, `"SystemsAge"` and `"Zhang2019"`, or
  `"all"` for every group. Repeated values are ignored, and an empty
  vector loads nothing.

- ext_data:

  A string. The path to the directory that holds the clock assets.
  Default is `NULL`, which uses the assets directory.

- ask:

  A boolean. Asks for confirmation before the assets directory changes.
  Default is `TRUE`. Pass `FALSE` to continue without asking, in a
  non-interactive session.

## Value

An `mc_assets` object. It holds the loaded asset for each requested
group, in the order of `groups`.

## Details

`load_mc_assets()` asks for confirmation before it downloads an asset
that is not already in the assets directory. It refuses to download in a
non-interactive session. `ask = FALSE` confirms the download in advance.

An asset in `ext_data` for a group that was not requested is not used,
and `load_mc_assets()` warns about it.

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

- Assets already in memory from `load_mc_assets()` are used directly.

## See also

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets already on disk.

- [`download_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/download_mc_assets.md)
  to write the assets to disk.

- [`clear_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/clear_mc_assets.md)
  to delete the assets from disk.

- [`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md)
  for the directory in effect.

- [`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md)
  to choose another directory.

## Examples

``` r
if (FALSE) { # interactive()
assets <- load_mc_assets("SystemsAge")
names(assets)
}
```
