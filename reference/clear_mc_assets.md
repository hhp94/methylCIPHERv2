# Clock Asset Removal

Deletes the clock assets for one or more groups from the assets
directory.

## Usage

``` r
clear_mc_assets(groups = "all", ask = TRUE)
```

## Arguments

- groups:

  A character vector. The asset groups to act on. One or more of
  `"PCBrainAge"`, `"PCClocks"`, `"SystemsAge"` and `"Zhang2019"`, or
  `"all"` for every group. Repeated values are ignored, and an empty
  vector selects nothing. Default is `"all"`.

- ask:

  A boolean. Asks for consent before a download or a delete. Default is
  `TRUE`. Pass `FALSE` to download or delete without asking, in a
  non-interactive session.

## Value

A character vector. The path to each deleted asset file, returned
invisibly. Empty when no file was removed.

## Details

`clear_mc_assets()` removes both the currently declared assets and every
superseded asset left behind by an earlier sync, for the requested
groups. It asks for consent before deleting, refuses in a
non-interactive session, and treats `ask = FALSE` as consent to proceed
without asking. If the assets directory holds no asset for the requested
groups, no file is removed.

## See also

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets on disk and their size.

- [`download_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/download_mc_assets.md)
  to write the assets to disk.

- [`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md)
  to read the assets into memory.

- [`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md)
  for the directory in effect.

- [`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md)
  to choose another directory.

## Examples

``` r
if (FALSE) { # interactive()
clear_mc_assets("SystemsAge")
}
```
