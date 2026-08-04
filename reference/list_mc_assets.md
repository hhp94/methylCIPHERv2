# Clock Asset Inventory

Lists the declared clock assets and their status in the assets
directory.

## Usage

``` r
list_mc_assets(groups = "all")
```

## Arguments

- groups:

  A character vector. The asset groups to act on. One or more of
  `"PCBrainAge"`, `"PCClocks"`, `"SystemsAge"` and `"Zhang2019"`, or
  `"all"` for every group. Repeated values are ignored, and an empty
  vector selects nothing. Default is `"all"`.

## Value

A data.frame. One row for each requested group, with its clock count,
CpG count, asset size, download status, and the count and total size of
its superseded assets.

## See also

- [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md)
  for the clocks in the catalog.

- [`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md)
  for the tags a `clocks` value accepts.

- [`clock_cpgs()`](https://hhp94.github.io/methylCIPHERv2/reference/clock_cpgs.md)
  for the CpGs a set of clocks needs.

- [`download_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/download_mc_assets.md)
  to write the assets to disk.

- [`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md)
  to read the assets into memory.

- [`clear_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/clear_mc_assets.md)
  to delete the assets from disk.

- [`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md)
  for the directory in effect.

- [`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md)
  to choose another directory.

## Examples

``` r
if (FALSE) { # interactive()
list_mc_assets()
list_mc_assets("SystemsAge")
}
```
