# Clock Asset Downloads

Downloads the clock assets for one or more groups into the assets
directory.

## Usage

``` r
download_mc_assets(groups = "all", ask = TRUE)
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

A character vector. The path to each requested asset file, named by
group id. Returned invisibly.

## Details

Only an asset that is missing from the assets directory is fetched.
`download_mc_assets()` asks for consent before a download, refuses in a
non-interactive session, and treats `ask = FALSE` as consent to proceed
without asking.

## See also

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  for the assets already on disk.

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
download_mc_assets("SystemsAge")
download_mc_assets("SystemsAge", ask = FALSE)
}
```
