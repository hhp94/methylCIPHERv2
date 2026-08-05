# Shared Parameters

Parameter definitions reused across the package.

## Arguments

- DNAm:

  A numeric matrix. The methylation beta values, with samples in the
  rows and CpGs in the columns.

- x:

  An `mc_result` object. The value returned by
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- clocks:

  A character vector. The clocks to score, named by clock id, group id,
  or tag. See
  [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md).

- pheno:

  A data.frame. The sample metadata, with one row for each sample.
  Default is `NULL`.

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

- all_columns:

  A boolean. Returns every column, including the ones the frame leaves
  out by default. Default is `FALSE`.

- groups:

  A character vector. The asset groups to act on. One or more of
  `"PCBrainAge"`, `"PCClocks"`, `"SystemsAge"` and `"Zhang2019"`, or
  `"all"` for every group. Repeated values are ignored, and an empty
  vector selects nothing. Default is `"all"`.

- long:

  A boolean. Returns one row for each sample and clock when `TRUE`, and
  one row for each sample, with one column for each clock, when `FALSE`.
  Default is `TRUE`.

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

## Clocks that use all the samples

Some clocks depend on information from all the samples, such as a
z-score. When `x` holds more than one batch, these clocks take their
value from every sample in `x`, and not from one batch alone. This is
the same calculation as
[`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md).
