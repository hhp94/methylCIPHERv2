# Package index

## Scoring

Score a beta matrix.

- [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
  : Epigenetic Clock Scores
- [`predict_sex()`](https://hhp94.github.io/methylCIPHERv2/reference/predict_sex.md)
  : Predicted Sex Karyotype

## Discovery

Read the clock catalog before a run.

- [`list_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clocks.md)
  : Epigenetic Clock Catalog
- [`list_clock_tags()`](https://hhp94.github.io/methylCIPHERv2/reference/list_clock_tags.md)
  : Clock Tag Registry
- [`clock_cpgs()`](https://hhp94.github.io/methylCIPHERv2/reference/clock_cpgs.md)
  : CpGs Required To Score Clocks

## Coverage

Report the panel coverage for each clock and for each sample.

- [`clocks_coverage()`](https://hhp94.github.io/methylCIPHERv2/reference/clocks_coverage.md)
  : Clock Coverage Counts
- [`samples_coverage()`](https://hhp94.github.io/methylCIPHERv2/reference/samples_coverage.md)
  : Sample Coverage Counts

## Working with an mc_result

Convert, combine and print the value that
[`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md)
returns.

- [`as.data.frame(`*`<mc_result>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/as.data.frame.mc_result.md)
  : Data Frame Method For An mc_result Object
- [`as.matrix(`*`<mc_result>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/as.matrix.mc_result.md)
  : Matrix Method For An mc_result Object
- [`rbind(`*`<mc_result>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/rbind.mc_result.md)
  : Combined Batches Of Scores
- [`refinalize_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/refinalize_clocks.md)
  : Scores Recomputed From All Samples
- [`print(`*`<mc_result>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/print.mc_result.md)
  : Print Method For An mc_result Object

## Analysis

Fit age acceleration and test associations.

- [`calc_accel()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_accel.md)
  : Age Acceleration Or Difference
- [`score_associations()`](https://hhp94.github.io/methylCIPHERv2/reference/score_associations.md)
  : Clock Age Associations

## Assets

Manage the downloadable weights that some clock groups need.

- [`list_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/list_mc_assets.md)
  : Clock Asset Inventory
- [`get_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/get_mc_assets_dir.md)
  : Active Assets Directory
- [`set_mc_assets_dir()`](https://hhp94.github.io/methylCIPHERv2/reference/set_mc_assets_dir.md)
  : Assets Directory Override
- [`download_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/download_mc_assets.md)
  : Clock Asset Downloads
- [`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md)
  : Loaded Clock Assets
- [`clear_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/clear_mc_assets.md)
  : Clock Asset Removal
- [`print(`*`<mc_assets>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/print.mc_assets.md)
  : Print Method For An mc_assets Object

## Citations

Collect the references for the clocks that a run scored.

- [`cite_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/cite_clocks.md)
  : Clock Citations
- [`as.data.frame(`*`<mc_citation>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/as.data.frame.mc_citation.md)
  : Data Frame Method For An mc_citation Object
- [`toBibtex(`*`<mc_citation>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/toBibtex.mc_citation.md)
  : Bibtex Method For An mc_citation Object
- [`print(`*`<mc_citation>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/print.mc_citation.md)
  : Print Method For An mc_citation Object

## Simulation

Generate a beta matrix for an example or a test.

- [`sim_DNAm()`](https://hhp94.github.io/methylCIPHERv2/reference/sim_DNAm.md)
  : Simulated Methylation Data
- [`print(`*`<mc_sim>`*`)`](https://hhp94.github.io/methylCIPHERv2/reference/print.mc_sim.md)
  : Print Method For An mc_sim Object
