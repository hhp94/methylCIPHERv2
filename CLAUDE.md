# CLAUDE.md

Guidance for Claude Code in this repo. This file holds **invariants** that rarely change.
Volatile detail (per-clock status, exact designs, dated reversals) lives in `dev/` -- see
"Source-of-truth docs" and prefer it for specifics.

## What this package is

`methylCIPHERv2` scores CpG-based DNA-methylation ("epigenetic clock") ages. One public scorer,
`calc_clocks()`, drives everything. The scoring contract (clock catalog + coefficient tensors) is
synced from the separate `methylCIPHER-meta` repo; fixtures are the scientific gate. Target is
**CRAN**, not Bioconductor. R (>= 4.4).

## Getting started (collaborators)

`R/sysdata.rda` (the compiled catalog) is **committed**, so you can develop, load, and test with no
meta repo and no downloads. Only regenerating it needs `sync()` (below), which you do not need to
contribute.

```r
# from the package root, in R:
install.packages("pak")
pak::local_install_deps(dependencies = TRUE)  # reads DESCRIPTION incl. GitHub-only Remotes
pkgbuild::compile_dll(".", force = TRUE)  # only after editing src/*.cpp -- load_all() reuses a stale dll
devtools::load_all()   # attach for interactive work
devtools::document()   # regenerates NAMESPACE + man/ from tags; needed for the Rcpp wiring
devtools::test()       # always-on tiers (cohort parity auto-skips if not staged)
```

`devtools::check()` / `R CMD check` is **maintainer-on-demand only** -- see the invariant below.

The package has compiled code (`src/`), so a working toolchain is required -- Rtools on Windows.

Soft deps (`betanorm`, `duckdb`, `DBI`) back specific paths only; tests skip when absent.

## Non-negotiable invariants

Do not reverse these without a `dev/DECISIONS.md` entry explaining why.

- **Never run `R CMD check` / `devtools::check()`. It is on-demand, maintainer-only.** It hangs in
  this environment: the check re-runs the full suite in a fresh install, and the suite is currently
  bloated enough that the run does not finish in a usable time. An agent that starts one blocks the
  session on something it cannot resolve. **Verify a change with `devtools::test()`** (always-on
  tiers) and say plainly that check was not run. This is the same rule as the parity tier -- see
  "Testing" -- and for the same reason: minutes-to-hours of maintainer wall-clock is the
  maintainer's call, not the agent's. Do not reach for it via `Rscript`, `pkgbuild`, or a background
  shell either; the prohibition is on the work, not on one entry point.

- **One engine + a finite, closed branch set.** Every work unit routes on the catalog pair
  `(weights_format, computation_type)` to shared `linear_score()` or a named branch (pre-transform,
  family orchestrator, sex-routed alias, external, custom). There is **no** recipe
  interpreter/walker, not even as a fallback.
- **Routing is total, and a gap is a hard stop.** `score_type()` returns a known tag for every
  catalog entry or `stop()`s naming the clock's group / `weights_format` / `computation_type`.
  There is no `"unsupported"` tag and nothing filters on one: a sync that adds a routing pair no
  branch claims must fail the always-on tier, never silently shrink it (see DECISIONS 2026-07-23).
- **A branch returns only its score** (`<n x 1 matrix>`), never a coverage record. The one exception
  is declared, not per-branch: a clock the catalog marks cohort-reducing (`cross_sample_at`, today
  only the two `DNAmPhysAge` clocks) returns its **per-sample intermediate** instead, which the
  scoring loop routes into `pending` rather than `scores`. `finalize_cross_sample()` performs the
  reduction once, after assembly, and every front end calls it unconditionally -- it is a no-op when
  `pending` is empty. Nothing tests for a clock id or for whether the run was chunked; the routing
  reads `spec$cross_sample` (DECISIONS 2026-07-27, "Phase 3"). Coverage/QC
  depends on no score, so it is computed once upstream of the scoring loop by `compute_coverage()`
  (`R/coverage.R`), keyed by clock id, and merged in `construct_mc_result()`. Per-sample miss is
  counted **once per distinct panel** (FitAge/GrimAge reuse panels) and kept **per panel role**:
  every clock has a score panel; a normalizing clock (only DunedinPACE) also has a norm panel. The
  counts follow the policy uniformly -- `vendor_mean` fills every absent CpG into the predictor
  (`used = present + imputed_full`), anything else drops them (`used = present`, `dropped = absent`).
  Routed members are masked to the samples they scored. `coverage_record()` (`R/coverage.R`) owns
  the per-clock record fields --
  including `normalizes` (**the one declared panel fact**; readers must not re-derive it from
  `norm_needed`) and the per-panel `score_imputed_partial` / `norm_imputed_partial`.
  **The record is one axis: every count in it is a CpG count**, the partial ones included --
  `score_imputed_partial` is how many of the panel's present CpGs were cohort-mean filled for any
  sample, so it reads against `score_needed` / `score_present` like its neighbours and is
  block-invariant (the fill columns are cohort facts from pass 1). Never store the cell sum
  `sum(sample_miss)` there: it collapses both axes at once, is neither a probe fact nor a sample
  fact, and can exceed `score_present`. **The sample axis is `$coverage$sample_miss`, and only
  that** (DECISIONS 2026-07-29).
  `$coverage$sample_miss` is `list(score = <n x k matrix>, norm = <n x k' matrix over just the
  normalizing columns>)`. **`$per_clock`, `$sample_miss` and `samples_coverage()` span one set --
  the clocks that read CpGs** -- so `k` is neither the returned columns nor the whole sequence:
  routing targets are in, pure composites are out even when returned (DECISIONS 2026-07-29,
  "pure composites"). `cached_cols()` / `count_sample_miss()` (integer) / `score_matrix()`
  remain the shared shape helpers (see DECISIONS 2026-07-24).
- **Result is an S3 record over `list`** (class `mc_result`): `$scores` (n x k double), `$pheno`,
  `$coverage`, `$provenance`. Never a `matrix` subclass (drops class + attrs on first subset).
  **Where a verb exists it is a method**, and the built surface today is exactly `print`,
  `as.matrix` and `cite_clocks` -- plus `rbind`, which exists in order to **refuse**. Coverage is
  deliberately not a method: it is the plain `clocks_coverage()` / `samples_coverage()`, **not**
  `summary()`. Citations dispatch as `cite_clocks()` -- a **package-owned** generic, because both
  `utils::citation` and `utils::cite` already exist as plain functions and taking either name
  masks it (DECISIONS 2026-07-23, 2026-07-24, 2026-07-25). Nothing else is promised: `as.data.frame`,
  `[`, `cbind`, `augment` and `codebook` were listed here for a year without being written, so
  they are **unbuilt ideas, not contracts** -- adding one is a new API decision, and until a human
  makes it there is no behaviour to match (DECISIONS 2026-07-27).
- **Scores only, and the record remembers its inputs.** `$scores` is scores -- no auto-appended
  phenotype columns. Separately, `$pheno` carries the *aligned* pheno narrowed to the id column
  plus the covariates the run actually required, so a saved record can answer what was fed in.
  Align pheno by sample id, never row order. If a frame conversion is ever built, it inherits this
  rule: scores plus the id column, and pheno stays off that path so it cannot leak by accident.
- **Imputation in one place, never crossing sources.** Partial NA on a present probe -> cohort mean
  (shared cache); a fully absent probe -> the clock's vendored ref, or drop by policy.
- **Accessors are the executable schema.** `calc_clocks` consumes accessors (`get_clock`,
  `clock_coefs`, ...), never raw nested catalog lists. No hand-written `schema.md`.
- **Never `$` in `R/`. Always `[[`.** `$` partial-matches on lists, so a missing exact field
  silently resolves to a longer one (`entry$covariates` -> `covariates_required`) and the caller
  gets a wrong value, not an error. The rule is **blanket, not scoped to catalog/pack/tensor
  reads**: a scoped rule needs a judgement call at every site about what the object is, and the
  list that is "obviously not a catalog entry" today is the one someone widens tomorrow. So it
  binds on every `$` in `R/` alike -- catalog entries, result records, `optim()` output, and
  environments, where `$` is exact and harmless but is not worth the exception.
  `tests/testthat/test-source-hygiene.R` enforces it by scanning **parse tokens**, not text, so a
  `$` inside a comment or a `"\\.qs2$"` regex does not count and a real access cannot hide; the
  only exempt file is the generated `R/RcppExports.R`. `options(warnPartialMatchDollar)` is
  **not** the fix -- a package cannot set a session global for its users and it does not fire
  under `R CMD check`.
- **Never `<<-` in `R/`. Mutable state is an explicit environment.** `<<-` does not name a target:
  it walks the enclosing frames and assigns into the first one that already binds that name,
  **creating a global** if none does. So renaming or deleting a local silently promotes a local
  update to a package/global one -- no error, and nothing at the read site says where the value
  came from. The replacement is a store you can point at: `st <- new.env(parent = emptyenv())`,
  write `st[["x"]] <- ...`, and have the mutator **return `st`** so readers take it as a value
  instead of inheriting it from a frame. Reads and writes on that env go through `[[` like
  everything else -- the rule above is blanket. Precedent: `miage_fit()` in `R/score_MiAge.R`,
  whose L-BFGS-B objective and gradient share one cached `b^(n-1)` that way. Outside `R/` this is
  a preference, not a rule: `data-raw/sync.R` and a condition collector in `tests/` still use
  `<<-` legitimately.
- **Accessors read declarations; they never search.** No `grep`/regex/fuzzy match over tensor
  names, clock ids, or file paths to find a payload. Resolve the declared pointer (component,
  `probe_sets[["scoring"]][["file"]]`, `imputation[["ref"]]`) and `stop()` when it is absent or
  ambiguous -- an accessor that cannot find its declaration has done its job by failing. Searching
  hides an upstream/sync gap and can silently return a sibling clock's tensor.
  **Three of those lists are keyed by sync, so the lookup is a lookup.** `components` is named by
  `name`, `probe_sets` by `role`, `recipe` by `out` (steps without one stay unnamed), and
  `key_declarations()` stops the build on a collision -- which is why the accessors reading them
  carry no multiplicity guard. The two remaining scans are the two whose predicate is not a key:
  `component_tensor()` on `row_key` and `stack_step()` on `op`, both of which genuinely collide in
  the shipped catalog, so `pick_one()` stays load-bearing there (DECISIONS 2026-07-28).
- **Coverage is never reported for a sample it is not true of, and a clock reports only what it
  counted itself.** A clock whose branch reads no betas -- it is assembled purely from other
  clocks' scores -- **has no coverage of its own**: `per_clock[[id]]` is `NULL`, it gets no
  `sample_miss` column and no `samples_coverage()` row, and its all-`NA` `clocks_coverage()` row
  says so. `clock_reads_cpgs()` (`R/score_cohort.R`) is the one source; it switches on
  `score_type()`, so it is a fact about the closed branch set, not a clock list. Today it selects
  the 7 sex-routed aliases, `GrimAgeV1` and `DNAmFitAge_{Sex}` -- `GrimAgeV2` keeps its record
  because its cox stack declares `internal` surrogates and it really does read its 1030 CpGs.
  **This loses nothing**: a clock that reads no betas can only be fed through its dependencies, so
  every CpG in its declared panel is already counted on a descendant that does read one (verified:
  `panel \ leaf-closure` is empty for all 10). The rule replaces the older "records coverage iff
  every component contributes to every sample", which let `DNAmFitAge` report 100% coverage and
  `score_dropped = 0` while 613 CpGs of its `GrimAgeV1` component were dropped -- its declared
  172-CpG panel is a strict subset of the 1201 that feed it (DECISIONS 2026-07-29). Do not "fill
  in" a `NULL` record with a merged figure, and do not restore a stitched per-sample count for an
  alias: read the descendants' rows for the denominators. **`samples_coverage()` drops NA-coverage
  rows**, which has exactly one source -- a routed member masked on a row its sex did not score --
  so the long frame carries one row per sample per family, under the model that scored it. A
  sample no model scored (unknown sex) therefore has no row at all, which is the same fact its
  `NA` score already carries. The frame was never a complete sample x clock grid; do not make it
  one by keeping rows that assert nothing.
  **The converse binds too: never score a CpG coverage did not count.** Coverage counts the
  *declared* panel, so a branch takes its CpGs from the resolved `cpgs` it is handed
  (`score_present` / `score_absent`) and never re-derives them against the block's cohort-wide
  `usable` set -- that set is every panel's union, so a stray coefficient would resolve against it
  and be scored silently. A branch holding a bare coef vector (GrimAgeV2 surrogates, PhysAge
  surrogates) goes through `component_present()` (`R/score_default.R`), which intersects with the
  declared panel and `stop()`s if the component names a CpG the panel does not -- the sync gap that
  would otherwise be invisible. The always-on smoke tier scores every bundled clock, so it exercises
  that guard over the shipped catalog without naming a clock (DECISIONS 2026-07-29).
- **The callable pool is not the catalog, and neither is the output.** Clocks that exist only as
  routing targets (the 14 sex-resolved DNAmFitAge members) are internal machinery: scored, kept
  for coverage, **never a score column**, and a hard error if requested by name, pointing at their
  alias. A sex-routed family returns exactly one column per alias, populated for every sample --
  never a male column, a female column and NAs. The pool, the refusal, its suggestion and the
  output filter (`drop_routed_members()`) all derive from one source (`sex_routed_members()`), so
  they cannot drift.
- **No network at install/build/check/CRAN test.** Double-precision coefficients only.
- **No commit SHA / pin as result provenance.** Correctness is proven by fixtures.
- **Correlation is never a numeric gate. Anywhere, for anything.** Not in parity, not in a unit
  test, not as a "sanity check" alongside a real bound. `cor()` is offset- and scale-invariant, so
  it cannot distinguish "we match the oracle" from "we are uniformly wrong": a constant intercept
  shift of any size and an arbitrary rescaling both score ~1.0. It also concentrates near 1 for any
  monotone-ish agreement, which hides catastrophic per-sample outliers in a cloud of correct
  points. A numeric agreement gate is always a **bounded per-element difference** -- absolute and
  relative, both, taken as a `max` -- so one bad sample fails the test. The same reasoning bans
  `median`/`mean` as the reducer over per-element differences: any statistic that averages away a
  minority of arbitrarily-wrong samples is the same bug wearing a different name.
- **Roxygen is on only where compiled code forces it; prose docs are still deferred.** The package
  went Rcpp, and `useDynLib` has no route into `NAMESPACE` except a roxygen tag -- so roxygen is
  enabled, `NAMESPACE` and `man/*.Rd` are **generated files**, and `devtools::document()` is a
  normal part of the workflow. That is the whole of the override: keep writing short `#` comments
  (see "Comments") and do **not** start authoring real `@param` / `@return` prose. `calc_clocks()`
  carries a placeholder block whose params are literally `x`; leave it that way until a human says
  otherwise (DECISIONS 2026-07-27).
- **Never hand-edit `NAMESPACE` or `man/*.Rd` -- own the tags, not the files.** They carry roxygen's
  "do not edit by hand" header and `document()` rewrites them from tags, **silently dropping**
  anything added by hand. This has bitten once: a hand-added
  `useDynLib(methylCIPHERv2, .registration = TRUE)` vanished on the next `document()` and took every
  compiled kernel down with it (`.Call()` -> "not available for .Call() for package"). So a new
  export or S3 method gets a bare `#' @export` beside the function (precedent in `R/mc_result.R`),
  and package-level wiring -- `@useDynLib`, `@importFrom` -- lives in
  `R/methylCIPHERv2-package.R`. Then run `document()` and check the diff is only what you intended.
  The maintainer still owns the exported *surface*: say in your summary which `export()` /
  `S3method()` entries a change implies rather than quietly widening it.
- **A `.cpp` edit needs an explicit rebuild.** `devtools::load_all()` happily reuses a stale
  `src/*.dll`, so a changed or newly added kernel silently does not exist. Run
  `pkgbuild::compile_dll(".", force = TRUE)` first; symptom of skipping it is the same
  "not available for .Call()" error above.

## sync.R workflow (`data-raw/sync.R`)

Pulls the scoring contract from `methylCIPHER-meta` into the package. Not run at build/check -- a
maintainer runs it and commits the regenerated `R/sysdata.rda`. **You do not need this to
contribute** (the catalog is committed). `sync()` needs read access to `methylCIPHER-meta`
(private, pre-release); `sync(upload = TRUE)` also needs a release-write token (maintainer-only).

- **Remote:** `https://github.com/hhp94/methylCIPHER-meta.git`.
- **Inputs R may read:** `manifest.json`, `weights/**`,
  `bibliography/{clock_citations.csv,clocks.bib}`. **Never** `control/`, `papers/`, `scripts/`, or
  `bibliography/papers.csv`.
- **Entry point:** `sync(source_git_sha = NULL, upload = FALSE, force = FALSE)`.
  1. Resolve + checkout meta at `source_git_sha` (clone under `data-raw/methylCIPHER-meta/`).
  2. **Always** rebuild catalog + accessor objects + small bundles -> `R/sysdata.rda` (~2s, no
     build-skip cache).
     - **One** small closed registry adapts the upstream contract package-side:
       `attach_sex_routed_aliases()` (one alias clock per `_group.meta.json` `routing.sex` stem).
       It runs inside the build so everything downstream sees ordinary catalog entries. Add to the
       registry; do not add a code path. (A second registry, `CUSTOM_GROUPS`, existed for MiAge's
       undeclared parameter blob; upstream now declares those tensors as ordinary `components`
       plus `code_deps`, so it was deleted.)
     - Verify a sync change by **dry-running the build in memory first** (build catalog +
       bundles, diff every panel against the committed `R/sysdata.rda`) before regenerating.
       `assert_declared_n_cpgs()` is the standing guard: every clock's derived scoring panel must
       equal its declared `n_cpgs`, with no exemption list.
  3. **External packs** (SystemsAge, PCClocks, PCBrainAge, Zhang2019): reuse when `force = FALSE`
     and `data-raw/assets/lockfile.rds` hits (every external clock's `bundle_hash` unchanged and
     every staged pack on disk); else rebuild the content-addressed `<group>-<payload_hash>.qs2`
     packs and rewrite the lockfile. `bundle_hash` (from `manifest.json`) moves iff that clock's
     meta or one of its declared artifacts moved -- unlike `source_git_sha`, which moved on every
     upstream commit and could not say which clock changed.
  4. `upload = TRUE` publishes packs to GitHub Releases; idempotent (content-address + remote
     "asset already present" skip mean unchanged weights are never re-uploaded).
- **Distribution tiers, and the unit is the clock, not the group.** Small groups ship **bundled** in
  `R/sysdata.rda`; the heavy packs ship **external** as release assets, cached at runtime in
  `tools::R_user_dir("methylCIPHERv2", "cache")`. No silent first-use download. `external_group` is
  a **per-clock** field, so a group may be on both sides of the split: `Zhang2019` bundles its
  514-CpG EN arm and packs its 319607-CpG BLUP arm. `EXTERNAL_GROUPS` and `EXTERNAL_CLOCKS`
  (`data-raw/sync.R`) both feed that one field, and everything downstream reads the field --
  `split_group_ids()` puts a mixed group in both buckets, `build_group_bundles(external =)`
  partitions its tensors by declaring member (DECISIONS 2026-07-29).
- **An external clock is not necessarily a pack-scored one.** `clock_is_external()` says where the
  weights live; `is_pack_scored()` says whether `score_pack_group()` computes the score. They agreed
  only while every external clock happened to be a batched weighted sum. Keep them apart: **needing a
  pack** (`pack_groups_needed()`, parity's `skip_if_no_pack()`) keys on externality, and **how the
  arithmetic runs** (`score_type()`'s group hooks, parity's relaxed `packs` tolerance) keys on the
  scoring path. `score_type()`'s external branch carries a group hook per group whose arithmetic is
  not a plain weighted sum -- `SystemsAge` (`center_scale`), `Zhang2019` (`sample_scale`) -- and
  those clocks reach their coefficients through `clock_coefs(id, packs)`, which reads the pack's
  raw tensor by the same `coef_path` the bundled arm reads out of `mc_bundles`. Do **not** hoist the
  `switch(gid, ...)` above the external check to avoid the hooks: that changes dispatch precedence
  for every clock in the catalog to solve a one-group problem (DECISIONS 2026-07-29).
- **"Assets" are the packs; the dir holding them is a cache only in the CRAN sense.** Public names
  say **assets** and nothing else, and **every one of them is `<verb>_mc_<noun>`** --
  `get_mc_assets_dir()` / `set_mc_assets_dir()` (the setter `NULL`-clears and returns the old value
  invisibly), `list_mc_assets()` (read-only table), `download_mc_assets()` (bytes -> disk),
  `load_mc_assets()` (-> RAM), `clear_mc_assets()` (delete). That is the whole pathing surface;
  there is no bare-noun accessor. **Every read-only question has a read-only answer** --
  `list_mc_assets()` reports size / `downloaded` / `superseded` per group without prompting,
  fetching or deleting, so no one has to call a mutating verb to find out what is on disk.
  The word "cache" is reserved for the unrelated internal `partial_cache` (cohort-mean fill). The
  dir itself must stay `tools::R_user_dir(..., "cache")` -- derived, reclaimable, never
  `which = "data"` (DECISIONS 2026-07-24).
- **Assets move in both directions under one consent rule.** `load_mc_assets()` /
  `download_mc_assets()` fill the dir and `clear_mc_assets()` empties it; all three take `ask`,
  prompt interactively, **refuse** non-interactively, and treat `ask = FALSE` as the explicit
  consent signal. Nothing is fetched or deleted unprompted -- CRAN requires a supported way to
  reclaim `R_user_dir()`, so `clear_mc_assets()` must stay a real delete, not a report.
  **Clear means clear** (`pak::cache_clean()` semantics): it removes the currently declared packs
  **and** every superseded one, with no opt-in flag. Filenames are content-addressed, so each sync
  that moves a `payload_hash` orphans the old file; leaving those behind made `clear` fail to
  reclaim and grew the dir without bound. The consent gate is what makes this safe -- the prompt
  counts the two kinds apart ("3 downloaded packs and 3 superseded packs") and lists every file
  before anything is deleted. The stale scan is not a search for a payload -- it never returns one,
  and the stem comes from the declared `file` field, so only the hash is a wildcard; a foreign stem
  or an uncontent-addressed file in the dir is never touched.
  **The gate argument fails closed.** `ask` is a strict flag: only `FALSE` consents, and anything
  that is not a single non-NA logical is an error, never permission. `ext_data` reaching
  `mc_resolve_assets_dir()` is a path or `NULL` only -- a loaded pack names no directory, so it
  stops rather than falling back to the default. Both were silent widenings of "permission" once
  (DECISIONS 2026-07-23); do not re-introduce an `isTRUE()`-style test on either.
- **One argument for the external data, one noun for the thing.** `ext_data` (on
  `load_mc_assets()`, `calc_clocks()`, `sim_DNAm()`) is `NULL` (open set, may download), a path
  (**closed set**, never downloads, missing is fatal), or loaded pack(s). Resolution order is
  `ext_data` > `mc.assets_dir` option > `MC_ASSETS_DIR` env > `R_user_dir` default.
  `download_mc_assets()` / `clear_mc_assets()` take **no** dir argument -- use the setter
  (DECISIONS 2026-07-24). The argument was `from` until 2026-07-28; it is `ext_data` because a
  four-letter English preposition cannot be grepped (DECISIONS 2026-07-28).
- **Identity key:** `payload_hash` (pack content-address) only -- it sets the pack filename and
  release tag, which is what makes re-upload of unchanged weights a no-op. It stays maintainer-side
  and never reaches a result record. Transfer integrity and bit rot are qs2's own
  `validate_checksum`; there is no second hash and no runtime re-hash of a loaded pack.
- **Gitignored, do not commit:** `data-raw/assets/` and `data-raw/methylCIPHER-meta/`.

## Testing

Three tiers. Pre-alpha and fast-moving, so tests guard **core functionality and observable
output**, not implementation detail (see "Test altitude").

- **Crash smoke (always):** `test-sim-smoke.R` scores every bundled clock in the **callable pool**
  (`resolve_clocks("all")`, not `names(mc_catalog)`) through `sim_DNAm()` + `calc_clocks()` with
  `expect_no_error`. External clocks excluded (pack-only); routing targets are covered as their
  alias's dependencies. The pool is **not** filtered by supported-ness -- building it calls
  `score_type()` on every clock, so an unroutable entry fails the tier here.
  **Its value is not "a clock stopped running"** -- parity proves that over a superset of these
  clocks. It is the only tier that runs `calc_clocks()` in the **default configuration** (both
  coverage gates on, full panels, no meta repo), and the only caller of `sim_DNAm()`. Parity
  scores with `min_clocks_coverage = 0, min_samples_coverage = 0` and is skipped on CRAN, so it
  can never stand in for this tier (DECISIONS 2026-07-24).
- **Value goldens (always, no meta dep):** hand-authored engine/machinery unit tests with goldens
  written in-test, one per scoring path (linear sum/mean, sex-split, imputation offset, bundled
  composites). External-pack scoring is smoke-only here; parity owns those goldens.
- **Cohort-gated parity fixtures** (science gate; only clock-golden source): run against **every
  registry cohort** -- `data-raw/methylCIPHER-meta/fixtures/{cohort}/beta.duckdb` for
  `cohort_EPICv1` and `cohort_450K` -- skipped unless BOTH `MC_PARITY=1` and that cohort is staged
  (`file.exists()`). Upstream ships one `fixtures[]` block per cohort; each (clock, cohort) pair is
  its own test. Run locally via the dev-only `test_parity()` (`R/dev-utils.R`). CRAN skips this
  tier; CI must stage the cohorts + set the flag.
  **Never run this tier unless the user explicitly asks for it.** `test_parity()` -- and any
  invocation that sets `MC_PARITY=1`, including `devtools::test()` / `test_file()` under that env
  var -- is minutes-long and reads the staged duckdb cohorts. It is not part of "run the tests":
  the default `devtools::test()` (parity auto-skipped) is. Verify a change against the always-on
  tiers, say that parity was not run, and let the maintainer ask for it.
  **Two axes, both gated, and only the scale-sensitive one ever varies.** Upstream retired its
  `parity` policy (weights_extraction.md sec 12). Every fixture must clear **`max_abs_diff` AND
  `max_rel_diff`** -- both, not either. They are not redundant: the absolute bound is the only one
  with meaning near zero, the relative bound the only one with meaning at large magnitude. Both use
  **`max`, never `median`** -- a median passes while half the samples are arbitrarily wrong.
  **`PARITY_REL_TOL` is `1e-10` everywhere, with no per-block exception** -- it is scale-free, so
  there is never a units-based reason to move it. Only `PARITY_ABS_TOL` is per-block, and relaxing
  it is a statement about units, never about correctness (DECISIONS 2026-07-25).
  **Four blocks, derived from the catalog.** `parity_block()` sends each clock to `horvath`
  (declared `fixtures[].oracle == "horvath_online"`), `packs` (`clock_is_external()`), `fitage`
  (group `DNAmFitAge`), or `core`; `parity_targets(block)` builds one loop per block over the
  shared `run_parity_target()` body. There is no clock list -- every block reads a declaration, so
  a regenerated fixture retires its block automatically. `PARITY_ABS_TOL` is
  `c(core = 1e-10, fitage = 1e-10, packs = 1e-6, horvath = 1e-10)`. The pack relaxation is
  measured, not guessed: those biomarkers reach ~3.3e6 (`PCB2M`), where `1e-10` absolute is below
  the float floor before any arithmetic, and the worst absolute miss over all 56 pack rows is
  2.05e-8 -- which is 6.3e-15 *relative*, i.e. noise. Every pack row sits within 7.5e-13 relative,
  130x inside the unchanged relative bound.
  **The `horvath` block is skipped, and that is a finding, not a shrug.** The oracle filled every
  completely-absent probe server-side with an unpublished per-probe constant (the staged submission
  has 2567/3208 all-NA rows and zero partial NAs) and BMIQ'd the 21k panel for `DNAmAge` only.
  Scored against the *submitted* matrix, the residual tracks the absent-probe count and nothing
  else: pairs with zero absent probes agree to ~1e-8 relative. These clocks are matmuls, so that
  agreement proves the tensors and the engine are right and the divergence is in the oracle's
  input. **Do not "fix" this with a tolerance** -- the residual spans 4.2e-08 to 2.7e-01, so any
  bound wide enough is vacuous (DECISIONS 2026-07-25).
  **`Horvath1` is the one exception to the absent-probe reading, because it is the one the oracle
  BMIQ'd.** 14 of the 15 declare `scheme = none`; `Horvath1` (= the oracle's `DNAmAge`) declares
  `bmiq`, and parity scores it with `normalize` at its opt-in default of **off**. Measured on
  `cohort_450K`, where it has **zero** absent probes on both panels: raw 7.715 abs / 2.70e-1 rel,
  and `normalize = c(Horvath1 = TRUE)` 0.114 abs / 1.93e-3 rel. So its gap is a normalization gap,
  not a fill gap, and "zero absent probes -> ~1e-8" is a claim about the other 14. What survives
  BMIQ is an EM implementation difference (mean +0.024 yr, sd 0.022, uncorrelated with age), and
  the block stays skipped anyway -- admitting it would need a third tolerance regime, which is the
  maintainer's call and has not been made (DECISIONS 2026-07-29).
  **A fixture is scored on the panel the oracle used, which is not always the scoring panel.** A
  recipe that declares a `sample_scale` op z-scores each sample over **every** probe in the input
  matrix, so feeding it the union of scoring panels moves each sample's mean/sd and the score with
  it (measured: 1.8e1 absolute, 82% relative). `needs_full_panel()` reads that op off the declared
  recipe -- never a clock list, today only the two `Zhang2019` arms -- and the target loads the whole
  array with `cohort_betas_full()` instead (DECISIONS 2026-07-25).
  **The blocks are generated, so a dropped fixture is silence, not a failure.** `parity_targets()`
  loops over `clock_fixtures()`; a fixture upstream drops emits **no test at all** -- not a skip,
  just two fewer passes in a green run. One ungated-by-cohort census test guards the generator:
  every catalog clock declares a fixture for **every** `PARITY_COHORTS` cohort, except the 7
  sex-routed aliases, which declare none because their 14 members carry them (both halves derived,
  never listed). It needs no duckdb and no staged cohort, so `test_parity()` runs it even where
  nothing is staged -- but it **is** behind `MC_PARITY`, so a plain `devtools::test()` does not
  catch a dropped fixture; CI does (DECISIONS 2026-07-26).
  **Standing state: 217 tests / 30 skip / 0 fail** -- census 1/1, core 130/130, fitage 28/28,
  packs 56/56, PhysAge 2/2, horvath 30 skipped. **The runner reports this as `PASS 657`**, because
  testthat counts *expectations*, not `test_that` blocks, and `expect_parity()` carries three
  (an all-finite check plus the abs and rel bounds): 214 x 3 + PhysAge 2 x 6 + census 3 = 657.
  Read a parity run by its **fail and skip** counts -- 0 and 30 -- and check the two against each
  other before concluding anything from the pass number.
  **These counts predate the 2026-07-29 Zhang2019 split and have not been re-measured.** The split
  added one clock, so `core` gains its two cohort targets (`Zhang2019BLUP` skips wherever its pack
  is not staged). Re-measure on the next parity run before trusting the totals; the fail count is
  still 0-or-bust.
  `KNOWN_PARITY_GAPS` (clock- or `clock@cohort`-keyed) holds only genuine skips and is **empty**
  today. `KNOWN_PARITY_GAP_GROUPS` (group-keyed) is empty too but stays a **separate** map,
  because group ids and clock ids share a namespace (`DNAmFitAge` is both) and one flat map could
  not say which a key meant.

### Test altitude -- keep tests loose enough to move fast

Assert what `calc_clocks()` *produces*, not how it is wired. A test that breaks on a no-behavior
refactor is too tight -- loosen or delete it.

- **Only `R/` is under test.** `tests/` covers package code and the data it ships, never
  `data-raw/`. `sync.R` is maintainer-side tooling against an upstream contract that upstream
  gates in its own suite; a downstream test of it duplicates that gate, and reaching it means
  sourcing a file the package does not ship and does not depend on. Do not source, parse or
  otherwise bind anything from `data-raw/` in a test.
- **Errors: assert *that*, not the wording.** `expect_error(expr)` with no regex. Pin a message or
  condition class only when a test must otherwise confuse two distinct failure modes.
- **Never `expect_identical()`. Always `expect_equal()`.** `expect_identical()` compares with
  `identical()`, which is bit-exact on doubles and also fails on differences that carry no meaning
  here -- integer vs double storage, an attribute that got dropped or added, a name reordering.
  Floating-point results that are correct to every digit anyone can act on still fail it, and the
  failure looks like a real numeric regression, so time gets spent chasing a last-bit difference in
  a summation order. `expect_equal()` applies a tolerance by default and is the right altitude for
  everything in this package, including counts.
- **No internal dispatch-tag tables.** Do not hard-code `clock_reduction()` / `score_type()` per
  clock; prove routing through output. The one allowed invariant: every catalog clock maps to a
  *known* tag -- and since `score_type()` stops otherwise, that test also proves the catalog
  routes.
- **No maintainer-side plumbing shapes.** Do not assert asset filenames, release tags, download
  URLs, or cache-dir order -- none reach a result. Test behavior (verifies on fetch, leaves no
  scratch, warns-not-stops on hash drift, closed set never downloads).
- **Re-derive a recipe in-test only until parity covers it.** Once a clock has a passing parity
  fixture, that fixture owns the numeric golden and only a smoke stays. In-test re-derivation is
  allowed where parity is still skip-listed -- the only numeric gate meanwhile.
- **Coverage counts and provenance flags are output** -- asserting
  `res$coverage$per_clock[["Hannum"]]$score_imputed_full` or `res$provenance$dependencies` is
  fair game.
- **Minimize test-helper files.** A fixture builder/mock lives atop the one test file that uses it;
  promote to `helper-*.R` only when >= 2 files genuinely share it (currently none). `sim_DNAm` /
  `random_betas` are package functions in `R/`, not test helpers.
- **Cohort/duckdb parity lives in one file** (`test-fixtures-parity.R`): one file-scoped read-only
  connection **per staged cohort** behind the `MC_PARITY` + `file.exists()` guard, torn down with
  `withr::defer(..., testthat::teardown_env())` -- not a module-global caching env.
- **Random inputs are unseeded.** Build DNAm with `random_betas()` (no seed); goldens are computed
  in-test from that same matrix, so they are seed-invariant. Derive the golden from the input, do
  not add a seed to pin a value.

## ASCII-only

Write **plain ASCII** in every file you create or edit -- no "smart" punctuation or symbols.
Use `--`, `->`, `<=` / `>=`, `x` (not em-dash, arrow, inequality/multiplication glyphs), and spell
out set notation.

- **Hard requirement** in package sources (`R/`, `man/`, `DESCRIPTION`, `NAMESPACE`, `tests/`,
  `data-raw/*.R`): non-ASCII triggers R CMD check warnings and breaks on Windows encodings.
- **Default everywhere else** (markdown, commit messages) too, for portability. Some old `dev/*.md`
  lines predate this rule -- do not add more, and prefer ASCII when editing them.

## CLI messages

`cli` is **front-door only**. Keep it for the public interactive surface; everything else is
plain `stop()` / `warning()` / `message()` with `call. = FALSE`.

**Keep `cli` in:**
- assets lifecycle (`R/mc_data.R`: consent, download, clear, path/`ext_data` validation)
- discovery printers (`list_tags`, `print.mc_citation`, `list_clocks` unknown-group)
- public S3 refusals (`rbind.mc_result`, `cite_clocks.default`)
- `calc_clocks` front door: `resolve_clocks` token errors (incl. did-you-mean), DNAm/pheno
  structure (`validate_inputs.R`), coverage gates, value gates / dead samples
  (`missingness.R`), missing pheno in `mc_cohort`, `sim_DNAm` unresolved panels

**Plain `stop()` everywhere else** -- accessors, score branches, pack dispatch, catalog/sync
bugs, normalize-arg validation, citation internals, soft-dep hints (`require_betanorm`), etc.

Rules that still apply on the keep set:

- **Bind every `{?}` plural marker with an explicit `cli::qty()` unless the quantity is the
  interpolation immediately before it.** cli resolves a marker against the *last interpolated
  value earlier in the same string*; with none it scans forward and needs exactly one candidate,
  and the quantity never carries across elements of a `c()` message vector. Get this wrong and
  the handler itself throws (`Cannot pluralize without a quantity`, `Multiple quantities for
  pluralization`) **in place of** the real diagnostic. Safe form:
  `"Add {cli::qty(need)}{?it/them} to {.arg pheno}."`
- The silent variant is worse than the crash: in
  `"{.val {id}} needs pheno column{?s} {.field {need}}"` the marker binds to `id`, so it is
  always singular no matter how many columns are missing. A marker that follows a styled
  `{.val {x}}` is bound to `x`, not to the vector you meant.
- **cli reflows whitespace.** A pre-aligned block (a manifest, a table) collapses onto one line
  when interpolated into a bullet. Use `cli::cli_verbatim()`, which emits lines as-is
  (`mc_manifest_lines()`); anywhere reflow is unavoidable -- inside `cli_abort()` / `cli_inform()`
  bullets -- carry no alignment at all and emit one self-contained bullet per row
  (`mc_manifest_bullets()`).
- **Never hand `askYesNo()` a multi-line prompt.** It passes the string straight to `readline()`,
  whose `prompt` is meant to be one short line; embedded newlines render malformed on Windows.
  Print the context with cli first, then ask a single-line question -- `mc_ask_yes_no()`
  (`R/mc_data.R`) is the one place that does this and every consent prompt goes through it. It
  takes its `header` **pre-formatted** via `cli::format_inline()` in the caller's frame and
  interpolates it as a value, so pluralization resolves against the caller's variables and the
  text is never re-parsed for braces.
- Tests assert *that* a message errors, never its wording -- see "Test altitude".

## Comments

- Plain `#` comments are the only in-source docs right now. Roxygen is enabled but reserved for
  machinery -- `@export`, `@useDynLib`, `@importFrom` -- not for prose (see invariants).
- Keep them **short**: 1-2 sentences on *what* the code does, not a rationale essay.
- The *why*, and every decision or reversal, goes only in `dev/DECISIONS.md`.

## Source-of-truth docs (`dev/`)

The `dev/` folder is local-only **except** these two, which are tracked:

- `dev/DECISIONS.md` -- append-only, newest-first, date-stamped log of *why* / reversals. Add an
  entry when a decision reverses a prior approach or is likely second-guessed; do not restate rules
  already stated here.
- `dev/id-streaming-plan.md` -- the one live design doc: chunking, binding, `prep()`. It covers
  work that is **not built yet**; everything already shipped is specified by this file's invariants
  and by the code.

`migration-plan.md` and `detail-plan.md` were retired on 2026-07-28 (DECISIONS). **Do not
reconstitute them.** Built behavior is specified by the invariants above plus the code; a separate
long-form spec of shipped behavior is a copy that rots. `sec N` citations in older DECISIONS
entries point at those retired files -- read them out of git history, not as live references.
What upstream declares (coef-path rule, declared-path set, tensor `row_key`/`col_key`, recipe
operand namespaces, the panel rule) is **not** restated in a `dev/` doc -- `data-raw/sync.R` is
self-documenting and is the only source for it. Read `sync.R` itself before touching `sync.R`.

The plan states **current truth only** -- superseded design is not annotated inline; its history
lives solely in `dev/DECISIONS.md`. When code and the plan disagree, the code is truth: fix the plan
and record the reconciliation in `dev/DECISIONS.md`.

Local-only (gitignored): `dev/legacy/` (frozen pre-rewrite sources), `dev/scratch.R`,
`dev/clock_tracker.csv`, and the `dev/*.py` build scripts.

## Contributing

- Branch off `main` and open a PR; do not push to `main`.
- Run `devtools::test()` before pushing. Run `devtools::document()` when you add or change a roxygen
  tag, and commit the regenerated `NAMESPACE` / `man/` alongside it.
- Reversing or second-guessing a design? Add a dated, newest-first `dev/DECISIONS.md` entry.
- Keep new or edited content ASCII.

## Environment and personal overrides

Keep **this** file environment-agnostic -- it is shared across operating systems and shells.

- The tracked `.Rprofile` auto-attaches `devtools` + `testthat` in interactive sessions. For a
  clean, profile-free parse or check, use `Rscript --vanilla` or `R CMD check`.
- Put machine-specific or personal notes (OS, shell, local paths, private scratch) in
  `CLAUDE.local.md` -- gitignored, loaded automatically, never reaches a collaborator.
