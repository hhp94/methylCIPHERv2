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

- **One beta entry point, and therefore no pre-flight check.** `calc_clocks()` is the only public
  surface that reads a beta matrix. Everything else reads the **catalog** (`list_clocks`,
  `clock_cpgs`, `list_clock_tags`) or a **finished record** (`clocks_coverage`, `samples_coverage`,
  `calc_accel`, `score_associations`, `refinalize_clocks`, `cite_clocks`). Of the exports only
  `calc_clocks` and `predict_sex` take a `DNAm` argument at all, and `predict_sex` touches it
  exclusively *through* `calc_clocks`; `sim_DNAm` generates a matrix rather than reading one.
  So a "dry run", a coverage preview, or a `report(DNAm)` arm is **a second beta reader**, and
  that is the thing being refused: a second reader takes its own independently-supplied matrix, so
  its verdict can be about a different object than the one that gets scored -- decoupling risk with
  nothing bought, since scoring is a matmul over an already-resident matrix and **both coverage
  gates are arguments**, so `min_clocks_coverage = 0, min_samples_coverage = 0` already yields the
  full report with no refusal. The pre-flight habit is inherited from upstream and does not
  transfer: an ENmix/minfi pipeline couples its steps through one shared object and caches an
  expensive IDAT parse, so checking before computing is both safe and necessary there. We have no
  pipeline object and no expensive parse. `predict_sex()` is **not** an exception -- it is a
  `calc_clocks()` call whose output feeds a later call, which is composition, not a pre-check.
  This rests on scoring staying cheap enough that running it is not a commitment; a streaming or
  chunked path would weaken that premise and would need an explicit answer rather than an
  inherited one (DECISIONS 2026-08-03).
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
  (`R/coverage.R`), keyed by clock id, and merged in `construct_mc_result()` -- which files it
  under the run's **batch label**, so the record reads `$coverage$per_clock[[batch]][[id]]` even
  for a single-pass run. Batch is a real axis for these counts and not a decoration: `rbind` binds
  records scored in separate `calc_clocks()` calls, and `score_imputed_partial` counts the panel
  CpGs in *that run's* partial cache, so two batches almost never agree (DECISIONS 2026-07-30,
  "Phase 4 gates"). Per-sample miss is
  counted **once per distinct panel** (FitAge/GrimAge reuse panels) and kept **per panel role**:
  every clock has a score panel; a normalizing clock (only DunedinPACE) also has a norm panel. The
  counts follow the policy uniformly -- `vendor_mean` fills every absent CpG into the predictor
  (`used = present + imputed_full`), anything else drops them (`used = present`, `dropped = absent`).
  **The norm panel fills or drops on the declared scheme, not on that policy**: quantile
  normalization needs the whole background panel and takes the target's value for an absent CpG,
  bmiq calibrates on what is present, so `norm_imputed_full` / `norm_dropped` key on
  `NORM_SCHEMES_FILL` (DECISIONS 2026-07-29). There is no `norm_used` -- it would be `norm_needed`
  or `norm_present`, both already reported.
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
  `$provenance` also carries the per-sample `mc_batch_id` (aligned to `sample_id`) and the retained
  `pending` intermediates that make an opt-in `refinalize_clocks()` exact.
  It carries **both coverage floors, keyed by batch label** like `$coverage$per_clock` --
  `rbind` keeps one per batch and reconciles nothing, and `samples_coverage()` finalizes them by
  taking the **most restrictive** (`max`) and re-warning under it, which is the only thing that
  makes a post-bind `coverage < threshold` filter well defined. There is no `below_min` column:
  the cell axis already exists, and a conditional one would mean different things per row after a
  bind. `min_clocks_coverage` is recorded but read by nothing -- it aborts, so a record's existence
  already proves it passed (DECISIONS 2026-08-03).
  **Where a verb exists it is a method**, and the built surface today is exactly `print`,
  `as.matrix`, `as.data.frame`, `cite_clocks` and `rbind`, plus the plain `calc_accel()`.
  Coverage is
  deliberately not a method: it is the plain `clocks_coverage()` / `samples_coverage()`, **not**
  `summary()`. `clocks_coverage()` is one row per **(clock, batch)**; `samples_coverage()` carries
  each sample's batch alongside its id. In both, `mc_batch_id` is the **last** column -- it is the
  key the two frames join on, but it is a hash, so it does not sit in front of `clock_id`.
  **And it reaches an exit frame only when the record spans more than one batch.** At one batch it
  is a single repeated hash: it carries no information, `clock_id` alone is already unique so the
  join still resolves, and prose docs are deferred so a user has nowhere to look up what it means.
  All four exits (`as.data.frame`, `calc_accel`, both coverage frames) share the **one** test in
  `is_multi_batch()` (`R/mc_result.R`), keyed on `length(unique(provenance[[mc_batch_id]]))` --
  the vector that fills the column, not `per_clock`'s names, so a frame can never be keyed on a
  count other than its own contents. **Every exit reads that test to decline building the column**
  rather than building one it will lose; `drop_single_batch()` still runs at all four and is the
  gate, now a no-op. The two forms produce `identical()` frames. The coverage frames thread the
  decision down to `panel_rows()` / `batch_coverage()` because they assemble by `rbind` and cannot
  add the column after the fact (DECISIONS 2026-08-01). A conditional schema is a real cost and is accepted
  knowingly: the four exits must appear and vanish **together** or the two coverage frames disagree
  about whether the join key exists (DECISIONS 2026-07-31, "the batch label is multi-batch only").
  Nothing internal is conditional -- `$provenance` always carries the per-sample vector, and
  `calc_accel()` always puts `mc_batch_id` in the formula namespace and always reserves it
  against `data =`.
  **The batch label is `mc_batch_id` everywhere the user can touch it** -- both coverage frames,
  both finalizer frames, `$provenance`, and the `calc_accel()` formula namespace. One string,
  and `mc_`-prefixed on purpose: a user's `batch` is their slides or plates, which is biology and
  is a covariate they may legitimately want in a model, while ours only says which samples shared
  a cohort-mean fill. A bare `batch` would both shadow theirs in a formula and collide on a join
  against their own metadata. `data =` supplying `mc_batch_id` is an error, not a precedence
  question (DECISIONS 2026-07-31, "one name for the batch label").
  **`rbind` binds and labels; it never reconciles.** Nothing is re-imputed, no denominator is
  merged, and no cross-sample column is recomputed unless `refinalize_clocks()` is called by hand
  -- but a multi-batch bind carrying a non-empty `pending` **says so once** (`say_pending()`,
  an `inform` not a `warn`), naming the columns from `names(pending)`, which is the catalog's
  declared `cross_sample` set. `refinalize_clocks()` **reads `pending` and never consumes it**, so
  it composes in any order and a second call is a no-op; never clear `pending` after re-finalizing.
  **`rbind` is the only verb that leaves `pending` unresolved. Every finalizer resolves it.**
  **A finalizer is any exit that takes an `mc_result` and returns something that is not one.** The
  test is mechanical, so the set is derived rather than listed and cannot go stale: it is
  `as.data.frame()`, `as.matrix()`, `calc_accel()` and `score_associations()` today. All four
  re-finalize on the way out and say so, under `say_pending()`'s exact guard -- non-empty
  `pending` **and** more than one batch. The set was enumerated instead of derived until
  2026-08-03 and drifted twice under it: `as.matrix()` was outside it until then, which meant it
  and `as.data.frame()` could hand back different numbers for the same multi-batch record, and
  `score_associations()` was re-finalizing without the enumeration or its own docs saying so
  (DECISIONS 2026-08-03). `rbind` is **not** a finalizer under the same test -- it returns an
  `mc_result` -- and must not become one, because `do.call(rbind, ...)` recurses and would
  re-finalize at every intermediate step; a finalizer is a leaf and hands back a value the record
  cannot be recovered from, so it must hand back the right numbers.
  **The two batch counts are cross-checked, not chosen between.** `n_batches()` (`R/mc_result.R`)
  derives the count from `provenance[[mc_batch_id]]`, the per-sample vector that fills the column,
  and `stop()`s if `length(per_clock)` disagrees. Every finalizer and both coverage frames route
  through it, and `finalized()` calls it **before** the `pending` test so `&&` cannot short-circuit
  past the check on the common path. Do not read either count directly for a batch total. A single-batch record is skipped because its reduction already spans its
  whole cohort -- re-finalizing is a numerical no-op there and the message would be pure noise
  (DECISIONS 2026-07-31, "one name for the batch label"). Wanting the per-batch reductions is
  still possible: finalize each record *before* binding.
  Its gates follow one line -- **record what batching forced, refuse what the caller chose
  differently** -- so a per-batch fill regime is recorded (that is what the batch axis is for)
  while overlapping ids, differing score columns, a differing `pheno_id` and a differing
  `normalize=` all throw. There is no `force =`. **Batch labels are derived, never assigned**:
  `construct_mc_result()` sets one to `batch_hash(pheno[[pheno_id]])` -- the full 16 hex of
  `xxhash64` over the pheno's **id column only**, never truncated -- and there is no `batch =`
  argument anywhere. So `rbind` carries
  no label policy at all: it mints nothing, renames nothing, renumbers nothing, and **drops**
  argument names (`unname(list(...))`). Do not refuse them instead: `split()` names its result, so
  refusing kills `do.call(rbind, lapply(split(...), ...))`, which is the blocking idiom the feature
  exists for. This is what makes re-association exact --
  `rbind(rbind(r1, r2), r3)` and `rbind(r1, r2, r3)` return `identical()` records. Do not hash the
  pheno *frame*: covariate values would fold into batch identity (correcting one subject's age
  renames the batch) and `digest` is storage-type sensitive. Do not restore an assigned label to
  fix a spelling -- a stored string cannot say whether a human chose it, which is exactly what sank
  `is_auto_label()` (DECISIONS 2026-07-30, "The batch label is derived"). **Hash the canonical
  form** -- `sort(..., method = "radix")` + `unname`/`as.character` + `serialize = FALSE` -- or the
  label follows the id *sequence*, the locale and the R serialization version instead of the id
  set. Two batches sharing a label then needs a 64-bit collision against id sets gate 1 just made
  disjoint (~2.7e-16 at 100 batches), so there is **no gate on labels** -- do not add one back.
  (An older DECISIONS entry says 12 hex / 48 bits and 1.8e-11; that was always wrong about the
  shipped width -- `batch_hash()` does not truncate -- and the conclusion only gets stronger.)
  `sim_DNAm(suffix =)` is **not** a batch
  argument -- it suffixes sample ids so two simulated blocks clear gate 1.
  Citations dispatch as `cite_clocks()` -- a **package-owned** generic, because both
  `utils::citation` and `utils::cite` already exist as plain functions and taking either name
  masks it (DECISIONS 2026-07-23, 2026-07-24, 2026-07-25). Nothing else is promised: `as.data.frame`,
  `[`, `cbind`, `augment` and `codebook` were listed here for a year without being written, so
  they are **unbuilt ideas, not contracts** -- adding one is a new API decision, and until a human
  makes it there is no behaviour to match (DECISIONS 2026-07-27).
  **Every `print.mc_*` method shares one grammar, built in `R/print.R`** -- a `<class> A x B` header,
  a `$component [what is shown]` line per list element, then `... N more <axis>`. The builders return
  strings, so a cli printer (`print.mc_citation`) and a `cat` printer emit identical text without
  moving the cli boundary. A new record class reuses the builders; it does not invent a fourth
  layout (DECISIONS 2026-07-29).
- **Scores only, and the record remembers its inputs.** `$scores` is scores -- no auto-appended
  phenotype columns. Separately, `$pheno` carries the *aligned* pheno narrowed to the id column
  plus the covariates the run actually required, so a saved record can answer what was fed in.
  **It is never `NULL`**: with no `pheno =` supplied, `resolve_pheno()` materializes the id column
  alone, so `$pheno` is one shape everywhere and the batch label always has a column to hash
  (DECISIONS 2026-07-30, "The batch label is derived"). Its columns therefore need no `rbind` gate
  -- they are `unique(c(pheno_id, covariates))`, both already pinned.
  **And it never carries row names.** `resolve_pheno()` has one exit and resets them there, on both
  branches, so a supplied pheno's row names never survive the id-join. Row names would be a second
  identity beside the id column, and the two can silently disagree under subsetting or `rbind`; the
  id column is the only key. `rbind` adds no reset of its own -- automatic row names stay automatic
  through `rbind`, so a second one would be guarding what `resolve_pheno()` already guarantees.
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
- **Never resolve a name by partial match. `$` is one instance of a general rule.** A token the
  user supplied is matched **exactly** against its closed set, so `pmatch()`, `charmatch()` and
  `match.arg()`'s abbreviation handling are all out at the front door. Two reasons, and the second
  is specific to this package. A partial match resolves silently to something the caller did not
  name, which is the same failure as `entry$covariates`. And **the closed sets here grow with every
  sync**: `"Sys"` resolves to `SystemsAge` today, and the sync that adds a second `Sys*` group turns
  working user code into an ambiguity error. An abbreviation cannot be supported without freezing
  the catalog.
  The tool is **`checkmate::assert_subset()`** -- exact, names the offending element **and** the
  valid set, and takes `empty.ok` as an explicit flag. Precedent: `list_clocks(tag =)` and
  `mc_resolve_groups()`. **`match.arg(several.ok = TRUE)` is not an alternative**, and was measured
  before being rejected (DECISIONS 2026-08-03): it returns `choices[1L]` for `NULL` (here, `"all"`
  -- a silent mass download), it errors only when *every* element fails so a typo beside a valid
  token is **dropped without a word**, and it does not deduplicate. It bans `character(0)`, which is
  the one property it has that we wanted, and that is one line to write ourselves.
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
- **`checkmate` asserts at the exported surface; internals get bare `stop()`.** Internals assume
  validated input -- by the time a helper runs the value has already crossed the front door, so a
  second assertion fires only when *we* have a bug, at which point a `checkmate` message aimed at a
  user is the wrong register. An export called by another export is a **trusted** caller
  (`predict_sex` -> `calc_clocks`), so the inner assertion is redundant, not protective; the
  exception is a value the inner one reads that the outer never validated, which is why
  `recorded_from_female()` keeps its assert. An internal guard that must stay -- a bounds check
  ahead of a kernel, say -- keeps its checks and drops to `stop()` with short greppable text
  (`check_moment_sets()` in `R/missingness.R` is the pattern).
  **`.var.name` is filled only where the deparse lies.** `checkmate` names the failing value by
  deparsing the expression, so at a boundary it already prints the caller's own word and a
  hand-written string is staleness risk. Fill it **iff the deparsed expression does not name
  something the caller can locate in their own call** -- which is not the same as "not a bare
  symbol": `assert_character(colnames(DNAm))` deparses to `colnames(DNAm)` and is left alone, while
  `check_pheno()`'s `pheno[[ID]]` names an internal parameter and carries
  `.var.name = paste0("pheno$", ID)`. **Never as a substitute for moving a check**: a wrong name in
  a `checkmate` message is usually evidence the check sits in the wrong frame -- that is how the
  missing front-door assert on `min_clocks_coverage` surfaced, as `Assertion on 'threshold' failed`
  (DECISIONS 2026-08-03).
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
- **Read `dev/WRITING.md` before editing any text a user can see. Reading it is the first step of
  the task, not a review at the end.** The scope is the whole user-facing surface and not only
  roxygen: doc blocks, `README.Rmd`, `vignettes/*.Rmd`, and every cli message. It is **not** a
  style guide to consult when unsure. It holds load-bearing, non-obvious rules, and close to every
  one of them is there because that exact failure has already happened in this package -- R1 to R8
  (the English rules), the cli mechanics (`sprintf` output must never become a cli template, cap
  before format, `cli::qty()` on every plural marker), the roxygen tag order and the `DOC_TYPES`
  param vocabulary, the shared-parameter donor's inheritance footgun, the **closed** `@seealso`
  set, the vignette-versus-README build split, and the line-break and ASCII rules for prose files.
  An agent that writes first and reads afterwards breaks a rule it did not know existed, and
  `lint_roxygen()` / `lint_seealso()` catch only the mechanical subset. Where a change and the file
  disagree, **update the file in the same pass** -- a rule the shipped files violate is worse than
  no rule.
- **The exported surface is documented, and `dev/WRITING.md` is how.** Roxygen went on because the
  package went Rcpp and `useDynLib` has no route into `NAMESPACE` except a tag, so `NAMESPACE` and
  `man/*.Rd` are **generated files** and `devtools::document()` is a normal part of the workflow.
  Prose docs were deferred behind that for a year and **shipped 2026-08-03**: all 26 user-facing
  topics now carry real `@param` / `@details` / `@returns` / `@examples`, the `@seealso` groups are
  closed, and `cite_clocks` merges its three methods onto one topic with `@rdname`. The tag order,
  the `DOC_TYPES` param vocabulary, the shared-parameter donor and its footgun, and the closed
  cross-reference set all live in `dev/WRITING.md`, under the invariant above.
  `lint_roxygen()` and `lint_seealso()` (`R/dev-utils.R`) must both come back
  empty; both were empty on 2026-08-03 (DECISIONS 2026-07-27, 2026-08-03).
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
  dir itself stays `tools::R_user_dir(..., "cache")` -- but **not because `which = "data"` would
  breach CRAN policy. It would not.** Policy permits data, configuration and cache alike, and its
  two conditions ("sizes as small as possible", "actively managed") attach to all three equally, so
  neither discriminates and `clear_mc_assets()` is required either way. This is a **platform call**:
  `"data"` is `%APPDATA%` on Windows, which roams -- a LAN copy at every logon under a roaming
  profile, or SMB reads on every scoring run under Folder Redirection -- while `"cache"` is
  `%LOCALAPPDATA%` and never roams. The macOS counter-risk is a purge of `~/Library/Caches`, whose
  cost is one consented re-download of a 43 MB corpus, and which the content-addressed filenames
  already produce routinely whenever a `payload_hash` moves. Precedent agrees: `BiocFileCache` and
  `ExperimentHub` -- the stack these users already run -- both resolve under `R_user_dir(..., "cache")`.
  **Durability is the job of `MC_ASSETS_DIR` / `set_mc_assets_dir()`, not of the default**
  (DECISIONS 2026-08-04, superseding the "policy violation" reading of DECISIONS 2026-07-24).
- **Assets move in both directions under one consent rule.** `load_mc_assets()` /
  `download_mc_assets()` fill the dir and `clear_mc_assets()` empties it; all three take `ask`,
  prompt interactively, **refuse** non-interactively, and treat `ask = FALSE` as the explicit
  consent signal. Nothing is fetched or deleted unprompted -- CRAN requires a supported way to
  reclaim `R_user_dir()`, so `clear_mc_assets()` must stay a real delete, not a report.
  **Clear means clear** (`pak::cache_clean()` semantics): it removes the currently declared packs
  **and** every superseded one, with no opt-in flag. Filenames are content-addressed, so each sync
  that moves a `payload_hash` orphans the old file; leaving those behind made `clear` fail to
  reclaim and grew the dir without bound. The consent gate is what makes this safe -- the prompt
  counts the two kinds apart ("3 downloaded assets and 3 superseded assets") and lists the files
  before anything is deleted -- capped at `MC_MSG_CAP` like every other list, with any remainder
  counted on its own line, because a message that renders an unbounded vector is the one failure
  mode the cap exists for and the counts above it are already exact (DECISIONS 2026-08-03). The
  stale scan is not a search for a payload -- it never returns one,
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
  **One golden in this tier is a third-party reference implementation, not a re-derivation**: the
  `DunedinPACE` package (Suggests, `danbelsky/DunedinPACE`) scores the degraded-coverage path that
  parity's clean fixtures cannot reach. It skips when the reference, `betanorm`, or the reference's
  own `preprocessCore` is absent, so it never becomes a hard dep (DECISIONS 2026-07-29).
- **Cohort-gated parity fixtures** (science gate; the only clock-golden source for a clean panel):
  run against **every registry cohort** --
  `data-raw/methylCIPHER-meta/fixtures/{cohort}/beta.duckdb` for
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
  **Standing state, re-measured 2026-08-02: 263 blocks / 32 skip / 0 fail** -- census 1/1,
  core 146, fitage 28, packs 56, PhysAge 2/2, horvath 30 skipped, Wang@cohort_450K 2 skipped.
  **The runner reports this as `PASS 699`**, because testthat counts *expectations*, not
  `test_that` blocks, and `expect_parity()` carries three (an all-finite check plus the abs and rel
  bounds): 228 targets run x 3 + PhysAge 2 x 6 + census 3 = 699.
  Read a parity run by its **fail and skip** counts -- 0 and 32 -- and check the two against each
  other before concluding anything from the pass number.
  `core` grew 130 -> 146 across the Zhang2019 split and the `DNAmSex_Wang` family; six of those
  are accounted for and the other ten arrived with an uncommitted `sysdata.rda` regeneration
  (DECISIONS 2026-08-02, "Wang parity"). The fail count is 0-or-bust either way.
  `KNOWN_PARITY_GAPS` (clock- or `clock@cohort`-keyed) holds only genuine skips. It holds
  **two** today, both `DNAmSex_Wang_*@cohort_450K`: that cohort's deposited matrix carries no
  sex-chromosome probes, so the panel is 0% present and the fixture is the oracle's empty-panel
  `0`. **Do not relax `check_coverage()`'s `ratio == 0` stop to make them pass** -- a 0 there is
  the `Female` quadrant of the sign map, not a small number.
  `KNOWN_PARITY_GAP_GROUPS` (group-keyed) is empty but stays a **separate** map,
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

**The line is audience, not transport.** A message about **input the user chose** is `cli`,
whatever function raises it. A message about a **package defect** -- a "cannot happen" condition,
a catalog/sync gap, a missing dispatch branch -- is a plain `stop()` with `call. = FALSE`. This
replaced an enumerated cli keep-set on 2026-08-03, because that list put four `resolve_normalize()`
messages about a `calc_clocks()` argument on the `stop()` side, where they could not carry markup
and so **could not meet the rules below at all**. A rule a whole class of user-facing messages is
structurally incapable of meeting is a two-tier system, not a rule (DECISIONS 2026-08-03).

So today: assets lifecycle, discovery printers, public S3 refusals, and the whole `calc_clocks`
front door (token resolution, `validate_inputs.R`, coverage gates, `missingness.R`, `mc_cohort`,
`clock_cpgs`, `resolve_normalize`) are cli. Accessors, score branches, pack dispatch, catalog/sync
bugs and citation internals are `stop()`. **`list_clock_tags()` is not a printer** -- it returns
the registry as a value and prints nothing (DECISIONS 2026-07-29).

A defect message is **hard-coded and greppable**: a fixed prefix, values appended after it, so a
bug report can be located from the pasted text with no stack trace asked of the user. Not a fully
constant string -- `check_moment_sets()`'s failing index is real debugging value -- but the fixed
part leads.

### How the text itself is written: `dev/WRITING.md`

**`dev/WRITING.md` is the single source, and this file does not restate it.** Read it before you
write or edit any text a user can see. It is tracked, so it resolves in a fresh clone. It holds:

- **R1 to R8**, the English rules. They bind **every message a user can see** -- cli message text
  and roxygen prose alike. They do **not** bind code comments, dev-facing `stop()` text,
  `data-raw/`, or `dev/` docs, so the ASCII section above is untouched and `--` stays required
  there (DECISIONS 2026-08-03).
- **The cli mechanics**: why `sprintf` output must never become cli input, why `bullets()` escapes
  braces, cap-before-format, `cli::qty()` on every plural marker, `cli_verbatim()` for anything
  pre-aligned, and never a multi-line `askYesNo()` prompt. Each is a bug that has actually
  happened here.
- **The roxygen template**: tag order, the `DOC_TYPES` param vocabulary, the shared-parameter
  donor and its one footgun, the closed `@seealso` groups, and the example rules.
- **`say_*` emits to the user; `note_*` records into the block's collector.** Do not use `note_`
  for something that prints.
- **The audit section**: the known-good exceptions an independent reader will otherwise report as
  defects, and the three CRAN shape rules the manual currently satisfies.

Two rules from it are repeated here only because they are enforced by the test suite and the
workflow rather than by prose: tests assert *that* a message errors and never its wording (see
"Test altitude"), and `lint_roxygen()` plus `lint_seealso()` must both be empty before a doc
change is done.

## Comments

- **Plain `#` comments are for the code; roxygen is for the manual.** Both are live. The exported
  surface carries real roxygen prose as of 2026-08-03, written to `dev/WRITING.md`.
- Keep `#` comments **short**: 1-2 sentences on *what* the code does, not a rationale essay.
- The *why*, and every decision or reversal, goes only in `dev/DECISIONS.md`.

## Source-of-truth docs (`dev/`)

The `dev/` folder is local-only **except** these three, which are tracked:

- `dev/DECISIONS.md` -- append-only, newest-first, date-stamped log of *why* / reversals
  (2026-07-30 and later). Add an entry when a decision reverses a prior approach or is likely
  second-guessed; do not restate rules already stated here.
- `dev/DECISIONS.old.md` -- full pre-2026-07-30 decision history. Dated citations earlier than
  that cut resolve here; do not restate that archive in the live log.
- `dev/WRITING.md` -- the single source for how user-facing text is written: the English rules,
  the cli mechanics, the roxygen template, the closed `@seealso` groups, and the manual's
  known-good exceptions. See "CLI messages" above; this file points there and does not restate it.

**There is no live design doc, and that is deliberate.** `migration-plan.md` and `detail-plan.md`
were retired on 2026-07-28, and `id-streaming-plan.md` -- which held the chunking / binding /
`prep()` design -- was deleted on 2026-08-02 once its shipped half was covered by the invariants
above. **Do not reconstitute any of them.** Built behavior is specified by the invariants plus the
code; a separate long-form spec of shipped behavior is a copy that rots. `sec N` citations in older
DECISIONS entries point at those retired files -- read them out of git history, not as live
references. What upstream declares (coef-path rule, declared-path set, tensor `row_key`/`col_key`,
recipe operand namespaces, the panel rule) is **not** restated in a `dev/` doc -- `data-raw/sync.R`
is self-documenting and is the only source for it. Read `sync.R` itself before touching `sync.R`.

So **the code is truth**, with no plan to reconcile it against. Unbuilt design lives in a
`dev/DECISIONS.md` entry stating the decision, or it is not written down yet.

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
