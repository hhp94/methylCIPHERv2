# CLAUDE.md

Guidance for Claude Code in this repo. This file holds **invariants** that rarely change.
Volatile detail (per-clock status, exact designs, dated reversals) lives in `dev/` -- see
"Source-of-truth docs" and prefer it for specifics. A `(DECISIONS <date>)` tag points at the entry
holding the measurements and the full argument: this file states the rule and the shortest reason,
and does not restate the evidence.

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
Soft deps back specific paths only and skip when absent. **A dep is declared for code, not for
tests**: `betanorm` is in `Suggests` because `R/` reads it (guarded by `require_betanorm()` in
`R/score_normalized.R`), while `duckdb`, `DBI` and `DunedinPACE` are **not declared at all** -- all
three are read only by `test-fixtures-parity.R`, which is `.Rbuildignore`d and so never reaches a
tarball for the unstated-dependency scan to see. `DBI` was dropped from `Suggests` on 2026-08-04
once that was true; it was never read by `R/`, contrary to what this file used to say.
`withr` stays in `Suggests` for the same static-scan reason but is effectively free -- testthat
`Imports` it, so anything that can run the suite already has it, and it needs no
`skip_if_not_installed()` (DECISIONS 2026-08-04).

## Non-negotiable invariants

Do not reverse these without a `dev/DECISIONS.md` entry explaining why.

- **Never run `R CMD check` / `devtools::check()`. It is on-demand, maintainer-only.** It hangs
  here: check re-runs the full suite in a fresh install and does not finish in usable time, so an
  agent that starts one blocks the session on something it cannot resolve. **Verify a change with
  `devtools::test()`** and say plainly that check was not run. Same rule and reason as the parity
  tier: minutes-to-hours of maintainer wall-clock is the maintainer's call. The prohibition is on
  the work, not one entry point -- not via `Rscript`, `pkgbuild`, or a background shell either.
- **One beta entry point, and therefore no pre-flight check.** `calc_clocks()` is the only public
  surface that reads a beta matrix; everything else reads the **catalog** (`list_clocks`,
  `clock_cpgs`, `list_clock_tags`) or a **finished record** (`clocks_coverage`, `samples_coverage`,
  `calc_accel`, `score_associations`, `refinalize_clocks`, `cite_clocks`). Only `calc_clocks` and
  `predict_sex` take a `DNAm` argument, and `predict_sex` touches it exclusively *through*
  `calc_clocks`; `sim_DNAm` generates a matrix rather than reading one. A "dry run", a coverage
  preview, or a `report(DNAm)` arm is **a second beta reader**, and that is what is refused: it
  takes its own independently-supplied matrix, so its verdict can be about a different object than
  the one that gets scored, and it buys nothing -- scoring is a matmul over an already-resident
  matrix and **both coverage gates are arguments**, so `min_clocks_coverage = 0,
  min_samples_coverage = 0` already yields the full report with no refusal. The pre-flight habit is
  inherited from ENmix/minfi, which couple their steps through one shared object and cache an
  expensive IDAT parse; we have neither. `predict_sex()` is composition, not a pre-check. This rests
  on scoring staying cheap; a streaming or chunked path would need an explicit answer rather than an
  inherited one (DECISIONS 2026-08-03).
- **One engine + a finite, closed branch set.** Every work unit routes on the catalog pair
  `(weights_format, computation_type)` to shared `linear_score()` or a named branch (pre-transform,
  family orchestrator, sex-routed alias, external, custom). There is **no** recipe
  interpreter/walker, not even as a fallback.
- **Routing is total, and a gap is a hard stop.** `score_type()` returns a known tag for every
  catalog entry or `stop()`s naming the clock's group / `weights_format` / `computation_type`. There
  is no `"unsupported"` tag and nothing filters on one: a sync that adds a routing pair no branch
  claims must fail the always-on tier, never silently shrink it (DECISIONS 2026-07-23).
- **A branch returns only its score** (`<n x 1 matrix>`), never a coverage record. One declared
  exception: a clock the catalog marks cohort-reducing (`cross_sample_at`, today the two
  `DNAmPhysAge` clocks) returns its **per-sample intermediate**, which the scoring loop routes into
  `pending` rather than `scores`. `finalize_cross_sample()` reduces once after assembly and every
  front end calls it unconditionally (a no-op when `pending` is empty); the routing reads
  `spec$cross_sample`, never a clock id and never whether the run was chunked (DECISIONS 2026-07-27).
  - Coverage depends on no score, so `compute_coverage()` (`R/coverage.R`) computes it once upstream
    of the scoring loop, keyed by clock id, and `construct_mc_result()` merges it under the run's
    **batch label**: `$coverage$per_clock[[batch]][[id]]`, even for a single-pass run. Batch is a
    real axis -- `score_imputed_partial` counts *that run's* partial cache, so two bound batches
    almost never agree (DECISIONS 2026-07-30).
  - Per-sample miss is counted **once per distinct panel** (FitAge/GrimAge reuse panels) and kept
    **per panel role**: every clock has a score panel, a normalizing clock (only DunedinPACE) also
    has a norm panel. The score panel follows the policy uniformly -- `vendor_mean` fills every
    absent CpG (`used = present + imputed_full`), anything else drops them (`used = present`,
    `dropped = absent`). **The norm panel fills or drops on the declared scheme, not on that
    policy** (quantile needs the whole background panel, bmiq calibrates on what is present), so
    `norm_imputed_full` / `norm_dropped` key on `NORM_SCHEMES_FILL`. There is no `norm_used`.
    Routed members are masked to the samples they scored.
  - `coverage_record()` owns the per-clock fields, including `normalizes` (**the one declared panel
    fact**; never re-derive it from `norm_needed`) and the per-panel `score_imputed_partial` /
    `norm_imputed_partial`. **The record is one axis: every count in it is a CpG count**, the
    partials included and block-invariant. Never store `sum(sample_miss)` there: it collapses both
    axes and can exceed `score_present`.
  - **The sample axis is `$coverage$sample_miss`, and only that**:
    `list(score = <n x k>, norm = <n x k' over just the normalizing columns>)`. **`$per_clock`,
    `$sample_miss` and `samples_coverage()` span one set -- the clocks that read CpGs** -- so
    routing targets are in and pure composites are out even when returned. `cached_cols()` /
    `count_sample_miss()` (integer) / `score_matrix()` are the shared shape helpers (DECISIONS
    2026-07-24, 2026-07-29).
- **Result is an S3 record over `list`** (class `mc_result`): `$scores` (n x k double), `$pheno`,
  `$coverage`, `$provenance`. Never a `matrix` subclass (drops class + attrs on first subset).
  `$provenance` carries the per-sample `mc_batch_id` (aligned to `sample_id`), the retained
  `pending` intermediates that make an opt-in `refinalize_clocks()` exact, and **both coverage
  floors keyed by batch label** like `$coverage$per_clock`. `rbind` keeps one floor per batch and
  reconciles nothing; `samples_coverage()` finalizes them by taking the **most restrictive** (`max`)
  and re-warning, which is the only thing that makes a post-bind `coverage < threshold` filter well
  defined. There is no `below_min` column -- it would mean different things per row after a bind.
  `min_clocks_coverage` is recorded but read by nothing: it aborts, so a record's existence proves
  it passed (DECISIONS 2026-08-03).
  - **Where a verb exists it is a method**, and the built surface is exactly `print`, `as.matrix`,
    `as.data.frame`, `cite_clocks` and `rbind`, plus the plain `calc_accel()`. Coverage is
    deliberately not `summary()`: it is `clocks_coverage()` (one row per **(clock, batch)**) and
    `samples_coverage()` (each sample's batch alongside its id), with `mc_batch_id` **last** in both
    -- it is the join key, but it is a hash, so it does not sit in front of `clock_id`. Citations
    dispatch as `cite_clocks()`, a **package-owned** generic, because `utils::citation` and
    `utils::cite` both exist as plain functions and taking either name masks it. `[`, `cbind` and
    `augment` are **unbuilt ideas, not contracts** (DECISIONS 2026-07-23/24/25, 2026-07-27).
    **`codebook()` is decided but unbuilt**: `data.frame(clock_id, description)` dispatching like
    `cite_clocks()`, blocked until upstream verifies a `description` per clock. Do not build it
    against a partly populated field (DECISIONS 2026-08-04).
  - **The batch column reaches an exit frame only when the record spans more than one batch** -- at
    one batch it is a repeated hash carrying no information. All four exits (`as.data.frame`,
    `calc_accel`, both coverage frames) share the **one** test in `is_multi_batch()`
    (`R/mc_result.R`), keyed on `length(unique(provenance[[mc_batch_id]]))`, the vector that fills
    the column and never `per_clock`'s names. Each **declines to build the column** rather than
    building one it will lose; `drop_single_batch()` still runs at all four as the gate, now a
    no-op, and the two forms produce `identical()` frames. The coverage frames thread the decision
    down to `panel_rows()` / `batch_coverage()`, which assemble by `rbind` and cannot add the column
    after the fact. The conditional schema is a real cost, accepted knowingly: the four exits must
    appear and vanish **together** or the two coverage frames disagree about whether the join key
    exists. Nothing internal is conditional (DECISIONS 2026-07-31, 2026-08-01).
  - **The label is `mc_batch_id` everywhere the user can touch it** -- both coverage frames, both
    finalizer frames, `$provenance`, and the `calc_accel()` formula namespace, where it is always
    present and always reserved against `data =` (supplying it there is an error, not a precedence
    question). `mc_`-prefixed because a user's `batch` is their slides or plates, a legitimate
    covariate that a bare name would shadow in a formula and collide with on a join (DECISIONS
    2026-07-31).
  - **Batch labels are derived, never assigned**: `construct_mc_result()` sets one to
    `batch_hash(pheno[[pheno_id]])`, the full 16 hex of `xxhash64` over the pheno's **id column
    only**, never truncated, and there is no `batch =` argument anywhere. **Hash the canonical
    form** (`sort(..., method = "radix")` + `unname`/`as.character` + `serialize = FALSE`) or the
    label follows the id *sequence*, the locale and the R serialization version instead of the id
    set. Do not hash the pheno *frame* (covariate values would fold into batch identity, and
    `digest` is storage-type sensitive), and do not restore an assigned label -- a stored string
    cannot say whether a human chose it, which is what sank `is_auto_label()`. Two batches sharing a
    label then needs a 64-bit collision against id sets gate 1 just made disjoint, so there is **no
    gate on labels**; do not add one back (DECISIONS 2026-07-30). `sim_DNAm(suffix =)` is **not** a
    batch argument -- it suffixes sample ids so two simulated blocks clear gate 1.
  - **`rbind` binds and labels; it never reconciles.** Nothing is re-imputed, no denominator merged,
    no cross-sample column recomputed unless `refinalize_clocks()` is called by hand -- but a
    multi-batch bind carrying non-empty `pending` **says so once** (`say_pending()`, an `inform` not
    a `warn`), naming the columns from `names(pending)`. `refinalize_clocks()` **reads `pending` and
    never consumes it**, so it composes in any order and a second call is a no-op. Its gates follow
    one line -- **record what batching forced, refuse what the caller chose differently** -- so a
    per-batch fill regime is recorded while overlapping ids, differing score columns, a differing
    `pheno_id` and a differing `normalize =` all throw. There is no `force =`. It mints, renames and
    renumbers nothing, and **drops** argument names (`unname(list(...))`) rather than refusing them,
    because `split()` names its result and refusing would kill
    `do.call(rbind, lapply(split(...), ...))`, the idiom the feature exists for. This is what makes
    re-association exact: `rbind(rbind(r1, r2), r3)` and `rbind(r1, r2, r3)` are `identical()`.
  - **A finalizer is any exit that takes an `mc_result` and returns something that is not one.** The
    test is mechanical, so the set is derived rather than listed and cannot go stale:
    `as.data.frame()`, `as.matrix()`, `calc_accel()`, `score_associations()`. All four re-finalize
    on the way out and say so, under `say_pending()`'s exact guard (non-empty `pending` **and** more
    than one batch); `rbind` is the only verb that leaves `pending` unresolved. It must **not**
    become a finalizer: `do.call(rbind, ...)` recurses and would re-finalize at every intermediate
    step, whereas a finalizer is a leaf handing back a value the record cannot be recovered from.
    Enumerating the set instead of deriving it drifted twice (DECISIONS 2026-08-03). A single-batch
    record is skipped because its reduction already spans its whole cohort; for the per-batch
    reductions, finalize each record *before* binding.
  - **The two batch counts are cross-checked, not chosen between.** `n_batches()` derives the count
    from `provenance[[mc_batch_id]]` and `stop()`s if `length(per_clock)` disagrees. Every finalizer
    and both coverage frames route through it, and `finalized()` calls it **before** the `pending`
    test so `&&` cannot short-circuit past the check. Never read either count directly.
  - **Every `print.mc_*` method shares one grammar, built in `R/print.R`**: a `<class> A x B` header,
    a `$component [what is shown]` line per list element, then `... N more <axis>`. The builders
    return strings, so a cli printer and a `cat` printer emit identical text without moving the cli
    boundary. A new record class reuses them (DECISIONS 2026-07-29).
- **Scores only, and the record remembers its inputs.** `$scores` is scores -- no auto-appended
  phenotype columns. Separately, `$pheno` carries the *aligned* pheno narrowed to the id column plus
  the covariates the run actually required. **It is never `NULL`**: with no `pheno =` supplied,
  `resolve_pheno()` materializes the id column alone, so `$pheno` is one shape everywhere and the
  batch label always has a column to hash (DECISIONS 2026-07-30). Its columns need no `rbind` gate
  -- they are `unique(c(pheno_id, covariates))`, both already pinned. **And it never carries row
  names**: `resolve_pheno()` has one exit and resets them there on both branches, because row names
  would be a second identity beside the id column and the two can silently disagree under subsetting
  or `rbind`. `rbind` adds no reset of its own. Align pheno by sample id, never row order.
- **Imputation in one place, never crossing sources.** Partial NA on a present probe -> cohort mean
  (shared cache); a fully absent probe -> the clock's vendored ref, or drop by policy.
- **A resolved panel is positions, and `facts[["usable_cols"]]` is the one axis.** Every
  `score_present_idx` / `norm_present_idx` is a position in that vector, `block[["usable_idx"]]`
  maps the same axis to the block's columns, and nothing in the scoring path looks a CpG up by name
  -- that lookup was 1.15s of a 3.0s run at 414k columns. A position is a bare integer, so **the
  wrong vector still yields valid positions**: code that sorts, uniques or subsets `usable_cols`
  between `resolve_cpgs()` and `mc_block()` scores the wrong CpGs and **returns a number**. Three
  rules keep that from happening, and a chunked front end is where they will be tested.
  `resolve_cpgs()` carries the vector it resolved against and `mc_block()` refuses one that is not
  `identical()` to it -- free, because R compares interned strings by pointer, and exact where a
  hash is not. `block_cols()` guards the positions themselves before indexing, because a `0`
  silently drops an element rather than erroring. And a subset filters names and positions with
  **one** mask, never two (`component_present()` in `R/score_default.R` is the pattern to copy);
  its `cols` follows the panel's order, not the coefficient vector's, and every consumer indexes by
  name off it (DECISIONS 2026-08-05).
- **Accessors are the executable schema.** `calc_clocks` consumes accessors (`get_clock`,
  `clock_coefs`, ...), never raw nested catalog lists. No hand-written `schema.md`.
- **Never `$` in `R/`. Always `[[`.** `$` partial-matches on lists, so a missing exact field
  silently resolves to a longer one (`entry$covariates` -> `covariates_required`) and the caller
  gets a wrong value, not an error. The rule is **blanket, not scoped to catalog reads** -- a scoped
  rule needs a judgement call at every site, and the list that is "obviously not a catalog entry"
  today is the one someone widens tomorrow -- so it binds on result records, `optim()` output and
  environments alike, where `$` is exact and harmless but not worth the exception.
  `tests/testthat/test-source-hygiene.R` enforces it by scanning **parse tokens**, not text, so a
  `$` in a comment or a `"\\.qs2$"` regex does not count and a real access cannot hide; the only
  exempt file is the generated `R/RcppExports.R`. `options(warnPartialMatchDollar)` is **not** the
  fix -- a package cannot set a session global for its users and it does not fire under R CMD check.
- **Never resolve a name by partial match. `$` is one instance of a general rule.** A token the user
  supplied is matched **exactly** against its closed set, so `pmatch()`, `charmatch()` and
  `match.arg()`'s abbreviation handling are all out at the front door. It resolves silently to
  something the caller did not name, and **the closed sets here grow with every sync**: `"Sys"`
  resolves to `SystemsAge` today, and the sync that adds a second `Sys*` group turns working user
  code into an ambiguity error. An abbreviation cannot be supported without freezing the catalog.
  The tool is **`checkmate::assert_subset()`** -- exact, names the offending element **and** the
  valid set, and takes `empty.ok` as an explicit flag (precedent: `list_clocks(tag =)`,
  `mc_resolve_groups()`). **`match.arg(several.ok = TRUE)` is not an alternative** and was measured
  before being rejected: it returns `choices[1L]` for `NULL` (here `"all"`, a silent mass download),
  it errors only when *every* element fails so a typo beside a valid token is dropped without a
  word, and it does not deduplicate. Its one wanted property (banning `character(0)`) is one line to
  write ourselves (DECISIONS 2026-08-03).
- **Never `<<-` in `R/`. Mutable state is an explicit environment.** `<<-` does not name a target:
  it walks the enclosing frames and assigns into the first one that already binds that name,
  **creating a global** if none does, so renaming or deleting a local silently promotes a local
  update to a package-level one, with nothing at the read site to say where the value came from. The
  replacement is a store you can point at: `st <- new.env(parent = emptyenv())`, write
  `st[["x"]] <- ...`, and have the mutator **return `st`** so readers take it as a value. Precedent:
  `miage_fit()` in `R/score_MiAge.R`. Outside `R/` this is a preference, not a rule --
  `data-raw/sync.R` and a condition collector in `tests/` still use `<<-` legitimately.
- **`checkmate` asserts at the exported surface; internals get bare `stop()`.** By the time a helper
  runs the value has crossed the front door, so a second assertion fires only when *we* have a bug,
  where a `checkmate` message aimed at a user is the wrong register. An export called by another
  export is a **trusted** caller (`predict_sex` -> `calc_clocks`), so the inner assertion is
  redundant; the exception is a value the inner one reads that the outer never validated, which is
  why `recorded_from_female()` keeps its assert. An internal guard that must stay -- a bounds check
  ahead of a kernel -- keeps its checks and drops to `stop()` with short greppable text
  (`check_moment_sets()` in `R/missingness.R` is the pattern).
  **`.var.name` is filled only where the deparse lies.** At a boundary the deparse already prints
  the caller's own word and a hand-written string is staleness risk. Fill it **iff the deparsed
  expression does not name something the caller can locate in their own call** -- not the same as
  "not a bare symbol": `assert_character(colnames(DNAm))` is left alone, while `check_pheno()`'s
  `pheno[[ID]]` names an internal parameter and carries `.var.name = paste0("pheno$", ID)`. **Never
  as a substitute for moving a check**: a wrong name is usually evidence the check sits in the wrong
  frame (DECISIONS 2026-08-03).
- **Accessors read declarations; they never search.** No `grep`/regex/fuzzy match over tensor names,
  clock ids, or file paths to find a payload. Resolve the declared pointer (component,
  `probe_sets[["scoring"]][["file"]]`, `imputation[["ref"]]`) and `stop()` when it is absent or
  ambiguous -- an accessor that cannot find its declaration has done its job by failing. Searching
  hides an upstream/sync gap and can silently return a sibling clock's tensor. **Three of those
  lists are keyed by sync, so the lookup is a lookup**: `components` by `name`, `probe_sets` by
  `role`, `recipe` by `out` (steps without one stay unnamed), with `key_declarations()` stopping the
  build on a collision -- which is why those accessors carry no multiplicity guard. The two
  remaining scans are the two whose predicate is not a key, `component_tensor()` on `row_key` and
  `stack_step()` on `op`, both of which genuinely collide in the shipped catalog, so `pick_one()`
  stays load-bearing there (DECISIONS 2026-07-28).
- **Coverage is never reported for a sample it is not true of, and a clock reports only what it
  counted itself.** A clock whose branch reads no betas -- assembled purely from other clocks'
  scores -- **has no coverage of its own**: `per_clock[[id]]` is `NULL`, it gets no `sample_miss`
  column and no `samples_coverage()` row, and its all-`NA` `clocks_coverage()` row says so.
  `clock_reads_cpgs()` (`R/score_cohort.R`) is the one source and switches on `score_type()`, so it
  is a fact about the closed branch set, not a clock list; today it selects the 7 sex-routed
  aliases, `GrimAgeV1` and `DNAmFitAge_{Sex}` (`GrimAgeV2` keeps its record -- its cox stack
  declares `internal` surrogates and it really does read its CpGs). **This loses nothing**: a clock
  that reads no betas can only be fed through its dependencies, so every CpG in its declared panel
  is already counted on a descendant that does. Do not "fill in" a `NULL` record with a merged
  figure and do not restore a stitched per-sample count for an alias -- read the descendants' rows
  for the denominators (DECISIONS 2026-07-29).
  - **`samples_coverage()` drops NA-coverage rows**, which has exactly one source: a routed member
    masked on a row its sex did not score. So the long frame carries one row per sample per family,
    under the model that scored it, and a sample no model scored (unknown sex) has no row at all --
    the same fact its `NA` score already carries. It was never a complete sample x clock grid.
  - **The converse binds too: never score a CpG coverage did not count.** Coverage counts the
    *declared* panel, so a branch takes its CpGs from the resolved `cpgs` it is handed
    (`score_present` / `score_absent`) and never re-derives them against the block's cohort-wide
    `usable` set -- that set is every panel's union, so a stray coefficient would resolve against it
    and be scored silently. A branch holding a bare coef vector (GrimAgeV2 surrogates, PhysAge
    surrogates) goes through `component_present()` (`R/score_default.R`), which intersects with the
    declared panel and `stop()`s if the component names a CpG the panel does not. The always-on
    smoke tier exercises that guard over the shipped catalog without naming a clock.
- **The callable pool is not the catalog, and neither is the output.** Clocks that exist only as
  routing targets (the 14 sex-resolved DNAmFitAge members) are internal machinery: scored, kept for
  coverage, **never a score column**, and a hard error if requested by name, pointing at their
  alias. A sex-routed family returns exactly one column per alias, populated for every sample --
  never a male column, a female column and NAs. The pool, the refusal, its suggestion and the output
  filter (`drop_routed_members()`) all derive from one source (`sex_routed_members()`).
- **No network at install/build/check/CRAN test.** Double-precision coefficients only.
- **No commit SHA / pin as result provenance.** Correctness is proven by fixtures.
- **Correlation is never a numeric gate. Anywhere, for anything.** Not in parity, not in a unit
  test, not as a "sanity check" alongside a real bound. `cor()` is offset- and scale-invariant, so
  it cannot distinguish "we match the oracle" from "we are uniformly wrong", and it concentrates
  near 1 for any monotone-ish agreement, hiding catastrophic per-sample outliers. A numeric
  agreement gate is always a **bounded per-element difference** -- absolute and relative, both,
  taken as a `max` -- so one bad sample fails the test. The same reasoning bans `median`/`mean` as
  the reducer over per-element differences.
- **Read `dev/WRITING.md` before editing any text a user can see. Reading it is the first step of
  the task, not a review at the end.** Scope is the whole user-facing surface, not only roxygen: doc
  blocks, `README.Rmd`, `vignettes/*.Rmd`, and every cli message. It is **not** a style guide to
  consult when unsure -- close to every rule in it is there because that exact failure has already
  happened here, and `lint_roxygen()` / `lint_seealso()` catch only the mechanical subset. Where a
  change and the file disagree, **update the file in the same pass**: a rule the shipped files
  violate is worse than no rule.
- **The exported surface is documented, and `dev/WRITING.md` is how.** Roxygen went on because the
  package went Rcpp and `useDynLib` has no route into `NAMESPACE` except a tag, so `NAMESPACE` and
  `man/*.Rd` are **generated files** and `devtools::document()` is a normal part of the workflow.
  Prose docs shipped 2026-08-03: all 26 user-facing topics carry real `@param` / `@details` /
  `@returns` / `@examples`, the `@seealso` groups are closed, and `cite_clocks` merges its three
  methods onto one topic with `@rdname`. `lint_roxygen()` and `lint_seealso()` (`R/dev-utils.R`)
  must both come back empty (DECISIONS 2026-07-27, 2026-08-03).
- **Never hand-edit `NAMESPACE` or `man/*.Rd` -- own the tags, not the files.** `document()`
  rewrites them from tags, **silently dropping** anything added by hand. This has bitten once: a
  hand-added `useDynLib(methylCIPHERv2, .registration = TRUE)` vanished on the next `document()` and
  took every compiled kernel down with it (`.Call()` -> "not available for .Call() for package").
  So a new export or S3 method gets a bare `#' @export` beside the function (precedent in
  `R/mc_result.R`), and package-level wiring (`@useDynLib`, `@importFrom`) lives in
  `R/methylCIPHERv2-package.R`. Then run `document()` and check the diff is only what you intended.
  The maintainer still owns the exported *surface*: say in your summary which `export()` /
  `S3method()` entries a change implies rather than quietly widening it.
- **A `.cpp` edit needs an explicit rebuild.** `devtools::load_all()` happily reuses a stale
  `src/*.dll`, so a changed or newly added kernel silently does not exist. Run
  `pkgbuild::compile_dll(".", force = TRUE)` first; the symptom of skipping it is the same
  "not available for .Call()" error above.

## sync.R workflow (`data-raw/sync.R`)

Pulls the scoring contract from `methylCIPHER-meta` into the package. Not run at build/check -- a
maintainer runs it and commits the regenerated `R/sysdata.rda`. **You do not need this to
contribute** (the catalog is committed). `sync()` needs read access to `methylCIPHER-meta` (private,
pre-release); `sync(upload = TRUE)` also needs a release-write token (maintainer-only).

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
       registry; do not add a code path.
     - Verify a sync change by **dry-running the build in memory first** (build catalog + bundles,
       diff every panel against the committed `R/sysdata.rda`) before regenerating.
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
  a **per-clock** field, so a group may be on both sides of the split (`Zhang2019` bundles its EN
  arm and packs its BLUP arm). `EXTERNAL_GROUPS` and `EXTERNAL_CLOCKS` both feed that one field, and
  everything downstream reads the field -- `split_group_ids()` puts a mixed group in both buckets,
  `build_group_bundles(external =)` partitions its tensors by declaring member (DECISIONS
  2026-07-29).
- **An external clock is not necessarily a pack-scored one.** `clock_is_external()` says where the
  weights live; `is_pack_scored()` says whether `score_pack_group()` computes the score. They agreed
  only while every external clock happened to be a batched weighted sum. Keep them apart: **needing
  a pack** (`pack_groups_needed()`, parity's `skip_if_no_pack()`) keys on externality, and **how the
  arithmetic runs** (`score_type()`'s group hooks, parity's relaxed `packs` tolerance) keys on the
  scoring path. `score_type()`'s external branch carries a group hook per group whose arithmetic is
  not a plain weighted sum -- `SystemsAge` (`center_scale`), `Zhang2019` (`sample_scale`) -- and
  those clocks reach their coefficients through `clock_coefs(id, packs)`, which reads the pack's raw
  tensor by the same `coef_path` the bundled arm reads out of `mc_bundles`. Do **not** hoist the
  `switch(gid, ...)` above the external check: that changes dispatch precedence for every clock in
  the catalog to solve a one-group problem (DECISIONS 2026-07-29).
- **"Assets" are the packs, and every public name is `<verb>_mc_<noun>`**: `get_mc_assets_dir()` /
  `set_mc_assets_dir()` (the setter `NULL`-clears and returns the old value invisibly),
  `list_mc_assets()` (read-only table), `download_mc_assets()` (bytes -> disk), `load_mc_assets()`
  (-> RAM), `clear_mc_assets()` (delete). That is the whole pathing surface; there is no bare-noun
  accessor. **Every read-only question has a read-only answer**: `list_mc_assets()` reports size /
  `downloaded` / `superseded` per group without prompting, fetching or deleting. The word "cache" is
  reserved for the unrelated internal `partial_cache` (cohort-mean fill).
  The dir stays `tools::R_user_dir(..., "cache")` -- but **not because `which = "data"` would breach
  CRAN policy. It would not**: policy permits data, configuration and cache alike under the same two
  conditions, so neither discriminates and `clear_mc_assets()` is required either way. It is a
  **platform call**: `"data"` is `%APPDATA%` on Windows and roams, `"cache"` is `%LOCALAPPDATA%` and
  never does; the macOS counter-risk costs one consented re-download. **Durability is the job of
  `MC_ASSETS_DIR` / `set_mc_assets_dir()`, not of the default** (DECISIONS 2026-08-04, superseding
  the "policy violation" reading of DECISIONS 2026-07-24).
- **Assets move in both directions under one consent rule.** `load_mc_assets()` /
  `download_mc_assets()` fill the dir and `clear_mc_assets()` empties it; all three take `ask`,
  prompt interactively, **refuse** non-interactively, and treat `ask = FALSE` as the explicit
  consent signal. Nothing is fetched or deleted unprompted -- CRAN requires a supported way to
  reclaim `R_user_dir()`, so `clear_mc_assets()` must stay a real delete, not a report. **Clear
  means clear** (`pak::cache_clean()` semantics): it removes the currently declared packs **and**
  every superseded one, with no opt-in flag, because content-addressed filenames orphan the old file
  on every `payload_hash` move and leaving those behind made `clear` fail to reclaim. The consent
  gate is what makes this safe -- the prompt counts the two kinds apart and lists the files first,
  capped at `MC_MSG_CAP` with any remainder counted on its own line (DECISIONS 2026-08-03). The
  stale scan is not a search for a payload: the stem comes from the declared `file` field, so only
  the hash is a wildcard and a foreign stem is never touched.
  **The gate argument fails closed.** `ask` is a strict flag: only `FALSE` consents, and anything
  that is not a single non-NA logical is an error, never permission. `ext_data` reaching
  `mc_resolve_assets_dir()` is a path or `NULL` only -- a loaded pack names no directory, so it
  stops rather than falling back to the default. Both were silent widenings of "permission" once
  (DECISIONS 2026-07-23); do not re-introduce an `isTRUE()`-style test on either.
- **One argument for the external data, one noun for the thing.** `ext_data` (on
  `load_mc_assets()`, `calc_clocks()`, `sim_DNAm()`) is `NULL` (open set, may download), a path
  (**closed set**, never downloads, missing is fatal), or loaded pack(s). Resolution order is
  `ext_data` > `mc.assets_dir` option > `MC_ASSETS_DIR` env > `R_user_dir` default.
  `download_mc_assets()` / `clear_mc_assets()` take **no** dir argument -- use the setter. The name
  is `ext_data` because a four-letter English preposition cannot be grepped (DECISIONS 2026-07-24,
  2026-07-28).
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
  `score_type()` on every clock, so an unroutable entry fails the tier here. **Its value is not "a
  clock stopped running"** (parity proves that over a superset): it is the only tier that runs
  `calc_clocks()` in the **default configuration** -- both coverage gates on, full panels, no meta
  repo -- and the only caller of `sim_DNAm()`. Parity scores with both gates at 0 and is skipped on
  CRAN, so it can never stand in for this tier (DECISIONS 2026-07-24).
- **Value goldens (always, no meta dep):** hand-authored engine/machinery unit tests with goldens
  written in-test, one per scoring path (linear sum/mean, sex-split, imputation offset, bundled
  composites). External-pack scoring is smoke-only here; parity owns those goldens. **Nothing in
  this tier reads a third-party package**: the one golden that did -- the `DunedinPACE` reference --
  moved to the parity tier on 2026-08-04, so what remains here runs anywhere.
  - **"Always" means every `devtools::test()`, not every CRAN run.** Since the 2026-08-04 trim
    (1284 -> 799 expectations) the internal half of this tier carries `skip_on_cran()`, so CRAN
    runs 391: the smoke tier, the front-door refusals, and the `mc_result` contract. **CRAN
    therefore applies no numeric gate at all** -- parity was already skipped there, so a green
    CRAN check proves the package loads, refuses correctly and returns a well-formed record, and
    proves nothing about the scores. Do not read it as more than that, and do not "restore
    coverage" by ungating: the arithmetic is parity's job.
  - **`skip_on_cran()` is not "does not run in pre-submission checks".** `NOT_CRAN` is unset on
    r-hub and on a GitHub Actions `R-CMD-check`, so the gated tier **does** run there; what the
    flag actually buys is CRAN's own machines not paying for it. Treat r-hub as the place the
    internal tier gets exercised across platforms, and do not reason about a gated test as though
    it only ever runs on the maintainer laptop. **`devtools::check()` is the opposite case** --
    it sets `NOT_CRAN=true`, so it runs the *whole* 799 and a `skip_on_cran()` there is inert.
    A local `devtools::check()` and a plain `R CMD check` on the tarball therefore run different
    suites; say which one a result came from.
  - **The gate is per `test_that`, first line, never at file level** -- testthat runs top-level
    code at collection, so a file-level skip reads as one skipped file instead of N skipped tests.
    Under `NOT_CRAN=false` the suite reports 89 skipped blocks across 24 files; a file-level gate
    would show 24. Anything calling a non-exported function is gated by default.
  - **Four goldens are kept against "parity owns it", because parity is structurally blind to
    them**: alias routing (fixtures exist on the 14 routed members, never on the 7 aliases, so
    *which sex's model scored which sample* is uncovered), DunedinPACE quantile normalization (the
    only always-on proof normalization is applied), the PhysAge mean-divisor fill offset (parity
    scores clean panels), and Wang mixed-request domain isolation (parity scores one clock per
    call). Do not delete these as parity-redundant; they are not (DECISIONS 2026-08-04).
    The `test-normalize.R` BMIQ golden was the fifth such keep for one day, then **moved into the
    parity tier** as the `parity (horvath normalized)` block later on 2026-08-04. What stays in
    `test-normalize.R` is the *record* half -- `provenance$normalized`, `cov$normalizes`, the
    `sample_miss$norm` column -- which parity does not look at.
  - `test-sim-smoke.R` is **ungated and untouched** by the trim.
  - **A test whose subject is the source tree does not ship: `.Rbuildignore` it, do not
    `skip_on_cran()` it.** Two files qualify and both are ignored as of 2026-08-04 --
    `test-fixtures-parity.R` and `test-source-hygiene.R` (precedent: `R/dev-utils.R`, ignored
    already). Any later `lint_roxygen()` / `lint_seealso()` test joins them. `skip_on_cran()` is
    the wrong tool here for two independent reasons, both measured against a real
    `devtools::check()` run:
    - **It does not fire under `devtools::check()`**, which sets `NOT_CRAN=true`. The hygiene
      scans ran in the installed package and **failed** -- `scan_sources()` returns `NULL`, not
      `character(0)`, over an empty file list, so they do not even pass vacuously.
    - **It cannot fix an unstated-dependency WARNING at all.** `checking for unstated
      dependencies in 'tests'` is a **static scan of the shipped sources**, so a `duckdb::` that
      never executes still counts. `duckdb` and `DunedinPACE` are deliberately undeclared
      (maintainer-gated tier), and removing the file from the tarball is the only thing that
      silences it. `DBI` is in `Suggests` and was never the problem.
- **Cohort-gated parity fixtures** (science gate; the only clock-golden source for a clean panel):
  run against **every registry cohort** -- `data-raw/methylCIPHER-meta/fixtures/{cohort}/beta.duckdb`
  for `cohort_EPICv1` and `cohort_450K` -- which need BOTH `MC_PARITY=1` and that cohort staged
  (`file.exists()`). Upstream ships one `fixtures[]` block per cohort; each (clock, cohort)
  pair is its own test. Run locally via the dev-only `test_parity()` (`R/dev-utils.R`). CRAN skips
  this tier; CI must stage the cohorts, set the flag, **and install `duckdb` itself** -- it is
  undeclared as of 2026-08-04, so DESCRIPTION will not pull it (DECISIONS 2026-08-04).
  - **Never run this tier unless the user explicitly asks for it.** `test_parity()`, and any
    invocation that sets `MC_PARITY=1`, is minutes-long and reads the staged duckdb cohorts. It is
    not part of "run the tests": the default `devtools::test()` (parity auto-skipped) is. Verify a
    change against the always-on tiers, say that parity was not run, and let the maintainer ask.
  - **The gate is on the generator, not on the generated test.** `staged_cohorts` is the file's one
    switch: `parity_targets()` returns nothing for a cohort that is not staged, and the PhysAge,
    census and Dunedin blocks are guarded the same way, so **a tier that cannot run emits one skip
    instead of 264**. What was skipped is still reported -- one line for the tier being off, one per
    unstaged cohort -- because a per-target skip is 263 copies of a reason that is the same every
    time, and testthat prints every one of them with its location. Do not put the flag test back
    inside `run_parity_target()`: the count of *generated* blocks is how a maintainer reads a parity
    run, and it must not silently include tests that never had a cohort to read (DECISIONS
    2026-08-04).
  - **The tier carries one non-cohort test**: the degraded-coverage golden against the `DunedinPACE`
    reference package (`danbelsky/DunedinPACE`), which builds its own holed panel, so it needs the
    flag but no duckdb and no staged cohort. It is the **only** gate on that path (DECISIONS
    2026-07-29). It is emitted only when the flag is set, with **no `skip_if_not_installed()` on the
    reference**: the tier only ever runs on a maintainer machine that has it, so a skip there would
    hide a silent non-run rather than protect anything (DECISIONS 2026-08-04).
  - **Two axes, both gated.** Every fixture must clear **`max_abs_diff` AND `max_rel_diff`**, not
    either: the absolute bound is the only one with meaning near zero, the relative bound the only
    one with meaning at large magnitude. Both use **`max`, never `median`**. **`PARITY_REL_TOL` is
    `1e-10` everywhere, with no per-block exception** -- it is scale-free, so there is never a
    units-based reason to move it. Only `PARITY_ABS_TOL` is per-block, and relaxing it is a
    statement about units, never about correctness (DECISIONS 2026-07-25).
  - **Four blocks, derived from the catalog.** `parity_block()` sends each clock to `horvath`
    (declared `fixtures[].oracle == "horvath_online"`), `packs` (`clock_is_external()`), `fitage`
    (group `DNAmFitAge`), or `core`; `parity_targets(block)` builds one loop per block over the
    shared `run_parity_target()` body. There is no clock list -- every block reads a declaration, so
    a regenerated fixture retires its block automatically. `PARITY_ABS_TOL` is
    `c(core = 1e-10, fitage = 1e-10, packs = 1e-6, horvath = 1e-10)`; the pack relaxation is
    measured, not guessed, and is a statement about magnitude, not about agreement.
  - **The `horvath` block is skipped, and that is a finding, not a shrug.** The oracle filled every
    completely-absent probe server-side with an unpublished per-probe constant and BMIQ'd its panel
    for `DNAmAge` only, so against the *submitted* matrix the residual tracks the absent-probe count
    and nothing else -- pairs with zero absent probes agree to ~1e-8 relative, which for a matmul
    proves the tensors and the engine are right and puts the divergence in the oracle's input.
    **Do not "fix" this with a tolerance**: the residual spans 4.2e-08 to 2.7e-01, so any bound wide
    enough is vacuous (DECISIONS 2026-07-25). **`Horvath1` is the one exception to that reading**,
    because it is the one the oracle BMIQ'd: of the 15, 13 declare `scheme = none` and `Horvath2`
    declares the inexpressible `noob`, leaving `Horvath1`'s `bmiq` as the only scheme we can apply.
    The four-block loop still scores it with `normalize` at its opt-in default of **off**, so its
    gap there is a normalization gap, not a fill gap.
  - **That exception is now its own generated block**, `parity (horvath normalized)`, admitted
    2026-08-04 -- the third tolerance regime CLAUDE.md previously said had not been decided.
    Membership derives from `is_normalized_horvath()` (horvath-online **and**
    `clock_norm_scheme() %in% NORM_SCHEMES`), never a clock list, and the block runs the clock with
    normalization **on**. **The cohort split is data-derived, not declared**: the test skips unless
    the cohort leaves *zero* scoring probes absent, because only then can the oracle's undisclosed
    fill not contaminate the comparison. Measured: normalizing takes `cohort_450K` (complete panel)
    from 7.7 years of disagreement to 0.11, while `cohort_EPICv1` is 19 probes short and stays at
    4.0 -- which is the fill gap talking, hence the guard rather than a second tolerance.
    `HORVATH_NORM_TOL` is a **snapshot of that residual, not an agreement target**, keyed
    `clock@cohort` and sitting just above the measurement (BMIQ is deterministic here -- verified
    bit-identical across three runs). A pair that clears the absent-probe guard with no entry in
    the map **fails**: a newly admissible pair needs its residual measured, never defaulted
    (DECISIONS 2026-08-04).
  - **A fixture is scored on the panel the oracle used, which is not always the scoring panel.** A
    recipe declaring a `sample_scale` op z-scores each sample over **every** probe in the input
    matrix, so feeding it the union of scoring panels moves each sample's mean/sd and the score with
    it. `needs_full_panel()` reads that op off the declared recipe -- never a clock list, today only
    the two `Zhang2019` arms -- and the target loads the whole array with `cohort_betas_full()`
    instead (DECISIONS 2026-07-25).
  - **The blocks are generated, so a dropped fixture is silence, not a failure.** `parity_targets()`
    loops over `clock_fixtures()`; a fixture upstream drops emits **no test at all**, just two fewer
    passes in a green run. One ungated-by-cohort census test guards the generator: every catalog
    clock declares a fixture for **every** `PARITY_COHORTS` cohort, except the 7 sex-routed aliases,
    whose 14 members carry them (both halves derived, never listed). It needs no duckdb, so
    `test_parity()` runs it even where nothing is staged -- but it **is** behind `MC_PARITY`, so a
    plain `devtools::test()` does not catch a dropped fixture; CI does (DECISIONS 2026-07-26).
  - **Standing state with both cohorts staged: 266 blocks / 0 fail**, the block count raised from
    264 by the two `parity (horvath normalized)` targets on 2026-08-04 (stage one fewer cohort and
    it drops, by design). testthat counts *expectations*, not `test_that` blocks, and
    `expect_parity()` carries three (all-finite, abs, rel): 228 targets x 3 + PhysAge 2 x 6 +
    census 3 + the Dunedin reference golden's 8 + the normalized-horvath 450K target's 3.
    **The skip count depends on what else is cached, so check it against a cause before reading
    anything into it.** With packs cached it was 32 skip / `PASS 707` on 2026-08-02, so 33 / 710
    now; measured 2026-08-04 on a machine with **no packs cached**: 266 blocks / 91 skip /
    `PASS 536` / 0 fail, the 91 being 56 packs + 30 horvath-online + 2 Wang gaps + 2 Zhang BLUP +
    the normalized-horvath EPICv1 guard. Read a parity run by its **fail and skip** counts,
    checked against each other, before concluding anything from the pass number; the fail count is
    0-or-bust (DECISIONS 2026-08-02, 2026-08-04).
  - `KNOWN_PARITY_GAPS` (clock- or `clock@cohort`-keyed) holds only genuine skips -- **two** today,
    both `DNAmSex_Wang_*@cohort_450K`, whose deposited matrix carries no sex-chromosome probes, so
    the panel is 0% present and the fixture is the oracle's empty-panel `0`. **Do not relax
    `check_coverage()`'s `ratio == 0` stop to make them pass**: a 0 there is the `Female` quadrant
    of the sign map, not a small number. `KNOWN_PARITY_GAP_GROUPS` (group-keyed) is empty but stays
    a **separate** map, because group ids and clock ids share a namespace (`DNAmFitAge` is both) and
    one flat map could not say which a key meant.

### Test altitude -- keep tests loose enough to move fast

Assert what `calc_clocks()` *produces*, not how it is wired. A test that breaks on a no-behavior
refactor is too tight -- loosen or delete it.

- **The suite has a budget, and adding to it is a cost.** It reached 1284 expectations and was cut
  to 801 on 2026-08-04; the failure mode is a routine change touching a dozen files without buying
  proportional safety. Before adding a test, ask what it catches that parity, the smoke tier and
  the existing block do not. **Delete on sight**: a golden for a clock that has a parity fixture, a
  loop restating one invariant over many clocks (assert it over a representative shape instead), a
  block whose only failure mode is R itself being broken, and a second entry point to a path that
  provably calls the first (DECISIONS 2026-08-04).
- **Only `R/` is under test.** `tests/` covers package code and the data it ships, never
  `data-raw/`. `sync.R` is maintainer-side tooling against an upstream contract that upstream gates
  in its own suite, and reaching it means sourcing a file the package does not ship. Do not source,
  parse or otherwise bind anything from `data-raw/` in a test.
- **Errors: assert *that*, not the wording.** `expect_error(expr)` with no regex. Pin a message or
  condition class only when a test must otherwise confuse two distinct failure modes.
- **Never `expect_identical()`. Always `expect_equal()`.** `identical()` is bit-exact on doubles and
  also fails on differences that carry no meaning here (integer vs double storage, a dropped
  attribute, a name reordering), so a result correct to every digit anyone can act on fails in a way
  that looks like a numeric regression. `expect_equal()` applies a tolerance and is the right
  altitude for everything here, counts included.
- **No internal dispatch-tag tables.** Do not hard-code `clock_reduction()` / `score_type()` per
  clock; prove routing through output. The one allowed invariant: every catalog clock maps to a
  *known* tag -- and since `score_type()` stops otherwise, that test also proves the catalog routes.
- **No maintainer-side plumbing shapes.** Do not assert asset filenames, release tags, download
  URLs, or cache-dir order -- none reach a result. Test behavior (verifies on fetch, leaves no
  scratch, warns-not-stops on hash drift, closed set never downloads).
- **Re-derive a recipe in-test only until parity covers it.** Once a clock has a passing parity
  fixture, that fixture owns the numeric golden and only a smoke stays.
- **Coverage counts and provenance flags are output** -- asserting
  `res$coverage$per_clock[["Hannum"]]$score_imputed_full` or `res$provenance$dependencies` is fair
  game.
- **Minimize test-helper files.** A fixture builder/mock lives atop the one test file that uses it;
  promote to `helper-*.R` only when >= 2 files genuinely share it. `helper-fixtures.R` is the one
  such file: `mc_pheno()` / `mc_ages()` (9 files), `mc_fake_cpgs()`, and the `grip_fixture()` /
  `gait_holed_fixture()` pair shared by `test-coverage-report.R` and `test-score-fitage.R`.
  `sim_DNAm` / `random_betas` are package functions in `R/`, not test helpers.
- **Cohort/duckdb parity lives in one file** (`test-fixtures-parity.R`): one file-scoped read-only
  connection **per staged cohort** behind the `MC_PARITY` + `file.exists()` guard, torn down with
  `withr::defer(..., testthat::teardown_env())` -- not a module-global caching env. The file also
  holds the one test that needs the tier flag but no connection (the Dunedin reference golden), so
  the two gates are separate: `parity_on` is the flag alone, `staged_cohorts` is the flag plus a
  connection, and each generator picks the one it needs.
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

**The line is audience, not transport.** A message about **input the user chose** is `cli`, whatever
function raises it. A message about a **package defect** -- a "cannot happen" condition, a
catalog/sync gap, a missing dispatch branch -- is a plain `stop()` with `call. = FALSE`. This
replaced an enumerated cli keep-set, which had put messages about a `calc_clocks()` argument on the
`stop()` side, where they could not carry markup and so could not meet the writing rules at all
(DECISIONS 2026-08-03).

So today: assets lifecycle, discovery printers, public S3 refusals, and the whole `calc_clocks`
front door (token resolution, `validate_inputs.R`, coverage gates, `missingness.R`, `mc_cohort`,
`clock_cpgs`, `resolve_normalize`) are cli. Accessors, score branches, pack dispatch, catalog/sync
bugs and citation internals are `stop()`. **`list_clock_tags()` is not a printer** -- it returns the
registry as a value and prints nothing (DECISIONS 2026-07-29).

A defect message is **hard-coded and greppable**: a fixed prefix leading, values appended after it,
so a bug report can be located from the pasted text with no stack trace asked of the user. Not a
fully constant string -- `check_moment_sets()`'s failing index is real debugging value.

### How the text itself is written: `dev/WRITING.md`

**`dev/WRITING.md` is the single source, and this file does not restate it.** Read it before you
write or edit any text a user can see. It is tracked, so it resolves in a fresh clone. It holds:

- **R1 to R8**, the English rules. They bind **every message a user can see** -- cli message text
  and roxygen prose alike. They do **not** bind code comments, dev-facing `stop()` text,
  `data-raw/`, or `dev/` docs, so the ASCII section above is untouched and `--` stays required there
  (DECISIONS 2026-08-03).
- **The cli mechanics**: why `sprintf` output must never become cli input, why `bullets()` escapes
  braces, cap-before-format, `cli::qty()` on every plural marker, `cli_verbatim()` for anything
  pre-aligned, and never a multi-line `askYesNo()` prompt. Each is a bug that has actually happened
  here.
- **The roxygen template**: tag order, the `DOC_TYPES` param vocabulary, the shared-parameter donor
  and its one footgun, the closed `@seealso` groups, and the example rules.
- **`say_*` emits to the user; `note_*` records into the block's collector.** Do not use `note_` for
  something that prints.
- **The audit section**: the known-good exceptions an independent reader will otherwise report as
  defects, and the three CRAN shape rules the manual currently satisfies.

Two rules from it are repeated here only because the test suite and the workflow enforce them rather
than prose: tests assert *that* a message errors and never its wording (see "Test altitude"), and
`lint_roxygen()` plus `lint_seealso()` must both be empty before a doc change is done.

## Comments

- **Plain `#` comments are for the code; roxygen is for the manual.** Both are live. The exported
  surface carries real roxygen prose as of 2026-08-03, written to `dev/WRITING.md`.
- Keep `#` comments **short**: 1-2 sentences on *what* the code does, not a rationale essay.
- The *why*, and every decision or reversal, goes only in `dev/DECISIONS.md`.

## Source-of-truth docs (`dev/`)

The `dev/` folder is local-only **except** these four, which are tracked:

- `dev/DECISIONS.md` -- append-only, newest-first, date-stamped log of *why* / reversals (2026-07-30
  and later). Add an entry when a decision reverses a prior approach or is likely second-guessed; do
  not restate rules already stated here.
- `dev/DECISIONS.old.md` -- full pre-2026-07-30 decision history. Dated citations earlier than that
  cut resolve here; do not restate that archive in the live log.
- `dev/WRITING.md` -- the single source for how user-facing text is written. See "CLI messages"
  above; this file points there and does not restate it.
- `dev/to-do.md` -- queued work, tracked since 2026-08-04. A **staging area, not a record**: an item
  that becomes a design commitment gets a DECISIONS entry when it lands, an item that becomes a rule
  moves here, and a shipped item is deleted rather than marked done. Read it before starting new
  work; the pre-alpha section is what blocks a public release.

**There is no live design doc, and that is deliberate.** `migration-plan.md` and `detail-plan.md`
were retired on 2026-07-28, and `id-streaming-plan.md` (the chunking / binding / `prep()` design) was
deleted on 2026-08-02 once its shipped half was covered by the invariants above. **Do not
reconstitute any of them.** Built behavior is specified by the invariants plus the code; a separate
long-form spec of shipped behavior is a copy that rots. `sec N` citations in older DECISIONS entries
point at those retired files -- read them out of git history, not as live references. What upstream
declares (coef-path rule, declared-path set, tensor `row_key`/`col_key`, recipe operand namespaces,
the panel rule) is **not** restated in a `dev/` doc: `data-raw/sync.R` is self-documenting and is the
only source for it. Read `sync.R` itself before touching `sync.R`.

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

- **There is no tracked `.Rprofile`, and there must not be one again.** It attached `devtools` +
  `testthat`, which is convenient locally and fatal everywhere else: R sources that file for every
  `Rscript` started in the repo root, including the one CI runs to set the library path **before a
  single dependency is installed**, where a bare `library()` is a hard error that fails the job.
  That is exactly how the pkgdown workflow died on its first run (2026-08-04). Guarding it on
  `interactive()` was the small fix and was rejected for the general one: a startup file that runs
  before the environment exists is a hazard whatever it contains, and a repo-wide one makes every
  contributor's R behave unlike a clean session. **A personal `.Rprofile` is gitignored** -- keep
  machine-specific startup there, same as `CLAUDE.local.md` -- and `^\.Rprofile$` stays in
  `.Rbuildignore` so a local one can never reach the tarball. Every R invocation in this repo now
  starts profile-free with or without `--vanilla`.
- Put machine-specific or personal notes (OS, shell, local paths, private scratch) in
  `CLAUDE.local.md` -- gitignored, loaded automatically, never reaches a collaborator.
