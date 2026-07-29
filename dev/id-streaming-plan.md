# Chunked scoring, prepared inputs, and record binding

The one tracked design doc, and it covers **only work that is not built yet**. Shipped behavior is
specified by `CLAUDE.md`'s invariants plus the code; `migration-plan.md` and `detail-plan.md` were
retired on 2026-07-28 and must not be reconstituted (DECISIONS).

Cited by four `dev/DECISIONS.md` entries **by phase number** -- see "Phase identifiers are stable"
below before renumbering anything. Section numbers are cited by nothing and may move freely.

Three problems share one seam (Phase 5) and are therefore designed together: chunking, binding, and
the prepared-input record `prep()` (Phase 7). Per-clock normalization was originally folded into
Phase 7 and has since been split out as a plain `normalize=` argument that shipped on its own
(sec 9.3).

## Status

**The seam is cut, the sample axis is classified, and every clock is chunk-safe. No chunk source is
built, so nothing here is reachable from a user session yet** -- `calc_clocks()` behaves exactly as it
did, and every internal below is unexported.

| | |
|---|---|
| Built | Phase 1 (by other means), **Phase 2** (sec 6.1), **Phase 3** (sec 6.2), **Phase 5** (sec 4), Phase 7's `normalize=` half (sec 9.3) |
| Design only | Phase 4 (sec 8), Phase 6 (sec 5), the `prep()` record (sec 9) |
| Proven | chunk invariance over **every** clock, always-on tier, `test-chunk-invariance.R` (sec 4.2, 10) |

**Picking this up cold?** Read sec 3 (phase table) and sec 4 (the seam) first. The seam exists now,
so Phase 6 adds a front end and two adapters rather than re-cutting anything. The normalization work
is done and is *not* a prerequisite for any of it; it shipped on the same seam only by coincidence of
ownership. Sec 9.4 lists what is left near `prep()`.

**No live correctness gap.** Phase 3 moved the two PhysAge reductions past assembly, so a chunked
front end may score the whole callable pool. `spec$cross_sample` stays non-empty -- it is the
catalog's declaration and those clocks really are cohort-reducing -- but it now names the clocks the
scoring loop *defers*, not clocks it gets wrong.

---

## 0. The two problems are not one problem

| | Chunking | Binding |
|---|---|---|
| Goal | cap resident set size | reuse the `mc_result` verbs over records scored separately |
| Speed | may be sacrificed (may even win near the thrash boundary) | irrelevant |
| Fill values | **one** cohort mean set, shared by every chunk | **one set per record**, deliberately different |
| Result | numerically equal to a single-pass run | a labelled union of differently-imputed batches |

Binding is required *by* chunking but is strictly wider than it. The motivating non-chunk case: three
time points, not nested, imputed per time point on purpose, then bound to use `[`, `cbind`,
`augment`, `clocks_coverage`. Both cases produce disjoint id sets and identical clock columns, so
**the id gates cannot tell them apart** -- only recorded batch provenance can (sec 8).

## 1. The load-bearing rule: `calc_clocks()` takes a complete cohort

`calc_clocks()` assumes its `DNAm` is the whole cohort, always. Cohort-completeness is a
**precondition, not a parameter**: no `chunked =` flag, no chunk-awareness on the public surface, no
scorer asking whether it is seeing everything. Chunking is a separate front end (sec 5).

Why it has to be this way: partial-vs-complete missingness is only definable against a column's NA
count over every sample, so any function that classifies missingness needs the cohort.

### 1.1 The cohort facts reach the scoring internal as arguments

Both front ends call a shared **scoring internal**, `score_cohort()`, and every cohort fact reaches
it inside `facts` -- chiefly `partial_fill`, a named numeric vector whose **names are the
cohort-partial columns** and whose **values are the cohort means**.

`partial_fill` does two jobs, and the second is the load-bearing one:

1. it supplies the fill values;
2. **its names are the column classification**, which is a cohort fact and must not be re-derived
   from a block.

There is exactly one producer and it is never optional: `mc_cohort()` derives `partial_fill` and
`score_cohort()` consumes it, on both front ends. What differs between them is only where
`mc_cohort()` got its numbers -- one traversal of a whole matrix today, accumulated block statistics
under Phase 6 -- and `score_cohort()` cannot tell. It never classifies a column itself, so there is
no local computation for a caller to correct afterwards and no "supplied vs derive it yourself"
branch to get wrong.

Correspondingly, `scan_missing_cpgs()` is called by `mc_cohort()` and by nothing else. Its row half
(the all-NA-sample abort, `row_miss == ncol`) is per-row and correct against a block; its column half
is a cohort fact and is why it belongs upstream.

Per needed column, as a block sees it:

| Cohort-wide state | In `partial_fill`? | Block sees | Same as single pass? |
|---|---|---|---|
| clean | no | raw, no NAs | yes |
| partial NA | yes | NAs filled from the cohort mean, into the cache | yes |
| all NA | no | absent from `usable_cols` -> vendor ref / drop | yes -- the all-NA set is cohort-invariant |
| partial cohort-wide, **all-NA within this block** | yes | every row filled from the cohort mean, every row counted as imputed | yes -- **this is the case the vector exists to get right** |

A column absent from `colnames()` in one block but present in another is not covered by any row of
that table, which is why colname agreement is an error and not a repair (sec 5.3).

**Pre-imputing the block outside the scorer is rejected**, and it is the one thing that would break
this. It is score-equivalent but it erases the record of what was imputed: `count_sample_miss()`
counts NAs in the **raw** matrix over the cached columns, so a pre-filled block reports
`score_imputed_partial = 0` for every clock and an all-zero `sample_miss`, and `check_row_coverage()`
degrades into a per-clock constant. That breaks "coverage is never reported for a sample it is not
true of" **silently**, which is the worst available failure mode. Handing the block over raw, with
the fill applied inside, keeps the NA pattern visible exactly where coverage reads it.

## 2. Imputation policy (unchanged, and now load-tested)

Partial NA on a present probe fills from **the user's own data only** -- the cohort mean, or the
chunk-set mean when the caller deliberately batched. A fully absent probe takes the clock's vendored
ref or is dropped by policy. This is already a CLAUDE.md invariant; what is new is that streaming
pressure tested it and it held.

The reason a cohort-independent fill is not an option even when it would make streaming trivial:
a vendored constant is a **different estimand**, not a cheaper estimate of the same one. Cohort-mean
fill leaves a technical missingness problem technical; a fixed external constant converts it into
biological error, silently, in a value that then flows into a clock. So the two-pass design is not a
performance preference -- it is the only exact way to honor the policy.

## 3. Phase identifiers are stable

Phase numbers are **names, not an order.** Four tracked citations refer to them; renumbering
silently invalidates prose in files that never mention this one.

| Phase | Meaning | Status |
|---|---|---|
| 1 | sample-id resolution | **done** by other means -- mandatory rownames (DECISIONS 2026-07-24) subsumed it |
| 2 | kind-1 / kind-2 split (chunk-safe vs cohort-reducing clocks) | **done** -- `split_cross_sample()` -> `spec$cross_sample` |
| 3 | every cross-sample op leaves the scoring loop | **done** -- `finalize_cross_sample()` |
| 4 | `rbind` gates | not started |
| 5 | the seam: split `calc_clocks()` so both front ends compose the same internals | **done** -- `R/score_cohort.R` |
| 6 | chunk sources + the two passes | not started |
| 7 | `prep()` -- the prepared-input record | not started; its `normalize=` half split out and shipped |

Execution order was **5 -> 2 -> 3 -> 6 -> 4**, with 7 available any time after 5. Phase 7 is
orthogonal to chunking but lands on the same seam, so it is designed here rather than in a doc of its
own. **Next is 6**, then 4. Phase 6 does not depend on Phase 4: the chunked front end assembles
*fragments* and never calls `rbind.mc_result`, so `rbind` can go on refusing while chunking ships
(sec 8).

## 4. Phase 5 -- the seam (built)

`calc_clocks()` is three pieces, divided by what each depends on, so both front ends compose the same
internals. They live in `R/score_cohort.R` and none of them is exported:

```
mc_spec(clocks, pheno_id, normalize, from, ask)  # data-independent
mc_cohort(DNAm_or_source, spec, pheno, ...)      # one scan -> the cohort facts; gates
score_cohort(DNAm, spec, facts)                  # a matrix and the facts -> scores + intermediates
finalize_cross_sample(scores, pending)           # the one cohort reduction, after assembly
```

`calc_clocks()` is those four plus `construct_mc_result()`; `calc_clocks_chunked()` is the same four
with `mc_cohort()` accumulating over blocks and `score_cohort()` called per block.
`score_cohort()` never learns which front end called it and `facts` is the only channel between them,
which is what makes sec 4.2 a property rather than a coincidence. The wrapper calls the internal,
**not** the public `calc_clocks()` -- otherwise it pays per-block resolution and asset loading and
emits per-block error messages.

`facts` carries `sample_id`, the aligned `pheno`, `usable_cols`, `cpg_list` and `partial_fill`.
`spec` carries the resolved clock ids, the compute `sequence`, `output_ids`, the `normalize`
decisions, the covariate union, `pheno_id`, the loaded `packs`, the `panels`, and `cross_sample`
(sec 6).

`score_cohort()` takes no pheno argument: it narrows `facts$pheno` to the rows in hand by matching
`rownames(DNAm)` against `facts$sample_id`, so a block can never be handed a misaligned pheno. That
match is the identity permutation on a single pass.

**Coverage does not move.** It is computed inside `score_cohort()` from the raw block plus the cache,
by the existing `compute_coverage()`, keyed by clock id, exactly as it is today. Nothing about the
coverage machinery is chunk-aware, and nothing is patched after the fact.

**Two ordering consequences, both live.** The samples gate is warn-only and reads assembled counts,
so `check_row_coverage()` now runs *after* the scoring loop rather than before it. And asset loading
is in the data-independent tier, so `load_mc_assets()` -- which may prompt and download -- now runs
before `check_DNAm()` and the pheno checks. A caller who asks for an external pack with a malformed
matrix is prompted first and errors second. That is the price of hoisting resolution out of the
per-block path; revisit it only with a check that does not re-warn when `mc_cohort()` repeats it.

`score_cohort()` factors as `score(prep(...))`; Phase 7 (sec 9) makes that boundary public rather
than inventing one.

### 4.1 The three tiers

**Data-independent** (`mc_spec()`, once, either front end): `resolve_clocks()`,
`resolve_clocks_sequence()`, `resolve_normalize()`, `drop_routed_members()`, `resolve_DNAm_extra()`,
the covariate union, `load_mc_assets()`, `clock_panels()`, `split_cross_sample()`.

**Per-row predicates** (block-local is already correct): rownames non-NULL / non-NA, the all-NA
sample abort. The missing-pheno abort sits here too -- it is a fact about the `pheno` argument, not
about any matrix, and `mc_cohort()` raises it before reading a block.

**Cohort-set predicates** (`mc_cohort()`; each of these is *wrong* when computed from a block):

| Fact | Single pass | Chunked |
|---|---|---|
| rowname uniqueness | within one matrix | across every block |
| pheno subset | `rownames(DNAm)` vs `pheno[[pheno_id]]` | union of all block ids, so the error names every missing id once instead of block 7's |
| column classification | `n_obs == 0` / `0 < n_obs < nr`, over the matrix given | the same test on `n_obs` accumulated over the cohort |
| column means | `sum / n_obs` over the matrix given | accumulated `(sum, n_obs)` |
| clocks gate | per matrix | cohort-wide, before pass 2 |
| colname agreement | n/a | across blocks (sec 5.3) |

Both column rows read the same `col_stats()` output (sec 5.2), which is what makes the chunked
column a change of accumulator rather than a second implementation.

The samples gate is warn-only, so it runs once after assembly on the concatenated counts.

### 4.2 The invariant this buys

> Holding the cohort facts fixed, scoring any row subset yields rows **equal** to scoring the whole
> cohort.

Testable, and the only thing that actually proves chunk-safety: score a cohort, score it in three
blocks, `expect_equal`. Lives in the always-on tier as `tests/testthat/test-chunk-invariance.R`
(sec 10), and needs no chunk source -- it splits a matrix in memory and calls the seam directly.

**Measured: every clock agrees, across a 3-block split, on a cohort carrying all three missingness
shapes at once, with no exclusion.** The assertion is `expect_equal` -- agreement within tolerance,
which is the claim being made. It is **not** bit-exactness: swept over all 101 sequence entries
(2026-07-28), 94 of them differ in the last bits, worst **9.3e-10 absolute** (`DNAmB2M`, whose scale
is 6.7e6) and worst **5.1e-13 relative** (`DNAmPhysAge`). That is ~200x inside `PARITY_REL_TOL`, so
chunking is sound; the earlier "agrees to the last bit" claim here was measured over a 6-clock list,
not the pool, and is retired (DECISIONS 2026-07-28).

**The drift is not the fill, and not the mean.** `observed_panel()` returns `identical()` `cols`
and `values` for a block and for the whole cohort -- the cohort facts do their job exactly. The
reassociation is downstream, in `%*%` (`linear_score()`): reference BLAS `dgemm` blocks by `M`, so
the row count changes the accumulation order. Reproducible with no package involved -- `A[1:k, ] %*% b`
differs from `(A %*% b)[1:k]` at some `k` and not others -- and bit-exact at every `k` under
`options(matprod = "internal")`.

**Consequence for Phase 6.** Computing each mean once in `mc_cohort()` is still right, and pass 1
should accumulate in one shared accumulator so the means themselves do not move. But that cannot buy
bit-exactness for the scores: the drift above exists with the means held fixed, so an explicit bound
is the only option, not the fallback. `expect_identical` is banned package-wide regardless
(CLAUDE.md, "Test altitude").

## 5. Phase 6 -- chunk sources and the two passes (not started)

### 5.1 An iterator contract, and two adapters over it

The contract is deliberately tiny, because both passes are **sequential full scans** -- no random
access, and pass 2's block order does not matter since assembly reorders by `sample_id`. A source
must only report dimnames without loading and yield sample-blocks.

| Source | Dependency | Role |
|---|---|---|
| in-memory block iterator | none | the **always-on test vehicle** and the supported manual path. Chunk invariance is a property of the algorithm, not of I/O, so splitting a matrix in memory proves it exactly as well as reading from disk |
| `DelayedArray` | Suggests (Bioc) | **the user-facing one.** HDF5 today via whatever backend the user brings; TileDB and Zarr free later, no code change. Friction is one line (`HDF5Array("f.h5", "beta")`), not an ETL |

Dispatch on the class of one argument; no constructor function.

**Rejected: an npy directory convention.** It measured best on every axis -- 1.00x storage, 0.018s
for a 100-sample block (6x HDF5), 0.1s for pass 1 (19x HDF5), readable with base R `seek()` +
`readBin()` and no package at all, since `RcppCNPy` cannot do partial reads anyway. It loses on
maintenance: npy carries no dimnames, so it needs sidecar files, which means *we* would own a format
spec, its validator, its versioning, and its edge cases (v1 vs v2 headers, `fortran_order`, endianness,
rejecting float32 because ~1e-7 precision blows `PARITY_REL_TOL`) forever. Hand-rolled correctness is
a permanent liability that a speed win does not buy off.

**Deferred: duckdb.** Genuinely attractive -- CRAN, already in Suggests, and the wide schema
(`beta` table, `cpg` column, samples as columns) is what the parity fixtures already use, with a
working reader in `test-fixtures-parity.R`. Cohort validation as SQL is its real strength: `DISTINCT`
ids, duplicate detection, pheno anti-joins, per-CpG NA counts, all exact and out of core. It loses on
two counts. Storage: measured 1.65x in-memory against HDF5's 0.53x, i.e. **3.1x HDF5**, and `DOUBLE`
is mandatory so it cannot be tuned down. And it needs an ETL -- "duckdb reads whatever you have" is
false; the user must freeze first, which is the expensive full read+write that chunking exists to
avoid. If it is ever added, the rule is **SQL owns identity and integers, R owns every float**:
duckdb's threaded hash aggregate is not reproducible run to run (measured ~1e-15 drift across thread
counts *and* across repeat runs at default threads; `SET threads=1` is deterministic, counts and
distincts always are).

**Rejected: `bigmemory`** -- weak dimnames support is disqualifying when dimnames *are* the identity.

Zarr is not a primary target: in R it is `Rarr` (Bioc, young) and `pizzarr`, neither widely installed,
with the format's center of gravity in Python. Its performance profile is HDF5's -- both chunked and
compressed -- and its real advantages (parallel writes, object stores) do not apply to a
single-machine sequential scan. Going through `DelayedArray` means never having to choose.

### 5.2 What we own, in C++, and what DelayedArray does

DelayedArray is used for **slicing only**: `dim()`, `dimnames()`, `[` + `as.matrix()`, and
`chunkdim()` for the geometry advisory. Four generics, about as stable as Bioconductor gets --
breakage exposure scales with surface area, and this is near the minimum.

Two kernels are ours, and they live in this package's own C++ (the package is going Rcpp so that the
`betanorm` normalization functions can be vendored in; these land in the same place):

- **`fill_imp_col(obj, mean_vec)`** -- apply externally-computed per-column fills. This is the
  primitive `partial_fill` is built on and the thing that lets a block be filled from cohort
  statistics. Measured motivation: a pure-R fill is 2-5x slower than the C++ path, so the kernel is
  load-bearing, not a convenience. **It mutates `obj` in place and returns nothing** -- Rcpp's
  `NumericMatrix` wraps the caller's SEXP without copying, so whatever is passed is rewritten.
  Its one caller, `build_partial_cache()`, always hands it a fresh `DNAm[, cols, drop = FALSE]`
  slice, and matrix subsetting always allocates, so the caller's own `DNAm` is never reachable. The
  safety property is therefore a caller contract, not a guarantee of the kernel (DECISIONS
  2026-07-27, "Phase 5: the seam is cut"). **It is stated only at the call site**
  (`build_partial_cache()`, `R/missingness.R`) -- the kernel source carries no comment, so a reader
  who arrives at `src/fill_imp_col.cpp` first sees an in-place mutation with nothing saying why that
  is safe. A second caller must re-establish the fresh-slice property with nothing in the file to
  warn them.
  Its `!R_finite` test also covers `+/-Inf`, but that branch is unreachable through the supported
  path: `col_stats()` scans the same columns earlier in `scan_missing_cpgs()` and aborts on the
  first one (sec 5.4), so only NA and NaN reach the fill.
- **`col_stats(obj, cols)`** -- one traversal, returning a **list**: `stats` (a 2 x ncol matrix with named
  rows `sum` and `n_obs`), the two global range flags `any_lt0` / `any_gt1`, and `inf_at`. The mean
  is `sum/n_obs` and the NA count is `nrow - n_obs`, so one pass classifies the columns, supplies the
  fill, **and** carries the value gates (sec 5.4). This is what `scan_missing_cpgs()` now runs on,
  replacing `slideimp::mat_miss(col = TRUE)` plus `slideimp::mean_imp_col()` -- two traversals
  collapsed into one, and the classification and the fill values can no longer disagree because they
  come off the same sweep. One accumulator called by both front ends is also what makes sec 4.2's
  agreement structural rather than lucky.
  **It fails fast on `Inf`:** the scan stops at the first one and returns its 1-based `c(row, col)`
  in `inf_at`, with `stats` set to `NULL` and the range flags covering only the scanned prefix. So
  `inf_at` is checked **first** and nothing else in the list is read in that branch -- a contract the
  caller has to honor, not something the shape enforces.
  The kernel **reports and does not decide**: `check_col_values()` (R) raises the abort and the
  warnings so all three reach the user through cli. Missing and bad are counted apart -- NA/NaN are
  "not observed" and fill from the cohort mean, `+/-Inf` is bad data no fill can repair.
  **`cols` is 1-based column indices into `obj`, and `NULL` scans every column** -- so the caller
  hands the kernel the whole matrix and an index rather than a materialized slice, and
  `scan_missing_cpgs()` no longer duplicates the panel union. Positions in `stats` and in
  `inf_at[2]` are relative to `cols`, not to `obj`, which is what keeps `check_col_values()`
  indexing `present_needed` unchanged. Out-of-range, `NA`, or non-positive indices `stop()`.
  **`row_obs` serves the row half too.** `count_sample_miss()` (`R/coverage.R`) takes a clock's
  cohort-mean-filled columns as `length(cached) - row_obs`, which is per row and correct against any
  block. That retired `slideimp::mat_miss(col = FALSE)`, the package's last `slideimp` call, so the
  Import is gone. `row_obs` is `NULL` on the `Inf` bail, which cannot happen downstream of the value
  gate -- `count_sample_miss()` `stop()`s on it rather than trusting the reasoning.

`DelayedMatrixStats` is deliberately **not** a dependency: its `rowMeans2()` traverses in its own
block order and lands ~5.6e-16 from a plain `colMeans()`, which is the same non-determinism as
duckdb's threaded `avg` arriving from another direction. `setAutoBlockSize()` is a session-wide
option a package must not set for its user, so owning the slicing is also what keeps block size an
argument of ours.

Block size is exposed as a **memory budget**, not a block count, since RSS is the actual constraint:
`n_samples_per_block = budget_bytes / (8 * n_probes_needed)`.

### 5.3 Orientation and geometry -- assume neither

Bioconductor stores DNAm as **probes x samples**; this package mandates **samples x probes** with
`rownames` as identity. The two conventions are *opposite*, so assuming either silently scores
garbage when wrong.

- **Infer orientation from dimnames, then verify** against the union scoring panel -- dimnames only,
  no data read. Error rather than guess when ambiguous; take an explicit argument as an override.
  This extends existing practice: `check_DNAm()` already warns "DNAm looks transposed -- CpG ids are
  in the rows". It is not "searching" in the banned sense, which governs resolving declared catalog
  pointers, not validating user input.
- Transpose **per block** -- never `t()` the array, which materializes the thing being chunked to
  avoid.
- **Slicing along samples is fixed by fiat**, and that is a choice about our traversal, not an
  assumption about the input: pass 2 scores per sample and the fill classification is per probe
  across all samples, so sample-blocks are the only slicing serving both.
- **Chunk geometry is an advisory, never an assumption.** Measured on 20k x 1k, 100-sample block read
  vs full pass 1: `20000 x 1` (sample-major) 0.047s / 1.05s; `500 x 50` (square-ish) 0.043s / 0.62s;
  `4472 x 223` (**`writeHDF5Array` auto**) 0.391s / 0.57s; `1 x 1000` (probe-major) 0.537s / 0.85s.
  Spread is 12x -- a message, not a refusal -- and note the auto default is one of the *worse* ones
  for this access pattern, while the two passes pull in opposite directions and a square-ish chunk is
  near-best on both. Read `chunkdim()`, recommend square-ish, do not block.
- Column agreement across blocks is an **error by default**. "Absent in block 2, present in block 1"
  is not the same fact as NA, and conflating them feeds the cohort means a denominator mixing
  structural absence with measurement failure. Offer an explicit intersect opt-in, never a silent one.

### 5.4 The two passes

**Pass 1** is one fused traversal over the needed panel, accumulating per block via `col_stats()`:
per-column `(sum, n_obs)` plus the two range flags and `row_obs`, the sample ids, the within-block
row predicates, and column agreement. Then once, at the end -- this is exactly what `mc_cohort()`
already derives today from a single traversal, so Phase 6 replaces the accumulator and nothing below
it:

| derived | from |
|---|---|
| `all_na_cols` | `n_obs == 0` |
| `partial_fill` | `sum / n_obs` over columns with `0 < n_obs < n_rows`; names are the classification |
| `usable_cols` | present-needed minus `all_na_cols` |
| `cpg_list` | `resolve_cpgs(usable_cols, panels)` |
| clocks gate | `check_coverage(cpg_list, min_clocks_coverage)` -- **throws here**, before any scoring |

plus global id uniqueness and the pheno subset against the union id set. The surviving worklist lets
pass 2 read fewer columns.

**The value gates already ride this sweep** (built ahead of Phase 6, because a latent bug forced it
-- see below). `check_col_values()` reads `inf_at` and the two range flags off `col_stats()`: any
`+/-Inf` **stops**, naming the sample and CpG it sits at; `any_lt0` and `any_gt1` each **warn**,
separately, because a matrix can trip both and neither implies the other. The range check is the
beta-vs-M-value tell -- `check_DNAm()` validates `mode = "double"` but not the range, so a minfi
`GenomicRatioSet` whose assay holds M-values would otherwise score plausible garbage. Exactly 0 and
exactly 1 pass, since those are ordinary saturated betas. It warns rather than stops because the
range is strong evidence, not a definition.

Deliberate asymmetry in what each gate can say: the abort carries an exact position, the warnings
carry none. An `Inf` is a single point defect worth locating, and the scan is already stopping
there; out-of-range values are a property of the whole matrix, where naming ten of 866k columns
would be noise. Two global booleans is the whole cost.

**An all-missing column is classified, not rejected.** `n_obs == 0` means no cohort mean exists,
which is an ordinary expected case: the column lands in `all_na_cols`, leaves `usable_cols`, and
takes the clock's vendored ref or the drop policy -- the standing imputation invariant. Only
`partial` columns are divided, and `n_obs > 0` holds for all of them by construction, so `0/0` never
arises. Erroring on it instead would refuse every matrix missing one probe.

**Why it could not wait for Phase 6.** `scan_missing_cpgs()` short-circuited on `anyNA(DNAm)`, and
`anyNA()` does not see an `Inf` -- but `col_stats()` treats every non-finite value as missing. So
whether an `Inf` was filled depended on whether an *unrelated* NA existed elsewhere in the matrix:
the same `Inf` gave a `-Inf` score with no NA present and a plausible 60.0 with one. (The
pre-Phase-5 `slideimp::mat_miss()` counted only NA, so it was consistently the former.) The scan now
runs unconditionally over the needed panel, which is what makes the gates reachable at all, and no
`anyNA()` short-circuit survives anywhere in `R/`.

**The dead-row half is scoped to the scoring panels, and that is a separate predicate.** The sweep
covers `panels_union()` -- score **and** norm -- because both halves need classification and fill
values, but a sample is dead iff it observes nothing on the **scoring** panels: a normalization CpG
cannot score anyone. Judging it on the union let a DunedinPACE sample with 0/173 scoring CpGs and
all 19,827 background CpGs pass, then score off a fully vendor-filled panel. `scan_missing_cpgs()`
now takes the scoring union as its own argument and reuses the union `row_obs` only when the two
sets coincide, which is every request with no normalizing clock (DECISIONS 2026-07-28).

**Pass 2**, per block: re-read **raw**, build the block's cache with
`fill_imp_col(block[, names(partial_fill)], values = partial_fill)`, call `score_cohort()` with the
block and the facts. Coverage for the block comes out of the unchanged `compute_coverage()`, which
counts NAs in the raw block over the cached columns. Keep the score matrices and the coverage
fragment.

**Assembly** is a concatenate and a sum, never a replace:

- score matrices and `pending` intermediates -- concatenate by clock id (disjoint rows)
- `sample_miss$score` / `$norm` -- concatenate the per-block vectors (disjoint rows)
- `score_imputed_partial` / `norm_imputed_partial` -- sum across blocks
- every other `per_clock` field is derived from the shared `cpg_list`, so every block already
  computed the identical value; take it

Then reorder by `sample_id`, call `finalize_cross_sample()` on the assembled intermediates (sec 6.2),
and run the samples gate once.

This is **not** `rbind.mc_result` and must not become it: assembly works on fragments before any
record exists, sums coverage counts, and finalizes a reduction that has not run yet, while `rbind`
binds finished records whose counts must stay per-record and whose reductions already ran (sec 8).

Naming: `calc_clocks_chunked()` follows the existing family. The wrapper stores a chunk label per
sample (sec 8).

### 5.5 Column narrowing conflicts with Zhang2019

The obvious RSS win is reading only the union panel, but a `sample_scale` clock takes its moments
over **every** probe in the input. Reuse the existing declared-recipe predicate `needs_full_panel()`
rather than a clock list -- it currently lives in `test-fixtures-parity.R` and moves to `R/` as part
of this phase. Whole-array column means cost ~7 MB even at 866k probes, so filling full width when
that predicate is true is affordable.

## 6. Phases 2 and 3 -- cohort reductions leave the scoring loop

### 6.1 Phase 2 -- the split, derived (built)

Exactly **2 of 129** catalog clocks are cross-sample: `DNAmPhysAge` and `DNAmPhysAge_years`
(`cross_sample_at = 11`). The fill does not touch them.

`clock_cross_sample_at()` / `clock_is_cross_sample()` read the declared field (`R/accessors.R`),
`split_cross_sample()` partitions a compute sequence, and `mc_spec()` stores the cohort-reducing half
as `spec$cross_sample`. Nothing carries a clock list: the census test walks `resolve_clocks("all")`
and asserts the partition is total, so a clock that starts declaring `cohort_zscore` upstream
classifies with no downstream edit.

**One sync defect fixed on the way.** `attach_sex_routed_aliases()` minted every alias with a
hardcoded `cross_sample_at = NA_integer_` while deriving every *other* alias field from its members.
An alias scores each row with exactly one member, so a cross-sample member makes the alias
cohort-dependent -- and the hardcoded NA would have told a chunked engine the alias was per-sample,
producing wrong scores rather than an error. It is now `alias_cross_sample_at(members)` (the min,
NA when neither reduces). Latent today: all 14 routed members are `DNAmFitAge` and per-sample, so
the derived value equals the old constant for all 7 aliases and `R/sysdata.rda` is unchanged
(verified, not assumed). The guard against it recurring is the census test, which checks each alias
against its members rather than against a constant.

The rest of the sync-side classification was audited at the same time and is correct: `cohort_zscore`
is the only cross-sample op in the whole upstream corpus (3 occurrences, 2 clocks); `sample_scale`
(Zhang2019) is a within-sample z-score; `center_scale` (SystemsAge) reads the vendored
`systems_pca_center` / `systems_pca_scale` tensors, so it is per-sample too.

### 6.2 Phase 3 -- the reduction moves after assembly (built)

`score_PhysAge()` reduces over an `n x k_surrogate` matrix of **derived surrogate scores**, not over
betas, plus an `n >= 2` guard. The two clocks do not have the same number of reductions, and Phase 3
has to replay each one exactly:

| clock | recipe | reductions |
|---|---|---|
| `DNAmPhysAge` | `cohort_zscore` (11) -> `row_sum` (12) -> `transform` (13) | one: `scale(raws)` |
| `DNAmPhysAge_years` | `cohort_zscore` (11) -> `row_sum` (12) -> `cohort_zscore` (13) -> `poly` (14) | two: `scale(raws)`, then `scale(phys)` before the polynomial |

`cross_sample_at` records only the **first** reduction, which is all a streaming engine needs (it
says where the chunk-safe prefix ends). It is not a step pointer for the finalize, which has to
follow the branch.

So the chunk-safe output is the surrogate raws (8 columns, ~640 KB at n = 10,000) and the reduction
is a cheap post-assembly step. It is required on any route: after a per-block run the surrogate raws
are gone and only the polynomial output survives, which cannot be inverted, so a per-block
`DNAmPhysAge` can never be repaired at assembly. Consequence accepted: under chunking `DNAmPhysAge`
is not a finished column until assembly, so the wrapper's last step is a **finalize**, not a bare
bind.

**What landed.** `score_PhysAge()` split into `physage_raws()` (per-sample, in the loop) and
`finalize_PhysAge()` (the reduction, following the branch rather than a step index, since the two
clocks reduce a different number of times). `score_cohort()` returns a third element, `pending`;
the loop routes a clock's output there instead of into `scores` **iff it is in `spec$cross_sample`**,
so the defer decision is the catalog's declaration and no clock list exists anywhere. The `n >= 2`
guard moved to the finalize, which is the correct altitude -- a 1-row *block* is fine, a 1-row
*cohort* is not, and both are now enforced where they are true.

**Sufficient statistics were considered and rejected as unnecessary.** Accumulating per-surrogate
`(sum, sumsq, n)` per block would let the z-score stream, but the raws are 8 doubles per sample, so
carrying them costs nothing and avoids a second accumulator to keep bit-agreeing with the first. It
also would not have covered `DNAmPhysAge_years`, whose *second* `cohort_zscore` reduces over `phys`,
itself a function of the first -- not streamable at all, trivial with the raws in hand.

**Two properties measured:** a 3-block split now finalizes to the single-pass answer for both clocks
(from 2.8 and 11 **years** off, to within the residual float drift of sec 4.2 -- `DNAmPhysAge` is
its worst relative case at 5.1e-13), and single-pass output is `identical()` to the pre-Phase-3
implementation -- so the standing parity run is undisturbed.

Both front ends call `finalize_cross_sample()` **unconditionally**; `pending` is empty for the other
127 clocks, so it is a no-op for a typical request. Nothing anywhere tests "is this PhysAge" or "are
we chunked".

## 7. What `augment()` owns (not started)

Cross-sample derivations that are not a clock's definition -- age acceleration / residuals,
user-requested z-scores -- happen in `augment()`, after binding, never in the scoring loop. That is
what dissolves the third `cbind` gate. `augment()` does not exist yet and is an **unbuilt idea, not
a contract** -- adding it is a new API decision (CLAUDE.md; DECISIONS 2026-07-27).

## 8. Phase 4 -- `rbind` (not started)

`rbind.mc_result` refuses today (`R/mc_result.R`). It can admit records under gates:

**Phase 3 did not make a score column batch-independent, and the earlier wording here claimed it
did.** It relocated the reduction from inside the scoring loop to after assembly *within one run*,
which is what buys chunk-safety. A finished record's `DNAmPhysAge` is still a z-score against that
run's cohort, so binding three time-point records yields three within-batch z-scores in one column.
That is the residual asymmetry stated at the end of this section, and the batch label is the answer
to it -- **not** a re-finalize at bind time. Re-running the reduction over the union would silently
rewrite numbers the user has already seen, and would standardize across batches that were
deliberately imputed differently (sec 0, sec 11).

1. **Disjoint** `sample_id` sets (`cbind` requires equal; this requires disjoint). Collision throws
   and names the ids -- this also catches lazy per-block `sample1..N` labelling, which is what the
   retired `synthetic_ids` gate was for (DECISIONS 2026-07-24).
2. Identical clock column sets; reorder to the first record, or throw.
3. Coverage denominators comparable -- same panels per clock, which follows from one catalog.
4. pheno consistency: same `pheno_id`, and no id appearing twice with different covariates (moot
   under gate 1, but the frames still have to concatenate).

**Record, never refuse, on differing fill regimes.** The time-point user is binding
differently-imputed batches on purpose. What makes that honest is a per-sample **batch label** plus a
per-batch imputation summary, so a saved record can answer which fill values produced any given row.
Chunk reassembly labels every row one batch; the time-point union labels three. Same shape, different
scientific object, and the label is the only thing that says which.

Phase 4 has a prerequisite this doc does not own: the verb surface. `[`, `cbind`,
`as.data.frame.mc_result` and `augment()` do not exist yet, so the batch-label-survives-subset
behavior cannot be built or tested until they do.

One residual asymmetry worth stating plainly: a record bound from per-time-point runs carries a
`DNAmPhysAge` that was z-scored **within** batch. That column is meaningful but is not the column a
single-cohort run would produce. The batch label is what lets a reader see that; do not try to
recompute it at bind time.

## 9. Phase 7 -- `prep()`, the prepared-input record (not started)

`prep(DNAm, clocks, ...)` returns an S3 record, list-like and keyed by clock id, whose element for a
clock is **the matrix that clock's coefficients multiply**: its panel, imputed clock-wise and
panel-wise, with its normalization applied or not. `calc_clocks()` accepts either a `DNAm` or a
`prep` record, because `score_cohort()` already factors as `score(prep(...))` -- Phase 7 names an
existing boundary and hands the user the object at it.

### 9.1 It is a plan, not a payload

Materialization is **lazy by default**: `p[["Horvath1"]]` builds one clock's matrix on demand and
`calc_clocks(p)` iterates with one resident at a time, so peak RSS is unchanged from today.
`materialize(p)` is the eager form, for callers who want the whole thing and can afford it.

Measured, and the reason laziness is not optional -- `Sigma |panel_i|` against the union:

| | `Sigma \|panel\|` | union | overlap | eager, n=1000 |
|---|---|---|---|---|
| 101 bundled clocks | 53,127 | 39,025 | 1.36x | 425 MB |
| 28 external (pack) clocks | 3,083,623 | 378,363 | **8.15x** | **24.7 GB** |

PCBrainAge alone declares 357,852 CpGs, the 13 SystemsAge organ clocks 125,175 each, the 14 PCClocks
78,464 each. A full EPIC `DNAm` at n=1000 is ~6.9 GB, so eager prep over everything is 3.6x the input
matrix. Eager is free for the bundled set and fatal for the packs, so the record holds a plan --
`DNAm` reference, the shared cohort-mean cache, the per-clock panel spec, the per-clock normalization
decision -- and materializes from it.

The shared cache survives this: the partial-NA cohort-mean fill is clock-independent, so it is built
once and every materialization draws from it. Only the absent-CpG vendored fill and the normalization
are per clock, and those happen at materialize time.

### 9.2 Why there is no single prepared matrix

`prep()` deliberately does **not** return one imputed-and-normalized `DNAm`, because no such object
is well defined:

- **Imputation is per clock.** A fully absent CpG has no single fill value -- it takes the clock's
  own vendored ref as a scoring CpG and the gold-standard mean as a normalization CpG
  (`score_Dunedin.R`, the `fill_ref` split). The same probe gets a different number depending on
  which clock is asking.
- **Normalization is per clock.** Different clocks declare different normalization panels and
  different targets, so one normalized matrix cannot serve two of them, and normalizing the shared
  input would silently change what every non-normalizing clock is scored on.

The line that does exist is cohort-global vs clock-specific: cohort means are hoisted (sec 1.1),
everything else stays per clock. `prep()` is a list *because* of that, not in spite of it.

### 9.3 Normalization is an argument, not part of this record

Normalization already lives inside the clock's branch, on the clock's own declared panel and target,
and never touches the shared `DNAm` -- `score_Dunedin()` builds a local matrix over the declared
`quantile_normalization_background` panel, normalizes it against the declared gold-standard tensor,
and takes only its scoring CpGs out of the result. There is nothing global to turn off, so the toggle
is a per-clock decision -- but it is **not** `prep()`'s to own. `normalize=` is an ordinary named
logical argument on `calc_clocks()` / `sim_DNAm()`, resolved once by `resolve_normalize()` before any
DNAm is read; `prep()` takes the same argument when it lands (DECISIONS 2026-07-27, "`normalize=` is
an argument, not `prep()`").

Current declared state: `normalization` is non-`none` on 4 clocks -- `BMIQ` on Horvath1 and Knight,
`noob` on Horvath2, `quantile` on DunedinPACE. Both BMIQ clocks ship a `bmiq_gold_standard`
probe_set with a file pointer: `goldstandard2.csv.gz`, 21,368 probes, `identical()` between the two
groups, with each scoring panel a strict subset (353/353, 148/148). `noob` is moot -- it is an
IDAT-intensity background correction, always done by the pipeline that generated the betas, and
unreachable from a beta matrix under any design.

Both facts the design needed are now settled:

1. **Expressible as a declared panel plus a vendored target?** Yes for both schemes, so the existing
   machinery generalizes -- `clock_norm_target()` is the scheme-agnostic target accessor and the
   branch becomes normalize-then-linear dispatching on `clock_norm_scheme()` (`quantile` ->
   `quantile_norm`, `bmiq` -> `bmiq_calibration`), the existing **pre-transform** branch kind and not
   a new entry in the closed branch set. A genuine whole-array BMIQ needing type I/II annotation
   would not be expressible in the catalog at all and stays the caller's own preprocessing -- which
   is exactly why upstream moved EpiTOC/EpiTOC2/HypoClock/Mayne/PedBE to `none`.
2. **Constitutive or published preprocessing?** Read off the scheme, not a per-clock field, which
   upstream did not ship: `quantile` is constitutive (DunedinPACE's QN is part of the definition and
   cannot be declined), `bmiq` is preprocessing and defaults **off**. A closed enum value, not a
   clock-id exception, so a future BMIQ clock inherits the answer.

`prep()` also takes the caller's statement about the input (already BMIQ-calibrated upstream, say),
which makes a clock skip a transform it would otherwise apply rather than double-calibrating. It
cannot police *which* gold standard the caller used; that is documentation.

Consequence, live once a caller opts in: Horvath1 goes from needing 353 CpGs to needing 353 to score
and ~21,000 to normalize -- the DunedinPACE shape (173 / 20,000), already handled by the norm-panel
coverage axis and by the gate split (scoring panel stops, normalization panel warns). Because the
default is off, `norm_needed` is empty for both BMIQ clocks unless asked, so nothing about the
required panel or the coverage surface moves for a caller who does not opt in. CLAUDE.md's
"(only DunedinPACE)" parenthetical about norm panels holds by default and goes stale on request.

Where BMIQ departs from QN: a fully absent norm CpG is **dropped**, not vendor-filled. QN fills from
the gold mean and normalizes the full panel; BMIQ estimates the sample's own mixture from the panel,
so filling with target-drawn values biases that fit toward the gold standard and shrinks the very
correction being computed. The thin-background warning's "filled from the reference mean" text is
therefore wrong for a BMIQ clock and has to become scheme-aware.

Normalization costs nothing in chunk-safety: `quantile_norm()` is per-sample against a fixed target
and is bit-identical on a row subset (measured, max abs diff 0), and BMIQ fits its mixture within a
sample. Neither adds to `cross_sample_at`.

### 9.4 Open

- What `p[["Horvath1"]]` returns for a normalizing clock: the normalized **scoring** panel (the thing
  the coefficients multiply), with the normalization panel transient inside materialization. The
  alternative -- returning the 21k panel -- makes element widths mean different things per clock and
  blows `materialize()` up by 60x on one clock. The normalization panel is what an auditor wants to
  see, so it likely deserves a separate accessor rather than being the default element.
- ~~The `normalize =` grammar.~~ Settled and shipped: a bare scalar is a policy that reaches only
  the clocks able to honor it, a named logical is a per-clock claim whose unhonorable entries error.
  No per-clock options list -- every BMIQ setting is fixed at Horvath's published choice.
- Deduplicating normalization work. BMIQ is a per-sample mixture fit and far more expensive than
  rank mapping; if several clocks share one gold panel and target, the normalized result should be
  memoized by (panel, target) identity. `resolve_cpgs()` already dedups the *panels* into
  `norm_parts` with a per-clock index -- the missing piece is caching the normalized output in the
  same place. Horvath1 and Knight are the live case: same 21,368 panel, `identical()` targets, so a
  run normalizing both currently fits the mixture twice per sample. Confirmed safe to memoize, not
  yet done. It matters more under chunking, where it would otherwise repeat every block.
  Sync-side, the same fact means the two vendored `goldstandard2.csv.gz` copies can collapse to one
  bundled tensor (hash both, **assert equality**, keep each clock resolving its own declared
  pointer) -- optional, and only worth doing next time `sync.R` is open.

**Row-level failure has no vocabulary yet.** BMIQ runs `on.sample.error = "continue"` with
`failed.sample = "NA"`, so a sample whose mixture fit fails scores `NA` while the cohort completes.
Measured consequence: that sample shows `coverage = 1.00` on **both** panels in `samples_coverage()`
and a zero `sample_miss`, because coverage is computed upstream of scoring and is correctly reporting
that every probe was present. Nothing in a saved record says why the row is `NA` -- missing pheno,
excessive NA and a failed calibration are indistinguishable after the fact, and the runtime warning
that does fire names no samples. `betanorm::bmiq_calibration()` returns `$success` (per sample) and
`$failures` (`sample_index`, `sample_name`, `stage`, `message`); `bmiq_panel()` currently keeps only
`$calibrated` and discards both. A `samples_coverage()` reason column is the intended fix.

**Not a gap:** BMIQ failing on `sim_DNAm()` output is expected, not a bug. `random_betas()` is
U(0,1), which has no 3-state structure, so the fitted components overlap and
`density_thresholds()` finds no crossing between adjacent means. Real methylation is multi-modal
and does not hit this. Tests that need a normalizing clock must build methylation-shaped input --
`test-normalize.R` jitters the gold standard, which is reliable where `random_betas()` is not.

## 10. Tests

Always-on, no new dependency -- the in-memory block iterator (sec 5.1) is what makes this possible,
and it is the whole reason that source exists. Every numeric comparison is `expect_equal`;
`expect_identical` is banned package-wide.

Landed with Phase 5, in `test-chunk-invariance.R`. One helper builds a cohort carrying all three
missingness shapes at once (ordinary partial NA, partial-but-all-NA-in-block-1, fully absent) and
every test reuses it, so the shapes cannot drift apart between assertions:

- **Chunk invariance** (sec 4.2): one cohort scored whole vs in three blocks, through the seam
  directly -- no chunk source needed. Both paths end in the same `finalize_cross_sample()`, and
  since Phase 3 the comparison runs over **every** clock in the sequence with no exclusion. The
  request is asserted to contain a cohort-reducing clock (`length(spec$cross_sample) > 0`), derived
  rather than named, so the test cannot silently stop exercising the deferred path.
- **The loop defers exactly the declared set**: `names(pending)` equals `spec$cross_sample`, those
  clocks leave the loop with a wider-than-one-column per-sample intermediate and no score, and the
  finalize turns each into an ordinary `n x 1` column.
- **The sample-axis split is derived**: `split_cross_sample()` partitions `resolve_clocks("all")`
  totally, and each sex-routed alias matches its members rather than a constant (the standing guard
  on the sync fix in sec 6).
- **The block-dependency case**: a column partial cohort-wide but all-NA within one block lands in
  `partial_fill` and counts as `imputed_partial`, not `imputed_full`, with every row of that block
  counted as imputed for it.
- **Coverage assembles by concatenate and sum**: `score_imputed_partial` sums across blocks,
  `sample_miss` concatenates, and panel-derived fields are already identical per block.
- **The clocks gate throws out of `mc_cohort()`**, before anything is scored.

Still to write, with Phase 6:

- **Gates fire before pass 2**: assert via a source that records how many times each block was read.
- **`rbind`**: disjoint ids bind; colliding ids throw; batch labels survive.

Landed early, in `test-value-gates.R` (sec 5.4): `+/-Inf` stops with and without an unrelated NA
present -- the two cases that once disagreed -- and the abort names the sample and CpG it sits at;
each range flag warns on its own and an M-value matrix trips both while still returning scores;
ordinary betas and all-NA columns pass in silence; and the kernel bails on `Inf` with `stats` unset,
which is the contract the caller depends on.

Skipped when Suggests are absent, per existing practice: the `DelayedArray` source, including one
transposed input and one deliberately hostile chunk geometry.

Parity is untouched -- it scores complete cohorts and is not part of "run the tests".

For Phase 7: `calc_clocks(DNAm, ...)` and `calc_clocks(prep(DNAm, ...))` agree; a lazy `prep` and a
`materialize()`d one agree; and a clock whose normalization is declined scores what the same clock
scores today.

## 11. Non-goals

- No random-access source contract. Two sequential scans is the whole requirement.
- No pre-imputation of a block outside the scorer (sec 1.1).
- No single prepared `DNAm` out of `prep()` -- neither the imputed nor the normalized form is a
  well-defined object (sec 9.2).
- No eager materialization by default, and no whole-array normalization inside `calc_clocks()`.
- No chunk-size autotuning. Take the block size from the caller (or the source's own geometry) and
  say what was used.
- No parallelism in the first cut. Chunking trades speed for RSS by design; adding cores changes the
  peak-memory story it exists to fix.
- No `rbind` across cohorts *reconciling* anything -- no re-imputation, no re-z-scoring, no merged
  coverage denominators. It binds and labels.
- No new identity key. Ids are the identity; `batch_set_id` was removed for good reasons
  (DECISIONS 2026-07-24) and a batch **label** is not its return.
