# Chunked scoring, prepared inputs, and record binding

Tracked design doc. Cited by `dev/detail-plan.md` sec 5.1 / 6 / 7.1 and by four `dev/DECISIONS.md`
entries, **by phase number** -- see "Phase identifiers are stable" below before renumbering anything.
Section numbers are cited by nothing and may move freely.

Three problems share one seam (Phase 5) and are therefore designed together: chunking, binding, and
the prepared-input record `prep()` (Phase 7). Per-clock normalization was originally folded into
Phase 7 and has since been split out as a plain `normalize=` argument that shipped on its own
(sec 9.3, DECISIONS 2026-07-27 latest).

Status: **no chunking is built.** Phase 1 landed by other means; Phase 2 landed halfway; Phase 7's
`normalize=` half shipped in full (sec 9.3). Phases 3, 4, 5, 6 and the `prep()` record itself are
still design only.

**Picking this up cold?** Read sec 3 (phase table) and sec 4 (the seam) first -- Phase 5 gates
everything except 7. The normalization work is done and is *not* a prerequisite for any of it; it
shipped on the same seam only by coincidence of ownership. Sec 9.4 lists what is left near `prep()`.

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

Both front ends call a shared **scoring internal**, and the cohort facts reach it as arguments --
chiefly `partial_fill`, a named numeric vector whose **names are the cohort-partial columns** and
whose **values are the cohort means**. Public `calc_clocks()` calls the internal with
`partial_fill = NULL`, meaning "derive it from the matrix you were handed", which is what makes the
precondition above true of the public function and false of nothing else.

`partial_fill` does two jobs, and the second is the load-bearing one:

1. it supplies the fill values;
2. **its names are the column classification**, which is a cohort fact and must not be re-derived
   from a block.

So when `partial_fill` is supplied, the column half of `scan_missing_cpgs()` is **bypassed, not
overridden** -- there is one source of truth, never a local computation that a caller corrects
afterwards. Its row half (the all-NA-sample abort, `row_miss == ncol`) is per-row, is correct against
a block, and stays.

Per needed column, with `partial_fill` supplied:

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
| 2 | kind-1 / kind-2 split (chunk-safe vs cohort-reducing clocks) | **half done** -- `cross_sample_at` is in the catalog; nothing acts on it |
| 3 | every cross-sample op leaves the scoring loop | not started |
| 4 | `rbind` gates | not started |
| 5 | the seam: split `calc_clocks()` so both front ends compose the same internals | not started |
| 6 | chunk sources + the two passes | not started |
| 7 | `prep()` -- the prepared-input record | not started; its `normalize=` half split out and shipped |

Execution order: **5 -> 2 -> 3 -> 6 -> 4**, with 7 available any time after 5. Phase 7 is orthogonal
to chunking but lands on the same seam, so it is designed here rather than in a doc of its own.

## 4. Phase 5 -- the seam

Split `calc_clocks()` into three pieces, divided by what each depends on, so both front ends compose
the same internals:

```
mc_spec(clocks, pheno_id, from, ask)        # data-independent
mc_cohort(DNAm_or_source, spec, pheno, ...) # one scan -> the cohort facts; gates
score_cohort(DNAm, spec, facts)             # a matrix and the facts -> score matrices
```

`calc_clocks()` is those three plus `construct_mc_result()`; `calc_clocks_chunked()` is the same
three with `mc_cohort()` accumulating over blocks and `score_cohort()` called per block.
`score_cohort()` never learns which front end called it and `facts` is the only channel between them,
which is what makes sec 4.2 a property rather than a coincidence. The wrapper calls the internal,
**not** the public `calc_clocks()` -- otherwise it pays per-block resolution and asset loading and
emits per-block error messages.

**Coverage does not move.** It is computed inside `score_cohort()` from the raw block plus the cache,
by the existing `compute_coverage()`, keyed by clock id, exactly as it is today. Nothing about the
coverage machinery is chunk-aware, and nothing is patched after the fact.

`score_cohort()` factors as `score(prep(...))`; Phase 7 (sec 9) makes that boundary public rather
than inventing one.

### 4.1 The three tiers

**Data-independent** (once, either front end): `resolve_clocks()`, `resolve_clocks_sequence()`,
`drop_routed_members()`, `resolve_DNAm_extra()`, the covariate union and its missing-pheno abort,
`load_mc_assets()`, `clock_panels()`.

**Per-row predicates** (block-local is already correct): rownames non-NULL / non-NA, the all-NA
sample abort.

**Cohort-set predicates** (`mc_cohort()`; each of these is *wrong* when computed from a block):

| Fact | Single pass | Chunked |
|---|---|---|
| rowname uniqueness | within one matrix | across every block |
| pheno subset | `rownames(DNAm)` vs `pheno[[pheno_id]]` | union of all block ids, so the error names every missing id once instead of block 7's |
| column classification | `col_miss == nr`, block `nr` | accumulated `n_non_na` over the cohort |
| column means | over the matrix given | accumulated `(sum, n_non_na)` |
| clocks gate | per matrix | cohort-wide, before pass 2 |
| colname agreement | n/a | across blocks (sec 5.3) |

The samples gate is warn-only, so it runs once after assembly on the concatenated counts.

### 4.2 The invariant this buys

> Holding the cohort facts fixed, scoring any row subset yields rows **equal** to scoring the whole
> cohort.

Testable, and the only thing that actually proves chunk-safety: score a cohort, score it in three
blocks, `expect_equal`. Belongs in the always-on tier (sec 10).

**Equal, not identical, and the reason is the mean.** A single pass sums a column's non-NA values in
one traversal; a chunked run sums block by block and adds the partials. Different association order,
so the two means can differ in the last bit, and every score built on them moves with it. Either the
single-pass path accumulates in the same block order -- one shared accumulator, agreement by
construction -- or the test carries an explicit bound. Prefer the shared accumulator and keep the
bound tight; `expect_identical` is banned package-wide regardless (CLAUDE.md, "Test altitude").

## 5. Phase 6 -- chunk sources and the two passes

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

- **`fill_imp_col(obj, values, subset = NULL)`** -- apply externally-computed per-column fills.
  This is the primitive `partial_fill` is built on and the thing that lets a block be filled from
  cohort statistics. It must **slice and fill in one traversal, returning a new matrix**: the slice
  is already a fresh allocation, so there is nothing to gain from mutating in place and everything to
  lose, since an in-place kernel that ever sees the caller's `DNAm` corrupts it invisibly.
  Measured motivation: a pure-R fill is 2-5x slower than the existing C++ `mean_imp_col()` path, so
  the kernel is load-bearing, not a convenience.
- **`col_stats(obj)`** -- per-column `(sum, n_non_na, n_na, min, max)` in **one** traversal. Pass 1
  needs about six quantities and each separate sweep is another read of the disk; a fused loop
  measured 1.31s against 2.33s for two `DelayedMatrixStats` calls. `min`/`max` come free and give the
  beta-vs-M-value range check (sec 5.4). One accumulator called by both front ends is also what makes
  sec 4.2's agreement structural rather than lucky.

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
per-column `(sum, n_non_na, n_na, min, max)`, the sample ids, the within-block row predicates, and
column agreement. Then once, at the end:

| derived | from |
|---|---|
| `all_na_cols` | `n_non_na == 0` |
| `partial_fill` | `sum / n_non_na` over columns with `0 < n_na`; names are the classification |
| `usable_cols` | present-needed minus `all_na_cols` |
| `cpg_list` | `resolve_cpgs(usable_cols, panels)` |
| clocks gate | `check_coverage(cpg_list, min_clocks_coverage)` -- **throws here**, before any scoring |

plus global id uniqueness and the pheno subset against the union id set. The value range from
`min`/`max` is a free beta-vs-M-value check: `check_DNAm()` validates `mode = "double"` but not the
[0, 1] range, so a minfi `GenomicRatioSet` whose assay holds M-values currently scores plausible
garbage. The surviving worklist lets pass 2 read fewer columns.

**Pass 2**, per block: re-read **raw**, build the block's cache with
`fill_imp_col(block[, names(partial_fill)], values = partial_fill)`, call `score_cohort()` with the
block and the facts. Coverage for the block comes out of the unchanged `compute_coverage()`, which
counts NAs in the raw block over the cached columns. Keep the score matrices and the coverage
fragment.

**Assembly** is a concatenate and a sum, never a replace:

- `sample_miss$score` / `$norm` -- concatenate the per-block vectors (disjoint rows)
- `score_imputed_partial` / `norm_imputed_partial` -- sum across blocks
- every other `per_clock` field is derived from the shared `cpg_list`, so every block already
  computed the identical value; take it

Then reorder by `sample_id`, run the samples gate once, and finalize any cohort reductions (sec 6).

Naming: `calc_clocks_chunked()` follows the existing family. The wrapper stores a chunk label per
sample (sec 8).

### 5.5 Column narrowing conflicts with Zhang2019

The obvious RSS win is reading only the union panel, but a `sample_scale` clock takes its moments
over **every** probe in the input. Reuse the existing declared-recipe predicate `needs_full_panel()`
rather than a clock list -- it currently lives in `test-fixtures-parity.R` and moves to `R/` as part
of this phase. Whole-array column means cost ~7 MB even at 866k probes, so filling full width when
that predicate is true is affordable.

## 6. Phases 2 and 3 -- cohort reductions leave the scoring loop

Exactly **2 of 129** catalog clocks are cross-sample: `DNAmPhysAge` and `DNAmPhysAge_years`
(`cross_sample_at = 11`). The fill does not touch them.

`score_PhysAge()` reduces over an `n x k_surrogate` matrix of **derived surrogate scores**, not over
betas -- `scale(raws)`, then `rowSums`, then a second `scale(phys)` before the polynomial. Two nested
cohort z-scores, plus an `n >= 2` guard.

So the chunk-safe output is the surrogate raws (negligible memory) and the reduction is a cheap
post-assembly step. Phase 2 is the engine acting on `cross_sample_at`; Phase 3 emits per-sample
intermediates from the scoring loop and finalizes after assembly. For a complete cohort
`calc_clocks()` composes both and its output is unchanged.

This is required on any route: after a per-block run the surrogate raws are gone and only the
polynomial output survives, which cannot be inverted, so a per-block `DNAmPhysAge` can never be
repaired at assembly. Consequence to accept: under chunking `DNAmPhysAge` is not a finished column
until assembly, so the wrapper's last step is a **finalize**, not a bare bind.

## 7. What `augment()` owns

Cross-sample derivations that are not a clock's definition -- age acceleration / residuals,
user-requested z-scores -- happen in `augment()`, after binding, never in the scoring loop. This is
what `dev/detail-plan.md` sec 7.1 means by the third `cbind` gate being "dissolved". `augment()` does
not exist yet; CLAUDE.md listing it among the verbs is stated intent (DECISIONS 2026-07-24).

## 8. Phase 4 -- `rbind`

`rbind.mc_result` refuses today (`R/generics.R`), and that is correct while any score column can be
cohort-dependent. After Phase 3 no score column is, and it can admit records under gates:

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

## 9. Phase 7 -- `prep()`, the prepared-input record

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

### 9.3 Normalization is decided here, from declarations

Normalization already lives inside the clock's branch, on the clock's own declared panel and target,
and never touches the shared `DNAm` -- `score_Dunedin()` builds a local matrix over the declared
`quantile_normalization_background` panel, normalizes it against the declared gold-standard tensor,
and takes only its scoring CpGs out of the result. There is nothing global to turn off, so the toggle
is a per-clock decision -- but it is **not** `prep()`'s to own. `normalize=` is an ordinary named
logical argument on `calc_clocks()` / `sim_DNAm()`, resolved once by `resolve_normalize()` before any
DNAm is read; `prep()` takes the same argument when it lands (DECISIONS 2026-07-27 latest).

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

- **Chunk invariance** (sec 4.2): one cohort scored whole vs in three blocks. Include a clock with
  partial NAs, one with a fully absent probe, and `DNAmPhysAge`.
- **The block-dependency case**: a column partial cohort-wide but all-NA within one block must land
  as `imputed_partial`, not `imputed_full`, in the chunked record's coverage -- and every row of that
  block must be counted as imputed for it.
- **Coverage is real under chunking**: `score_imputed_partial` and `sample_miss` from a chunked run
  equal the single-pass values.
- **Gates fire before scoring**: a cohort failing `min_clocks_coverage` throws without reading
  pass 2 (assert via a source that records how many times each block was read).
- **`rbind`**: disjoint ids bind; colliding ids throw; batch labels survive.
- **Value range**: an M-value matrix (roughly [-10, 10]) is rejected by pass 1, not scored.

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
