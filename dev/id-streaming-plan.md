# Chunked scoring and record binding

Tracked design doc. Cited by `dev/detail-plan.md` sec 5.1 / 6 / 7.1 and by four `dev/DECISIONS.md`
entries, **by phase number** -- see "Phase identifiers are stable" below before renumbering anything.

Status: **nothing here is built.** Phase 1 landed by other means; Phase 2 landed halfway. The rest is
design.

---

## 0. The two problems are not one problem

| | Chunking | Binding |
|---|---|---|
| Goal | cap resident set size | reuse the `mc_result` verbs over records scored separately |
| Speed | may be sacrificed (may even win near the thrash boundary) | irrelevant |
| Fill values | **one** cohort mean set, shared by every chunk | **one set per record**, deliberately different |
| Result | bit-identical to a single-pass run | a labelled union of differently-imputed batches |

Binding is required *by* chunking but is strictly wider than it. The motivating non-chunk case: three
time points, not nested, imputed per time point on purpose, then bound to use `[`, `cbind`,
`augment`, `clocks_coverage`. Both cases produce disjoint id sets and identical clock columns, so
**the id gates cannot tell them apart** -- only recorded batch provenance can (sec 9).

## 1. The load-bearing rule: `calc_clocks()` takes a complete cohort

`calc_clocks()` assumes its `DNAm` is the whole cohort, always. Cohort-completeness is a
**precondition, not a parameter**: no `chunked =` flag, no chunk-awareness in any scorer, no branch
asking whether it is seeing everything. Chunking is a separate front end (sec 6).

Why it has to be this way: partial-vs-complete missingness is only definable against a column's NA
count over every sample, so any function that classifies missingness needs the cohort. Making that a
precondition means the existing machinery is correct by assumption rather than correct by
inspection.

### 1.1 There is no `cohort_stats` argument

The wrapper does **not** inject cohort statistics into `calc_clocks()`. It **pre-imputes the chunk**
before handing it over, which makes the scoring path a pure function of the matrix it is given. The
equivalence, per needed column:

| Cohort-wide state | Wrapper does | `calc_clocks()` sees | Same as single pass? |
|---|---|---|---|
| clean | nothing | clean, raw | yes |
| partial NA | fills with the cohort mean | clean column, `partial_cache = NULL` | yes -- `observed_panel()` only *orders* cached-before-raw, never branches on provenance |
| all NA | cannot fill | all-NA -> `all_na_cols` -> vendor ref / drop | yes -- the all-NA set is cohort-invariant |
| partial cohort-wide, all-NA *within this chunk* | fills | clean column | yes -- **this is the case the fill exists to erase** |

A column absent from `colnames()` in one chunk but present in another is not covered by any row of
that table, which is why colname agreement is an error and not a repair (sec 6.2).

Rejected alternative, so nobody re-proposes it: an argument on `calc_clocks()` carrying cohort means
plus the column classification. It works, but it puts a chunk-aware parameter on the public
one-cohort function, and it still would not have fixed coverage (sec 5). A hoist plus two thin front
ends gets the same sharing with no public surface.

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
| 5 | hoist the shared checks out of `calc_clocks()` | not started (new) |
| 6 | chunk sources + the two-pass wrapper | not started (new) |

Execution order: **0 -> 5 -> 2 -> 3 -> 6 -> 4.** Phase 0 is a standalone prerequisite (sec 4).

## 4. Phase 0 -- `partial_cache` becomes a mean vector

`build_partial_cache()` (`R/missingness.R`) returns `slideimp::mean_imp_col()`, i.e. a materialized
`n x k` imputed matrix held beside `DNAm`. It should be a **named numeric vector of column means**,
applied at read time.

- `build_partial_cache()` -> accumulable `(sum, n_non_na)` per column, reduced to means.
- `cached_cols()` (`R/utils.R`) keys off `names(means)` instead of `colnames(partial_cache)`.
- `observed_panel()` fills NAs in `DNAm[, cached]` from the vector instead of `cbind`-ing a
  pre-filled block.

Independent of chunking and worth doing regardless: it removes a second large object from the
ordinary in-memory path. Not a chunking blocker, because sec 1.1 means the wrapper hands over an
already-filled chunk and `partial_cache` is `NULL` there.

## 5. Phase 5 -- the hoist, and why coverage is the whole point

Pre-imputation is score-equivalent but it **erases the record of what was imputed.**
`panel_sample_miss()` (`R/coverage.R`) counts NAs in the *raw* matrix over cached columns; with
`partial_cache = NULL` every count is zero. A chunked run would therefore report
`score_imputed_partial = 0` for every clock, an all-zero `$coverage$sample_miss$score`, and
`check_row_coverage()` would never fire -- breaking "coverage is never reported for a sample it is
not true of" **silently**, which is the worst available failure mode.

Note the asymmetry: the **clocks** gate survives (`check_coverage()` reads `usable_cols`, and the
all-NA set is cohort-invariant); only the **samples** gate breaks.

So coverage must be computed from the pre-imputation NA pattern, which only pass 1 sees. That is not
a duplication of engine logic -- coverage is *already* computed once upstream of the scoring loop and
keyed by clock id -- it is a hoist. And it delivers the feature that motivated it: with coverage
built in pass 1, `min_samples_coverage` throws **before** any scoring.

### 5.1 The three tiers

Split `calc_clocks()` so both front ends compose the same internals.

**Data-independent** (once, either front end): `resolve_clocks()`, `resolve_clocks_sequence()`,
`drop_routed_members()`, `resolve_DNAm_extra()`, the covariate union and its missing-pheno abort,
`load_mc_assets()`, `clock_panels()`.

**Per-row predicates** (chunk-local is already correct): rownames non-NULL / non-NA, the all-NA
sample abort (`row_miss == ncol`, a per-row test).

**Cohort-set predicates** (pass 1; each of these is *wrong* when handed a chunk today):

| Fact | Today | Chunked |
|---|---|---|
| rowname uniqueness | within one matrix | across every chunk |
| pheno subset | `rownames(DNAm)` vs `pheno[[pheno_id]]` | union of all chunk ids, so the error names every missing id once instead of chunk 7's |
| column classification | `col_miss == nr`, chunk `nr` | cohort `n` |
| column means | over the matrix given | accumulated |
| both coverage gates | per matrix | cohort-wide, pre-scoring |
| colname agreement | n/a | across blocks (sec 6.3) |

The scoring loop becomes an internal taking a pre-imputed matrix and nothing else. The wrapper calls
that internal, **not** the public `calc_clocks()` -- otherwise it pays for a coverage computation it
must discard and gets per-chunk error messages.

### 5.2 The invariant this buys

> Holding the cohort facts fixed, scoring any row subset yields **bit-identical** rows to scoring the
> whole cohort.

Testable, and the only thing that actually proves chunk-safety: score a cohort, score it in three
chunks, `expect_identical`. Belongs in the always-on tier (sec 10).

## 6. Phase 6 -- chunk sources and the wrapper

### 6.1 An iterator contract, and two adapters over it

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

### 6.2 What we own, and what DelayedArray does

DelayedArray is used for **slicing only**: `dim()`, `dimnames()`, `[` + `as.matrix()`, and
`chunkdim()` for the geometry advisory. Four generics, about as stable as Bioconductor gets --
breakage exposure scales with surface area, and this is near the minimum.

`DelayedMatrixStats` is deliberately **not** a dependency. We own the accumulator, for three reasons
that are independent of dependency hygiene:

1. **Determinism.** `rowMeans2()` traverses in its own block order and lands ~5.6e-16 from a
   plain-matrix `colMeans()`. One accumulator called by both the single-cohort and chunked paths makes
   them agree *by construction*, which is what keeps sec 5.2's `expect_identical` viable. (Same
   failure mode as duckdb's threaded `avg`, arriving from a different direction.)
2. **One traversal, not many.** Measured: `rowMeans2()` + `rowSums2(is.na())` is two sweeps at 2.33s;
   a fused loop -- slice, `matrixStats` kernels on the plain block, accumulate -- is **1.31s**,
   agrees to 5.6e-16, and yields the value range and the id vector in the same pass. Pass 1 needs
   about six quantities; each DelayedMatrixStats call is another sweep of the disk.
3. **No global mutation.** `setAutoBlockSize()` is a session-wide option; a package must not set one
   for its user. Owning the slicing means block size comes from our own argument.

`matrixStats` is already in Imports, so the per-block kernels cost nothing new.

Block size is exposed as a **memory budget**, not a block count, since RSS is the actual constraint:
`n_samples_per_block = budget_bytes / (8 * n_probes_needed)`.

### 6.3 Orientation and geometry -- assume neither

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

### 6.4 The two passes

Pass 1 is **one fused traversal**, accumulating per block: per-column `(sum, n_non_na, n_na)` over the
needed panel; per-sample panel miss; the sample ids; the within-block row predicates; column
agreement; and the **value range**, which is a free beta-vs-M-value check -- `check_DNAm()` validates
`mode = "double"` but not the [0, 1] range, so a minfi `GenomicRatioSet` whose assay holds M-values
currently scores plausible garbage. Then once: global id uniqueness, pheno subset against the union
id set, cohort column classification, means, the full coverage record, and both coverage gates. A
gate failure throws **here**, before any scoring, and the surviving worklist lets pass 2 read fewer
columns.

Pass 2, per block: re-read, fill partial NAs from the pass-1 means, call the scoring internal, keep
only the score matrix. Assemble, reorder by `sample_id`, attach the pass-1 coverage, finalize any
cohort reductions (sec 7).

Naming: `calc_clocks_chunked()` follows the existing family. The wrapper stores a chunk label per
sample (sec 9).

### 6.5 Column narrowing conflicts with Zhang2019

The obvious RSS win is reading only the union panel, but a `sample_scale` clock takes its moments
over **every** probe in the input. Reuse the existing declared-recipe predicate
(`needs_full_panel()`, already in the parity tier) rather than a clock list. Whole-array column means
cost ~7 MB even at 866k probes, so filling full width when that predicate is true is affordable.

## 7. Phases 2 and 3 -- cohort reductions leave the scoring loop

Exactly **2 of 129** catalog clocks are cross-sample: `DNAmPhysAge` and `DNAmPhysAge_years`
(`cross_sample_at = 11`). Pre-imputation does not touch them.

`score_PhysAge()` reduces over an `n x k_surrogate` matrix of **derived surrogate scores**, not over
betas -- `scale(raws)`, then `rowSums`, then a second `scale(phys)` before the polynomial. Two nested
cohort z-scores, plus an `n >= 2` guard.

So the chunk-safe output is the surrogate raws (negligible memory) and the reduction is a cheap
post-assembly step. Phase 2 is the engine acting on `cross_sample_at`; Phase 3 emits per-sample
intermediates from the scoring loop and finalizes after assembly. For a complete cohort
`calc_clocks()` composes both and its output is unchanged.

Consequence to accept: under chunking, `DNAmPhysAge` is not a finished column until assembly, so the
wrapper's last step is a **finalize**, not a bare bind.

## 8. What `augment()` owns

Cross-sample derivations that are not a clock's definition -- age acceleration / residuals,
user-requested z-scores -- happen in `augment()`, after binding, never in the scoring loop. This is
what `dev/detail-plan.md` sec 7.1 means by the third `cbind` gate being "dissolved". `augment()` does
not exist yet; CLAUDE.md listing it among the verbs is stated intent (DECISIONS 2026-07-24).

## 9. Phase 4 -- `rbind`

`rbind.mc_result` refuses today (`R/generics.R`), and that is correct while any score column can be
cohort-dependent. After Phase 3 no score column is, and it can admit records under gates:

1. **Disjoint** `sample_id` sets (`cbind` requires equal; this requires disjoint). Collision throws
   and names the ids -- this also catches lazy per-chunk `sample1..N` labelling, which is what the
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

One residual asymmetry worth stating plainly: a record bound from per-time-point runs carries a
`DNAmPhysAge` that was z-scored **within** batch. That column is meaningful but is not the column a
single-cohort run would produce. The batch label is what lets a reader see that; do not try to
recompute it at bind time.

## 10. Tests

Always-on, no new dependency -- the in-memory block iterator (sec 6.1) is what makes this possible,
and it is the whole reason that source exists:

- **Chunk invariance** (sec 5.2): one cohort scored whole vs in three chunks, `expect_identical`.
  Include a clock with partial NAs, one with a fully absent probe, and `DNAmPhysAge`.
- **The chunk-dependency case that used to be wrong**: a column partial cohort-wide but all-NA within
  one chunk must land as `imputed_partial`, not `imputed_full`, in the chunked record's coverage.
- **Coverage is real under chunking**: `score_imputed_partial` and `sample_miss` from a chunked run
  equal the single-pass values. This is the regression test for sec 5's silent-zero failure.
- **Gates fire before scoring**: a cohort failing `min_samples_coverage` throws without reading
  pass 2 (assert via a source that records how many times each block was read).
- **`rbind`**: disjoint ids bind; colliding ids throw; batch labels survive; `[` after `rbind` keeps
  them.

- **Value range**: an M-value matrix (roughly [-10, 10]) is rejected by pass 1, not scored.

Skipped when Suggests are absent, per existing practice: the `DelayedArray` source, including one
transposed input and one deliberately hostile chunk geometry. Optional and cheap if wanted later --
a cross-check that our fused accumulator agrees with `DelayedMatrixStats::rowMeans2()` within a
bounded diff (they sit ~5.6e-16 apart, so this is a bounded gate, never `identical`). That test is
the only conceivable call site for `DelayedMatrixStats`, which is otherwise not a dependency.

Parity is untouched -- it scores complete cohorts and is not part of "run the tests".

## 11. Non-goals

- No random-access source contract. Two sequential scans is the whole requirement.
- No chunk-size autotuning. Take the block size from the caller (or the source's own geometry) and
  say what was used.
- No parallelism in the first cut. Chunking trades speed for RSS by design; adding cores changes the
  peak-memory story it exists to fix.
- No `rbind` across cohorts *reconciling* anything -- no re-imputation, no re-z-scoring, no merged
  coverage denominators. It binds and labels.
- No new identity key. Ids are the identity; `batch_set_id` was removed for good reasons
  (DECISIONS 2026-07-24) and a batch **label** is not its return.
