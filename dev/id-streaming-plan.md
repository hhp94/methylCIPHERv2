# Chunked scoring, prepared inputs, and record binding

The one tracked design doc. It covers **only work that is not built yet**, plus the built seam that
work composes against. Shipped behavior is specified by `CLAUDE.md`'s invariants and by the code;
this file never restates them.

**Current truth only.** Superseded design is not annotated here -- its history lives solely in
`dev/DECISIONS.md`. When code and this file disagree, the code is truth: fix this file and record
the reconciliation in DECISIONS.

Four DECISIONS entries cite **phase numbers**; see sec 3 before renumbering. Section numbers are
cited by nothing and may move freely.

**Picking this up cold?** Read sec 3 (phase table), sec 4 (the seam, built), then sec 9 (Phase 7,
next). **Sec 5 is parked -- read it as a measurement record, not as a plan** (see the pin at its
head).

## Status

**The seam is cut, the sample axis is classified, every clock is chunk-safe, and `rbind` assembles
separately-scored blocks. No chunk source is built, so nothing in sec 5 is reachable from a user
session** -- the supported path for a cohort that does not fit is project with `clock_cpgs()`,
block, score, bind (sec 8). Phase 4 added `refinalize_clocks()`, a derived batch label and a batch
axis on `$coverage$per_clock`; every other internal named here is unexported.

| | |
|---|---|
| Built | Phase 1 (by other means), **2** (sec 6.1), **3** (sec 6.2), **4** (sec 8), **5** (sec 4), Phase 7's `normalize=` half (sec 9.3); both `src/` kernels (sec 5.5) |
| Designed, not built | the `prep()` record (sec 9) -- **next** |
| **Parked** | **Phase 6** (sec 5) -- the streaming front end (DECISIONS 2026-07-30) |
| Proven | chunk invariance over **every** clock, always-on tier, `test-chunk-invariance.R` (sec 4.2, 10) |
| Measured, not built | the access engine and input canonicalization (sec 5.7) |

**No live correctness gap.** Phase 3 moved the two PhysAge reductions past assembly, so a chunked
front end may score the whole callable pool. `spec$cross_sample` names the clocks the scoring loop
*defers*, not clocks it gets wrong.

**Panel sizes below predate the 2026-07-30 sync**, which took the catalog to 134 clocks (120
callable, 91 bundled callable). Sec 5 and sec 9.1's measured tables, and sec 4.2's sequence count,
were taken before it and have not been re-measured -- read them as dated measurements, not as
current counts.

---

## 0. Chunking and binding are two problems

| | Chunking | Binding |
|---|---|---|
| Goal | cap resident set size | reuse the `mc_result` verbs over records scored separately |
| Speed | may be sacrificed | irrelevant |
| Fill values | **one** cohort mean set, shared by every chunk | **one set per record**, deliberately different |
| Result | numerically equal to a single-pass run | a labelled union of differently-imputed batches |

The chunking column describes Phase 6, which is parked -- **no built path produces it**, so do not
read "numerically equal to a single-pass run" as a claim about `rbind`. A user who blocks by hand
and binds gets the binding column, per-batch fill and all.

The two are still separate problems, and binding is the wider one. Its motivating non-chunk case:
three time points, not nested, imputed per time point on purpose, then bound. Chunk reassembly and
that union produce disjoint id sets and identical clock columns alike, so **the id gates cannot
tell them apart** -- only the recorded batch can, which is why the label exists and why a
re-finalize is opt-in rather than automatic (sec 8).

## 1. `calc_clocks()` takes a complete cohort

Cohort-completeness is a **precondition, not a parameter**: no `chunked =` flag, no chunk-awareness
on the public surface, no scorer asking whether it is seeing everything. Partial-vs-complete
missingness is only definable against a column's NA count over every sample, so any function that
classifies missingness needs the cohort. A cohort that does not fit is blocked by the caller and
bound (sec 8), not handed to `calc_clocks()` whole.

`calc_clocks()` **hard-errors on a non-matrix array-like** (a `DelayedMatrix`, an on-disk handle):
`check_DNAm()`'s `assert_matrix` already refuses one. What it does not do is say what to do
instead, and the answer moved with the park -- project with `clock_cpgs()`, block, score, bind. A
dedicated branch carrying that message is worth having and is **not** scheduled.

### 1.1 The cohort facts reach the scoring internal as arguments

Both front ends call `score_cohort()`, and every cohort fact reaches it inside `facts` -- chiefly
`partial_fill`, a named numeric vector whose **names are the cohort-partial columns** and whose
**values are the cohort means**.

`partial_fill` does two jobs, and the second is load-bearing:

1. it supplies the fill values;
2. **its names are the column classification**, a cohort fact that must not be re-derived from a
   block.

There is exactly one producer and it is never optional: `mc_cohort()` derives `partial_fill` and
`score_cohort()` consumes it, on both front ends. What differs is only where `mc_cohort()` got its
numbers, and `score_cohort()` cannot tell. It never classifies a column itself, so there is no
local computation for a caller to correct afterwards.

Correspondingly, `scan_missing_cpgs()` is called by `mc_cohort()` and by nothing else. Its row half
(the all-NA-sample abort) is per-row and correct against a block; its column half is a cohort fact.

Per needed column, as a block sees it:

| Cohort-wide state | In `partial_fill`? | Block sees | Same as single pass? |
|---|---|---|---|
| clean | no | raw, no NAs | yes |
| partial NA | yes | NAs filled from the cohort mean, into the cache | yes |
| all NA | no | absent from `usable_cols` -> vendor ref / drop | yes -- the all-NA set is cohort-invariant |
| partial cohort-wide, **all-NA within this block** | yes | every row filled from the cohort mean, every row counted as imputed | yes -- **this is the case the vector exists to get right** |

A column absent from `colnames()` in one block but present in another is not covered by any row of
that table. Column agreement is therefore an **error**, never a repair; an explicit intersect
opt-in is allowed, a silent one is not.

**The block is handed over raw and filled inside the scorer.** `count_sample_miss()` counts NAs in
the **raw** matrix over the cached columns, so a pre-filled block reports an all-zero `sample_miss`
and `check_row_coverage()` degrades into a per-clock constant -- silently. Pre-imputing a block
outside the scorer is a non-goal (sec 11).

## 2. Imputation policy

Partial NA on a present probe fills from **the user's own data only** -- the cohort mean, or the
chunk-set mean when the caller deliberately batched. A fully absent probe takes the clock's vendored
ref or is dropped by policy. A cohort-independent fill is not an option: a vendored constant is a
different estimand, so cohort-mean fill is the only exact way to honor the policy, which is what
makes the two-pass design mandatory rather than a preference.

## 3. Phase identifiers are stable

Phase numbers are **names, not an order.** Four tracked citations refer to them.

| Phase | Meaning | Status |
|---|---|---|
| 1 | sample-id resolution | **done** by other means -- mandatory rownames (DECISIONS 2026-07-24) |
| 2 | kind-1 / kind-2 split (chunk-safe vs cohort-reducing clocks) | **done** -- `split_cross_sample()` -> `spec$cross_sample` |
| 3 | every cross-sample op leaves the scoring loop | **done** -- `finalize_cross_sample()` |
| 4 | `rbind` gates | **done** -- `R/bind.R`, `refinalize_clocks()` (sec 8) |
| 5 | the seam: both front ends compose the same internals | **done** -- `R/score_cohort.R` |
| 6 | the streaming front end: input canonicalization + the two passes | **parked** (sec 5) |
| 7 | `prep()` -- the prepared-input record | **not started -- next**; its `normalize=` half shipped |

Phase 6 does not depend on Phase 4: the chunked front end assembles *fragments* and never calls
`rbind.mc_result`, so unparking it changes nothing about how `rbind` behaves. That independence is
what let Phase 6 park without blocking anything, and it still holds now that Phase 4 has shipped --
a fragment is not a record. Phase 7 is orthogonal to chunking but lands on the same seam.

## 4. Phase 5 -- the seam (built)

`calc_clocks()` is four pieces plus `construct_mc_result()`, divided by what each depends on. They
live in `R/score_cohort.R` and none is exported:

```
mc_spec(clocks, pheno_id, normalize, ext_data, ask)  # data-independent
mc_cohort(DNAm, spec, pheno, min_clocks_coverage)    # one scan -> the cohort facts; gates
score_cohort(DNAm, spec, facts)                      # a matrix and the facts -> scores + intermediates
finalize_cross_sample(scores, pending)               # the one cohort reduction, after assembly
```

The chunked front end is the same four, with `mc_cohort()` accumulating over blocks and
`score_cohort()` called per block. `score_cohort()` never learns which front end called it and
`facts` is the only channel between them, which is what makes sec 4.2 a property rather than a
coincidence. **The wrapper calls the internals, not the public `calc_clocks()`** -- otherwise it
pays per-block resolution and asset loading and emits per-block error messages.

`facts` carries `sample_id`, the aligned `pheno`, `usable_cols`, `cpg_list`, `partial_fill` and
`sample_moments`. `spec` carries the resolved clock ids, the compute `sequence`, `output_ids`, the
`normalize` decisions, the covariate union, `pheno_id`, the loaded `packs`, the `panels`,
`cross_sample` and `needs_moments`.

`score_cohort()` takes no pheno argument: it narrows `facts$pheno` to the rows in hand by matching
`rownames(DNAm)` against `facts$sample_id`, so a block can never be handed a misaligned pheno.

**Coverage does not move.** It is computed inside `score_cohort()` from the raw block plus the cache
by `compute_coverage()`, keyed by clock id. Nothing about the coverage machinery is chunk-aware and
nothing is patched after the fact.

**Two live ordering consequences.** The samples gate is warn-only and reads assembled counts, so
`check_row_coverage()` runs *after* the scoring loop -- as does `check_score_values()`, which the
chunked front end must call **once over the bound scores**, never per block. And asset loading is in
the data-independent tier, so `load_mc_assets()` -- which may prompt and download -- runs before
`check_DNAm()` and the pheno checks.

`score_cohort()` factors as `score(prep(...))`; Phase 7 makes that boundary public.

### 4.1 The three tiers

**Data-independent** (`mc_spec()`, once, either front end): `resolve_clocks()`,
`resolve_clocks_sequence()`, `resolve_normalize()`, `drop_routed_members()`, the covariate union,
`load_mc_assets()`, `clock_panels()`, `split_cross_sample()`, `note_full_panel_clocks()`.

**Per-row predicates** (block-local is already correct): rownames non-NULL / non-NA, the all-NA
sample abort. The missing-pheno abort sits here too -- it is a fact about the `pheno` argument, not
about any matrix, and is raised before a block is read.

**Cohort-set predicates** (`mc_cohort()`; each is *wrong* when computed from a block):

| Fact | Single pass | Chunked |
|---|---|---|
| rowname uniqueness | within one matrix | across every block |
| pheno subset | `rownames(DNAm)` vs `pheno[[pheno_id]]` | union of all block ids, so the error names every missing id once |
| column classification | `n_obs == 0` / `0 < n_obs < nr` | the same test on `n_obs` over the cohort |
| column means | `sum / n_obs` | same, per probe (sec 5.3) |
| clocks gate | per matrix | cohort-wide, before pass 2 |
| colname agreement | n/a | across blocks |

The samples gate is warn-only, so it runs once after assembly on the concatenated counts.

### 4.2 The invariant this buys

> Holding the cohort facts fixed, scoring any row subset yields rows **equal** to scoring the whole
> cohort.

Lives in the always-on tier as `tests/testthat/test-chunk-invariance.R`, and needs no chunk source
-- it splits a matrix in memory and calls the seam directly.

**Measured: every clock agrees across a 3-block split, on a cohort carrying all three missingness
shapes at once, with no exclusion.** The assertion is `expect_equal`, not bit-exactness: over all
101 sequence entries, 94 differ in the last bits, worst **9.3e-10 absolute** (`DNAmB2M`, scale
6.7e6) and worst **5.1e-13 relative** (`DNAmPhysAge`) -- ~200x inside `PARITY_REL_TOL`.

The drift is not the fill and not the mean: `observed_panel()` returns `identical()` `cols` and
`values` for a block and for the whole cohort. It is `%*%` in `linear_score()` -- reference BLAS
`dgemm` blocks by `M`, so the row count changes accumulation order. Bit-exact under
`options(matprod = "internal")`. **Consequence for Phase 6:** an explicit bound is the only option,
not a fallback, and under real chunking the *means* also move in the last bits because block
boundaries differ from a single pass.

## 5. Phase 6 -- the streaming front end (parked)

> **Pinned 2026-07-30. This section is a measurement record, not a plan** -- read it for the numbers
> and the kernel contracts (sec 5.5, which describes built code), not as work to pick up. Nothing
> below is deleted and nothing below is scheduled.
>
> **What the front end is worth, stated plainly:** exactly one case survives -- a cohort too big to
> hold **even projected to the requested panel**, *and* carrying NA. The NA is the load-bearing
> half. Cohort-mean fill is the only thing sample-blocking can get wrong, so on a complete matrix a
> user-side loop is already exact and this front end buys nothing but convenience.
>
> **Why that is niche.** The projection lever below is stronger than the section assumed: the 87
> bundled callable clocks union to 20,430 scoring CpGs (38,540 with the BMIQ norm panel), which is
> 1.5 GB at n=1e4 against 64.5 GB for the whole 866k array. A request for every bundled clock is
> resident past any cohort that exists. And sec 5.2's store rested on gzip being unseekable --
> but a cohort genuinely too big to hold is not sitting in a `.csv.gz`, it is in HDF5/Zarr/TileDB,
> where two passes are free and the store has no justification left.
>
> **What covers the rest: Phase 4** (sec 8). `clock_cpgs()` already ships, so a user with a
> random-access store can project and block today; assembly is the only missing piece.
>
> Full reasoning, and the two findings that outlive this park (Bioconductor's HDF5 layout, the
> float32 non-goal), are in DECISIONS 2026-07-30.

The target: a small box -- **2-4 cores, 4-8 GB RAM** -- against a tall `.csv.gz` far larger than RAM,
with no ETL the user has to run. On that box the run is **deflate-bound and serial** (gzip inflate is
whole-stream and single-threaded), not compute-bound, which is why the design minimizes full passes
over the source and does not parallelize scoring (sec 11). If the data fits in memory, none of this
section applies.

**The resident set is the projected panel, never the file.** Only `needed_cpgs x n_samples` has to be
in memory, so the sizing rule is

```
chunking is required iff  3 * 8 * |needed_cpgs| * n_samples  >  budget_bytes
```

(the `3` is sec 5.3's copy multiplier). Against a 1e6-probe array at n=1e4 -- 80 GB as doubles -- with
a 4 GB budget:

| request | \|panel\| | projected, n=1e4 | pass-2 block | blocks |
|---|---|---|---|---|
| all 101 bundled clocks | 39,025 | 3.1 GB | 4,271 | 3 |
| PCClocks (14) | 78,464 | 6.3 GB | 2,124 | 5 |
| SystemsAge (13 organs, one shared panel) | 125,175 | 10.0 GB | 1,331 | 8 |
| PCBrainAge | 357,852 | 28.6 GB | 466 | 22 |

So **requesting a PC-scale panel is what forces chunking, not cohort size**: the bundled set is
resident to n ~ 4,300 while PCBrainAge chunks from n ~ 470. `resolve_cpgs()` already dedups panels,
which is why 13 SystemsAge organ clocks read one panel and not 13.

### 5.1 Input canonicalization

The stored shape and the package shape are opposite, and both are fixed:

- **Stored / accepted: probes x samples ("tall").** This is what the ecosystem produces -- the
  Horvath server csv, GEO series matrices, Bioconductor assays, HDF5 written by
  `writeHDF5Array()`. Accept it as the canonical input orientation.
- **Internal: samples x probes**, `rownames` as identity. Unchanged.
- **Transpose on the small projected result, never the array.** A pass-2 block is a panel-projected
  sample block, so `t()` lands on something the size of the block.

**dtype is float64, closed.** float32's ~1e-7 relative error is three orders coarser than
`PARITY_REL_TOL = 1e-10`.

**Null spellings must be declared, and the assert is ours.** Both readers default to making a
column `character` the moment they see a spelling they do not know, which surfaces as
`check_DNAm()`'s `mode = "double"` assertion talking about a matrix with nothing pointing at the
file, the column, or the token. Measured defaults:

| spelling | `NA` | `NaN` | `nan` | empty | `null` | `None` | `#N/A` | `inf` / `-inf` / `Infinity` |
|---|---|---|---|---|---|---|---|---|
| `data.table::fread` reads as NA | yes | no | no | no | no | no | no | no |
| duckdb `read_csv_auto` reads as NULL | **no** | no | no | yes | no | no | no | no |

So: pass an explicit null-string set
(`NA`, `NaN`, `nan`, empty, `null`, `None`, `NULL`, `#N/A`), then **assert every data column came
back numeric and name the offending column plus token when it did not**. With the set declared,
both readers resolve correctly and `inf` / `-inf` / `Inf` / `Infinity` all survive as `+/-Inf`,
which the `any_inf` diagnostic depends on -- **never map `Inf` to `NA` at read time.**

**Identity: the id column is column 1 of a tall file when header field 1 is empty or does not look
like a CpG id.** Verified against `write.csv` (blank quoted), pandas `to_csv()` (blank), pandas
`index_label=`, polars, and GEO `ID_REF`; a blank first field reads as `V1`. Sample ids are
`header[-1]`.

**Orientation detection is the coverage check, not a heuristic.** Resolve the request's panel first,
read only the header, and test it as CpG ids; on failure read only column 1 and test that. If
neither clears `min_clocks_coverage`, throw the ordinary coverage error -- it already says what is
wrong. Take an explicit override argument. No data is read to decide.

### 5.2 The store is mandatory, and it is a duckdb file

`duckdb` and `DBI` are already in Suggests. duckdb owns **file reading, projection and chunked
access**; the package owns every float that reaches a score. `data.table` does **not** become a
dependency.

The user hands us `.csv.gz`; we ingest it **once** into a duckdb file and every pass reads that. The
store is not an optimization and not request-scoped scratch -- it is what makes more than one pass
affordable.

**Why: gzip is whole-stream, a columnar store is per-page.** A `.csv.gz` cannot be seeked, so *every*
pass over it re-inflates and re-parses the whole file no matter how little of it that pass needs.
Cohort-mean fill makes two passes non-negotiable (sec 2), so the choice is between paying the deflate
twice per request and ingesting once. A per-request spill of just the resolved panel was considered
and rejected: it is smaller (3.1-28.6 GB against ~80 GB) and free as a byproduct of pass 1, but its
panel depends on the clock request, so it re-deflates the source once per request -- and the realistic
session is a sequence of requests over one file (DECISIONS 2026-07-30).

**duckdb, not parquet:** it handles the wide schema (1e4+ sample columns) better than parquet's
per-column-per-row-group footer metadata, and a `PRIMARY KEY` on the CpG id gives an index plus an
ingest-time duplicate-probe check. Not measured against parquet at width -- a decision, not a
benchmark.

Ingest contract:

- one `CREATE TABLE AS SELECT * FROM read_csv(...)` straight off the `.gz`, no extraction step
  (measured 7.26s for the 190 MB reference cohort);
- the CpG id column is the `PRIMARY KEY`, so **duplicate probe ids surface at ingest, once**;
- a `rid` (`row_number() OVER ()`) is added at load and every ordered read uses `ORDER BY rid`. A
  plain `LIMIT/OFFSET` scan was measured to return file order, but order without `ORDER BY` is not
  guaranteed and `rid` costs nothing;
- sample ids are the header (sec 5.1), so id uniqueness and the pheno subset are settled before any
  data row is parsed (sec 5.4).

**The store is user data, not a reclaimable cache** -- ~80 GB for a 1e6 x 1e4 cohort, since duckdb
costs about what uncompressed doubles cost and full-precision betas do not compress. It must not go
to `R_user_dir(..., "cache")` beside the packs; lifecycle and consent are open (sec 5.8).

### 5.3 The two passes have two different axes

This is the core of Phase 6. SQL is good at one axis of a tall table and bad at the other, so each
axis is served by the tool that fits it.

**Pass 1 chunks by probes.** One `rid BETWEEN` query yields all samples x a probe range; **transpose
it to samples x probes** (sec 5.1 -- the transpose lands on the chunk, never the array) and run the
existing `col_stats()` kernel on it. That is the orientation the kernel was built for, so none of its
contract changes:

- `stats` is per-probe `sum` / `n_obs`. Because probe blocks are **disjoint in probes**, each probe's
  mean and observed count are **complete within its block** -- there is **no cross-block accumulator
  for the CpG axis at all**, which is the whole reason pass 1 cuts this way.
- `row_mean` / `row_m2` / `row_obs` are per-sample, and `cols` plus the complement sweep behave exactly
  as sec 5.5 describes. These are the only values that accumulate across chunks.
- The value-gate flags (`any_lt0`, `any_gt1`, `any_inf`) and `overflow_col` keep today's meaning, with
  `overflow_col` still naming a CpG.

**Reading the chunk untransposed does not work.** On a probes x samples block `stats` becomes
per-sample `sum` / `n_obs` -- no second moment, so the `sample_scale` sd is unavailable -- and `cols`
would select samples rather than the panel.

**Do not write row-wise aggregates in SQL** -- and the reason is arithmetic, not width. A
`sum(S1)+sum(S2)+...` chain is one expression tree a node deep per sample and hits duckdb's
`Max expression depth limit of 1000` at 2,000 columns. *Returning* 1e4 columns is fine (measured to
20,000); writing an expression over them is not. So the limit never applies here: pass 1 issues a bare
projection, and `col_stats()` does every sum in C++. This is not a workaround for the depth limit, it
is sec 5.2's rule -- no SQL aggregate ever produces a number that reaches a score.

**The probe filter is on iff `spec$needs_moments` is FALSE.** The two Zhang2019 arms z-score each
sample over **every** probe in the input matrix, so when one is requested pass 1 reads unfiltered and
`cols` is the needed subset within each chunk -- exactly `col_stats()`'s split (sec 5.5), `stats`
panel-scoped while the row accumulators cover everything. Projecting first would silently redefine that
z-score: substituting the panel union was measured at 1.8e1 absolute / 82% relative.

**Per-sample moments come from the kernel, not from SQL.** `avg()` plus `stddev_samp()` over the sample
columns is the cheap axis (1.99s over 20,000 columns) and would work, but duckdb does not treat
`+/-Inf` as missing while `col_stats()` does (sec 5.5), so the two disagree on a file carrying an `Inf`
-- and the abort that would make this moot is itself computed in pass 1, so the SQL route needs the
value gate ordered ahead of it. Pass 1 already scans every probe whenever `needs_moments` is TRUE, so
the kernel gives the moments in the same traversal under the package's own semantics. The divergence is
recorded because it is a real property of the store, not because the design relies on it.

**The per-sample sd needs a cross-chunk merge that does not exist yet.** `scan_missing_cpgs()`
(`R/missingness.R`) reads `row_mean` / `row_m2` / `row_obs` off a single sweep. Under chunking they
arrive per chunk and must be combined pairwise (Chan's parallel Welford) -- a few lines, but it is the
**only** accumulator pass 1 needs, since per-probe stats are complete within a probe block.

Derived once at the end of pass 1 -- the same values `mc_cohort()` derives today, so Phase 6
replaces the accumulator and nothing below it:

| derived | from |
|---|---|
| `all_na_cols` | `n_obs == 0` |
| `partial_fill` | per-probe mean over columns with `0 < n_obs < n_samples`; names are the classification |
| `usable_cols` | present-needed minus `all_na_cols` |
| `cpg_list` | `resolve_cpgs(usable_cols, panels)` |
| clocks gate | `check_coverage(cpg_list, min_clocks_coverage)` -- **throws here**, before any scoring |

Id uniqueness and the pheno subset are **not** on this list -- they are settled from the header before
pass 1 begins (sec 5.4).

**Pass 2 chunks by samples.** Per block: one panel-projected query (semi-join the panel id table x
the block's sample columns, `ORDER BY rid`), transpose to samples x probes, build the block cache
with `fill_imp_col()`, call `score_cohort()`. Coverage comes from the unchanged
`compute_coverage()` on the raw block.

**Block size is a memory budget, not a block count, and the divisor is three copies:**

```
n_samples_per_block = budget_bytes / (3 * 8 * n_probes_needed)
```

Peak per block is the block itself, plus `build_partial_cache()`'s slice over the partial columns, plus
`pack_design()`'s `n x |panel|` copy (~2.9 GB for PCBrainAge at n=1000, sec 5.6). **The cache is not a
sliver**: at any realistic NA rate nearly every column carries at least one NA, so it is a full second
copy of the block. Report what was used; do not autotune.

**Assembly is a concatenate, never a replace and never a sum:**

- score matrices and `pending` intermediates -- concatenate by clock id (disjoint rows)
- `sample_miss$score` / `$norm` -- concatenate the per-block vectors (disjoint rows)
- **every** `per_clock` field is a CpG count off the shared cohort facts, so every block already
  computed the identical value -- take it

Then reorder by `sample_id`, call `finalize_cross_sample()` on the assembled intermediates, run the
samples gate once, and run `check_score_values()` once.

This is **not** `rbind.mc_result` and must not become it (sec 8).

Naming: `calc_clocks_chunked()`.

### 5.4 Every cohort gate splits into accumulate-then-report-once

This is the bulk of the work and it is shared by every input shape. Today's gates report inline from
one shot; under chunking each would fire per block:

- `check_DNAm()`'s transposed heuristic -- once per block
- `check_col_values()`'s three warnings -- once per block, per flag
- the dead-row abort -- names block 7's samples, not the cohort's
- rowname uniqueness and the pheno subset -- within-matrix only

So `mc_cohort()` becomes `cohort_acc_new()` / `cohort_acc_add(block)` / `cohort_acc_finalize(pheno)`,
and **the single-pass path calls all three on one block** so there is literally one accumulator and
sec 4.1's table is structural rather than aspirational.

Fail-fast vs accumulate: an overflowed column sum **aborts immediately**, naming block and position.
Dead rows, duplicate ids and missing pheno ids **accumulate and abort once** after pass 1. The three
range/Inf flags OR-accumulate and warn once.

Structural pheno checks stay at `cohort_acc_new()` (before any block is read); the id-dependent ones
-- the covariate-NA warning, the subset abort -- move to finalize.

**Two of these are not accumulations at all under the tall orientation.** Sample ids are the header
(sec 5.1), so **rowname uniqueness and the pheno subset are settled by a header-only read** before any
data row is parsed, and duplicate *probe* ids surface as the `PRIMARY KEY` violation at ingest
(sec 5.2). What genuinely accumulates is the dead-row count and the three range/Inf flags. Sec 4.1's
"union of all block ids" is therefore just the header on this input shape.

**`pheno` is not the id source.** It is `NULL` for any request with no covariates (`mc_cohort()`), so
the sample set comes from the header and pheno only ever narrows it.

### 5.5 Kernels we own

Both live in `src/` and both are already built:

- **`fill_imp_col(obj, mean_vec)`** -- applies externally-computed per-column fills. **Mutates
  `obj` in place and returns nothing**: Rcpp's `NumericMatrix` wraps the caller's SEXP without
  copying. Its one caller, `build_partial_cache()` (`R/missingness.R`), always hands it a fresh
  `DNAm[, cols, drop = FALSE]` slice, and matrix subsetting always allocates. **The safety property
  is a caller contract, not a guarantee of the kernel, and it is stated only at the call site** --
  a second caller must re-establish the fresh-slice property.
- **`col_stats(obj, cols, row_moments)`** -- one traversal returning `stats` (2 x ncol, rows `sum`
  and `n_obs`), `row_obs`, `row_obs_complement`, `row_mean`, `row_m2`, the flags `any_lt0` /
  `any_gt1` / `any_inf`, and `overflow_col`.
  - **Every non-finite value is missing** -- NA, NaN and `+/-Inf` alike are skipped, not summed.
    `any_inf` says an `Inf` was among them: a data bug worth naming, not a different fate.
  - **It bails only on an overflowed column sum**, returning `stats` / `row_obs` / `row_mean` /
    `row_m2` as `NULL` with `overflow_col` set. So `overflow_col` is checked **first** and nothing
    else is read in that branch (`check_overflow_col()`).
  - The kernel **reports and does not decide**: `check_col_values()` raises the abort and warnings.
  - `row_moments = TRUE` adds Welford per-row mean/`m2`. **The two halves have different column
    domains and the kernel owns the split:** pass 1 sweeps `cols` fused with the moments, pass 2
    sweeps only the complement into the row accumulators. So `stats` and the value gates stay
    `cols`-scoped while the row accumulators cover all of `obj`, and each column is read exactly
    once. **`cols` must be unique** for this.
  - `cols` is 1-based indices into `obj`; `NULL` scans every column and skips validation. Positions
    in `stats` and `overflow_col` are relative to `cols`.
  - `row_observed()` `stop()`s on a `NULL` `row_obs` rather than trusting that the value gate ruled
    it out.

### 5.6 Sources

| Source | Dependency | Role |
|---|---|---|
| in-memory matrix + budget | none | the **always-on test vehicle** and the supported manual path. Real win even for a resident matrix: `pack_design()` copies an `n x \|panel\|` block, ~2.9 GB for PCBrainAge at n=1000 |
| csv / csv.gz (tall) | duckdb (Suggests) | **the low-friction one.** Ingested to the store (sec 5.2) straight off the `.gz` -- no extraction step and no ETL the *user* runs |
| a list of files | duckdb (Suggests) | ingested the same way; cross-block colname agreement and duplicate ids are ingest-time checks |
| HDF5 / DelayedArray | Suggests (Bioc) | already random-access. Whether it stays a separate path or is simply ingested is open (sec 5.8) |

Dispatch on the class of one argument; no user-facing constructor.

`DelayedMatrixStats` is deliberately **not** a dependency, and `setAutoBlockSize()` is a
session-wide option a package must not set for its user -- owning the slicing is what keeps block
size an argument of ours.

### 5.7 Measured baseline

All on 16 threads. Reference cohort: 50,000 probes x 500 samples tall, 5% NA, methylation-shaped,
full precision -- 25M cells, raw doubles **190.7 MB**.

**Storage.** duckdb is consistently **0.99-1.14x raw doubles**; csv.gz varies entirely with data
compressibility (0.92x raw here, 0.27x on synthetic correlated data). HDF5 is 0.82-0.85x raw, so it
carries **no material storage win over csv.gz** -- its value is random access alone.

| | csv.gz | plain csv | duckdb store | HDF5 |
|---|---|---|---|---|
| size | 175.9 MB | 397.5 MB | 189.5 MB | ~160 MB |

**Access.** `fread` is far faster for one full read and unusable for chunked access, because
`skip=` is O(offset): 0.02 / 0.12 / 0.27s at offsets 0 / 25k / 49.9k.

| operation | fread | duckdb |
|---|---|---|
| one full read of the plain csv | **0.19s** | 2.86s (view) |
| ingest to native store | n/a | 3.25s (plain) / 7.26s (from `.gz`) |
| sweep 250 chunks of 10k rows x 10 cols | 69.0s | **0.58s** |
| sweep 500 chunks of 100 rows x all cols | 68.4s | **4.33s** |

**The real pass-2 query** -- 10,000-CpG panel semi-join x a 10-sample block -- is below timer
resolution; all 50 blocks projected **and** transposed sweep in **0.32s**. Pass-1 per-probe stats in
SQL: 0.81s all probes, 0.37s panel-only. Per-sample moments for all 500 samples: one query, 0.09s.

**Float agreement**: duckdb `sum()` vs R `sum(na.rm = TRUE)` is **1.04e-15 relative**, five orders
inside `PARITY_REL_TOL`; `count()` matches R's non-NA count exactly.

**Column-count scaling** (cells held at ~5M) -- the ceiling on the wide-tall schema:

| sample columns | 500 | 2,000 | 10,000 | 20,000 |
|---|---|---|---|---|
| ingest | 1.95s | 3.97s | 42.7s | 96.5s |
| per-sample moments | 0.07s | 0.14s | 0.80s | 1.99s |
| block read | 0.00s | 0.02s | 0.00s | 0.02s |
| row-wise `+` chain | ok | **fails** | **fails** | **fails** |

Ingest carries roughly **10 ms per column** of fixed overhead, so a 10k-sample table has a ~45s
ingest floor before any data volume.

**csv parse cost scales with columns, not cells** -- header-only parse is 0.01 / 0.05 / 0.20 / 0.50
/ 0.90s at 10k / 50k / 200k / 500k / 866k columns, and for a wide file header-only ~= full read.
This is the measurement behind sec 5.1's tall canonical orientation.

### 5.8 Open

- **The store at 1e4 sample columns -- the blocking benchmark.** Ingest carries ~10 ms per column of
  fixed overhead and 20,000 columns took 96.5s at only **5M** cells, so a 1e4-column ingest sits in an
  overhead-dominated regime nothing has been measured in. **Every n=1e4 figure in this section is
  arithmetic off a 500-sample baseline.** Measure ingest wall-clock and peak RSS, a probe-range scan,
  and a sample-column projection before building. If ingest is bad the fix is partitioning the store by
  sample block (e.g. 8 tables of 1,250 columns), which changes the schema -- hence "before". The
  fallback shape is long/tidy (`probe, sample, value`): constant width, both axes an ordinary
  `GROUP BY`, but ~8.5B rows for a 10k-sample EPIC cohort.
- **Store lifecycle, now the headline open**, because the store is mandatory (sec 5.2) and large. It is
  user data, so it must not land in `R_user_dir(..., "cache")` beside the packs, and it should probably
  have **no default location at all** -- the caller names a path or nothing is written. Naming
  (`<verb>_mc_<noun>`), `ask`-gated consent and a `clear_*` verb mirroring the assets surface are
  unresolved.
- **Whether the CpG `PRIMARY KEY` index beats a zonemap-assisted hash semi-join** for a 39k-value panel
  lookup. duckdb may plan a scan either way, and the index costs ingest time.
- **Probe-axis accumulation, which would remove sample-major blocking for exactly the clocks that need
  it.** A pack-scored clock is a batched weighted sum, and a weighted sum is additive over the
  contracted axis: accumulating `n x k` partial scores over probe chunks is ~10 MB resident at n=1e4
  instead of a 28.6 GB projected panel, and per-sample coverage counts are additive too. What needs a
  sample's whole panel at once is the non-additive set -- quantile/BMIQ normalization and the per-sample
  custom branches -- whose panels are small (bmiq gold is 21,368). The blocker is the **branch
  contract**, not the arithmetic: "a branch returns only its score" would become "a branch returns a
  partial score plus a merge rule", and every branch in the closed set would need an additivity
  classification. It sits beside sample-blocking rather than replacing it, so it is recorded, not
  scheduled.
- **Whether HDF5 / DelayedArray stays a separate path** or is ingested into the same store.
- **Provenance.** Whether a chunked record records its own geometry. Recording nothing keeps a
  chunked record fully `expect_equal` to a single-pass one, which makes the headline test
  unqualified; Phase 4's per-sample batch label may subsume it anyway.

## 6. Phases 2 and 3 -- cohort reductions leave the scoring loop (built)

### 6.1 Phase 2 -- the split, derived

Exactly **2 of 134** catalog clocks are cross-sample: `DNAmPhysAge` and `DNAmPhysAge_years`
(`cross_sample_at = 11`). The fill does not touch them.

`clock_cross_sample_at()` / `clock_is_cross_sample()` read the declared field,
`split_cross_sample()` partitions a compute sequence, and `mc_spec()` stores the cohort-reducing
half as `spec$cross_sample`. Nothing carries a clock list: the census test walks
`resolve_clocks("all")` and asserts the partition is total.

A sex-routed alias derives its axis from its members via `alias_cross_sample_at(members)` (the min,
`NA` when neither reduces) -- an alias scores each row with exactly one member, so a cross-sample
member makes the alias cohort-dependent. All 14 routed members are per-sample today, so the derived
value equals the old constant for all 7 aliases. The census test checks each alias against its
members rather than against a constant.

`cohort_zscore` is the only cross-sample op in the upstream corpus. `sample_scale` (Zhang2019) is a
within-sample z-score; `center_scale` (SystemsAge) reads vendored tensors, so it is per-sample too.

### 6.2 Phase 3 -- the reduction moves after assembly

`score_PhysAge()` reduces over an `n x k_surrogate` matrix of **derived surrogate scores**, not over
betas. The two clocks do not have the same number of reductions:

| clock | recipe | reductions |
|---|---|---|
| `DNAmPhysAge` | `cohort_zscore` (11) -> `row_sum` (12) -> `transform` (13) | one: `scale(raws)` |
| `DNAmPhysAge_years` | `cohort_zscore` (11) -> `row_sum` (12) -> `cohort_zscore` (13) -> `poly` (14) | two: `scale(raws)`, then `scale(phys)` |

`cross_sample_at` records only the **first** reduction -- it says where the chunk-safe prefix ends.
It is not a step pointer for the finalize, which follows the branch.

`physage_raws()` runs per-sample in the loop; `finalize_PhysAge()` performs the reduction.
`score_cohort()` returns a third element, `pending`; the loop routes a clock's output there instead
of into `scores` **iff it is in `spec$cross_sample`**, so the defer decision is the catalog's
declaration. The `n >= 2` guard lives in the finalize -- a 1-row *block* is fine, a 1-row *cohort*
is not.

Both front ends call `finalize_cross_sample()` **unconditionally**; `pending` is empty for the other
127 clocks. Nothing anywhere tests "is this PhysAge" or "are we chunked".

## 7. What `augment()` owns (not started)

Cross-sample derivations that are not a clock's definition -- age acceleration / residuals,
user-requested z-scores -- happen in `augment()`, after binding, never in the scoring loop.
`augment()` does not exist and is an **unbuilt idea, not a contract**; adding it is a new API
decision (CLAUDE.md).

## 8. Phase 4 -- `rbind` (built)

**This is now the front end for a cohort that does not fit.** With sec 5 parked, a user with a
random-access store projects with `clock_cpgs()`, blocks, scores, and binds -- so `rbind` is not a
convenience over chunking, it is the supported path. It lives in `R/bind.R`
(DECISIONS 2026-07-30, "Phase 4 gates").

### 8.1 One line decides every gate

> **Record what batching forces. Refuse what the caller chose differently.**

A per-batch fill regime is *forced*: cohort means are per-run by construction, so a batched user
cannot avoid it and refusing would refuse the whole feature. Everything else that can differ between
two records is a free argument with one obviously right answer across batches, so a difference is a
mistake and the bind says so.

**Per-batch imputation is what gets recorded, and only one of its two effects matters.** As noise it
is second-order: a block mean is unbiased for the cohort mean with SE `sd/sqrt(n_block)`, so on a
typical panel the linear predictor moves by ~1e-4. The effect worth naming is a probe **all-NA
within one block but partial cohort-wide** -- that block takes the vendored ref while the others take
the cohort mean, so every sample in it shifts by the same amount. That is a systematic per-block
offset, it is exactly the case `test-chunk-invariance.R` already isolates, and no bind-time
operation can undo it. The batch label plus the per-batch coverage (sec 8.3) is what makes it honest
rather than hidden.

### 8.2 The gates

`rbind.mc_result` admits records under exactly four:

1. **Disjoint ids.** `anyDuplicated()` over the concatenated `$provenance$sample_id` -- **never** the
   pheno id column, which is `NULL` for any run with no covariates. Collision throws, names the ids,
   and says to prefix them per batch; that also catches lazy per-block `sample1..N` labelling.
2. **Identical score columns**, reordered to the first record or thrown. `output_ids` is the compute
   sequence minus routed members and the sequence is a deterministic function of the requested ids,
   so this **already implies identical panels** -- there is no separate denominator gate.
3. **Identical `pheno_id`.**
4. **Identical `$provenance$normalized`.** Same clocks with different `normalize=` give identical
   columns and identical scoring panels over genuinely different numbers. The experimenter comparing
   normalized against raw is already refused by gate 1 (same samples, colliding ids) and wants two
   columns side by side anyway; the only caller who reaches gate 4 is one whose loop varied the
   argument across chunks.

There is no `force =`, and record *contents* are never compared to detect that two arguments are the
same record -- the id sets are the whole test. Pack-version drift is undetectable by design
(`payload_hash` never reaches a record) and is not defended against.

### 8.3 The batch label

`$provenance$batch` is a character vector aligned to `sample_id`, present on **every** record. It
is **derived, not assigned**: `construct_mc_result()` sets it to `batch_hash(pheno[[pheno_id]])`,
12 hex of `xxhash64` over the pheno's id column. There is no `batch =` argument on `calc_clocks()`
and no way to choose a label (DECISIONS 2026-07-30, "The batch label is derived").

`$pheno` is always materialized -- the id column alone when the caller supplied none -- so the hash
always has a column to read. The id column and **not** the whole frame: hashing covariate values
would make correcting one subject's age rename the batch.

Because the label is a function of the record's own ids, `rbind` has no label policy. It mints
nothing, renames nothing, renumbers nothing, and **drops** argument names -- `split()` names its
result by factor level, so refusing them would kill
`do.call(rbind, lapply(split(seq_len(n), g), score))`, which is the blocking idiom sec 8 exists to
support. Re-association is therefore exact, not merely order-preserving:
`rbind(rbind(r1, r2), rbind(r3, r4))` and `rbind(r1, r2, r3, r4)` return `identical()` records.
The hashed value is canonicalized first (sorted with `method = "radix"`, unnamed, `serialize =
FALSE`) so the label follows the id *set* rather than its order, the locale or R's serialization
version. Gate 1 makes the id sets disjoint, so two batches sharing a label needs a 48-bit collision
-- 1.8e-11 at 100 batches, 1.8e-9 at 1000. There is no fifth gate on labels, and adding one back
would guard something gate 1 already forbids.

`sim_DNAm(suffix =)` is **not** a batch argument, which is why it no longer shares the word. It
**suffixes the sample ids** (`sample1_T1`) so two simulated blocks are disjoint by construction and
clear gate 1 without the caller relabelling. It defaults to `NULL`: a default would rename every
sample. It is reported back on `$suffix` and is never fed into `calc_clocks()`.

**Coverage nests by batch, on every record.** `$coverage$per_clock` is batch -> clock -> record, and
`clocks_coverage()` is one row per (clock, batch) with `batch` **last** -- it is the key the frame
is on, but a hash reads as noise in front of `clock_id`. Merging is not available:
`score_imputed_partial` counts the panel CpGs in *that run's* partial cache, and two independently
scored batches almost never share an NA pattern, so a merged figure would be wrong on nearly every
bind. Every count stays a CpG count on the CpG axis; the axis is now indexed by batch.
`$coverage$sample_miss` is already per-sample and simply concatenates.

### 8.4 Cross-sample re-finalization

A finished record's `DNAmPhysAge` is a z-score against **that run's** cohort, so binding three
time-point records yields three within-batch z-scores in one column. The batch label is the answer
to that. **Never re-finalize by default** -- it would rewrite numbers the user has already seen and
standardize across batches that were deliberately imputed differently, which for the time-point user
destroys the object they built on purpose.

`refinalize_clocks()` is the opt-in form, and it is **exact rather than approximate**:
`physage_raws()` runs per-sample, so the per-sample intermediates retained on `$provenance$pending`
reproduce the single-pass number over the bound cohort. Keeping them there leaves the record's four
top-level elements intact. The chunk-reassembly user is the one who wants the re-finalize and the
batch label is what tells them apart from the time-point user; the routing keys on
`spec$cross_sample`, so no clock id is named, and the call is a no-op that says so when a record
deferred nothing.

**`rbind` says when it matters.** A bind that produces more than one batch *and* carries a
non-empty `pending` emits one `cli_inform` naming the affected columns and pointing at
`refinalize_clocks()`. The names come from `names(pending)`, which **is** the catalog's declared
`cross_sample` set arriving through `spec$cross_sample` -- so the message names clocks without a
clock list existing anywhere in `R/bind.R`. It stays an `inform`, not a `warn`: a per-batch
reduction is what batching forces, and for the time-point user it is the intended result. Single-
batch binds and ordinary clocks are silent. `Reduce(rbind, ...)` emits once per intermediate bind,
since each intermediate really is in that state; `do.call(rbind, ...)` emits once.

**Re-finalizing early is never wrong, only superseded.** `refinalize_clocks()` reads `pending` and
never consumes it, so it is a pure function of an accumulating store and composes in any order:
`refinalize(rbind(r3, refinalize(rbind(r1, r2))))` reproduces the single-pass column over all
three blocks (measured: 3.6e-15 absolute), `pending` grows 12 -> 18 across the interleaving, and a
second call on a settled record is a no-op. The mid-flight record between the two calls is a
genuine mix -- 12 re-finalized rows against 6 within-batch ones -- which is why the outer call is
still needed and why the message fires on the bind that created it. Do not "optimize" this by
clearing `pending` after a re-finalize: that is what makes the composition work.

Prerequisite this doc does not own: `[`, `cbind`, `as.data.frame.mc_result` and `augment()` do not
exist, so batch-label-survives-subset cannot be built or tested until they do.

## 9. Phase 7 -- `prep()`, the prepared-input record (not started)

`prep(DNAm, clocks, ...)` returns an S3 record, list-like and keyed by clock id, whose element for a
clock is **the matrix that clock's coefficients multiply**: its panel, imputed clock-wise and
panel-wise, with its normalization applied or not. `calc_clocks()` accepts either a `DNAm` or a
`prep` record, because `score_cohort()` already factors as `score(prep(...))`.

### 9.1 It is a plan, not a payload

Materialization is **lazy by default**: `p[["Horvath1"]]` builds one clock's matrix on demand and
`calc_clocks(p)` iterates with one resident at a time, so peak RSS is unchanged.
`materialize(p)` is the eager form.

Measured -- `sum |panel_i|` against the union:

| | `sum \|panel\|` | union | overlap | eager, n=1000 |
|---|---|---|---|---|
| 101 bundled clocks | 53,127 | 39,025 | 1.36x | 425 MB |
| 28 external (pack) clocks | 3,083,623 | 378,363 | **8.15x** | **24.7 GB** |

PCBrainAge alone declares 357,852 CpGs, the 13 SystemsAge organ clocks 125,175 each, the 14
PCClocks 78,464 each. Eager is free for the bundled set and fatal for the packs, so the record holds
a plan -- `DNAm` reference, the shared cohort-mean cache, the per-clock panel spec, the per-clock
normalization decision -- and materializes from it. The shared cache survives: the partial-NA
cohort-mean fill is clock-independent, so it is built once and every materialization draws from it.

### 9.2 There is no single prepared matrix

- **Imputation is per clock.** A fully absent CpG takes the clock's own vendored ref as a scoring
  CpG and the gold-standard mean as a normalization CpG. The same probe gets a different number
  depending on which clock is asking.
- **Normalization is per clock.** Different clocks declare different panels and targets, so one
  normalized matrix cannot serve two of them.

The line that does exist is cohort-global vs clock-specific: cohort means are hoisted (sec 1.1),
everything else stays per clock.

### 9.3 Normalization is an argument, and its half shipped

`normalize=` is an ordinary named logical on `calc_clocks()` / `sim_DNAm()`, resolved once by
`resolve_normalize()` before any DNAm is read. `prep()` takes the same argument when it lands. A
bare scalar is a policy that reaches only the clocks able to honor it; a named logical is a
per-clock claim whose unhonorable entries error. There is no per-clock options list.

Declared state: `normalization` is non-`none` on 4 clocks -- `bmiq` on Horvath1 and Knight, `noob`
on Horvath2, `quantile` on DunedinPACE. Both BMIQ clocks ship a `bmiq_gold_standard` probe_set
(`goldstandard2.csv.gz`, 21,368 probes, `identical()` between the two groups, each scoring panel a
strict subset). `noob` is unreachable from a beta matrix under any design.

Read constitutive-ness off the **scheme**, not a per-clock field: `quantile` is constitutive
(DunedinPACE's QN cannot be declined), `bmiq` is preprocessing and defaults **off**.
`clock_norm_target()` is the scheme-agnostic target accessor and the branch is normalize-then-linear
dispatching on `clock_norm_scheme()` -- the existing **pre-transform** branch kind, not a new entry
in the closed branch set.

Where BMIQ departs from QN: a fully absent norm CpG is **dropped**, not vendor-filled. QN fills from
the gold mean and normalizes the full panel; BMIQ estimates the sample's own mixture from the panel,
so filling with target-drawn values biases the fit toward the gold standard. The thin-background
warning's fate text is therefore scheme-aware.

Normalization costs nothing in chunk-safety: `quantile_norm()` is per-sample against a fixed target
and bit-identical on a row subset (measured, max abs diff 0); BMIQ fits within a sample. Neither
adds to `cross_sample_at`.

### 9.4 Open

- What `p[["Horvath1"]]` returns for a normalizing clock: the normalized **scoring** panel, with the
  normalization panel transient inside materialization. The 21k panel likely deserves a separate
  accessor rather than being the default element.
- **Deduplicating normalization work.** BMIQ is a per-sample mixture fit and far more expensive than
  rank mapping. Horvath1 and Knight share one 21,368 panel with `identical()` targets, so a run
  normalizing both fits the mixture twice per sample. `resolve_cpgs()` already dedups the *panels*
  into `norm_parts` with a per-clock index -- the missing piece is caching the normalized output in
  the same place. Confirmed safe to memoize, not done. It matters more under chunking, where it
  would otherwise repeat every block.
- **Row-level failure has no vocabulary yet.** BMIQ runs `on.sample.error = "continue"`, so a sample
  whose mixture fit fails scores `NA` while the cohort completes -- and shows `coverage = 1.00` on
  both panels with a zero `sample_miss`, because coverage is computed upstream of scoring and is
  correctly reporting that every probe was present. Missing pheno, excessive NA and a failed
  calibration are indistinguishable after the fact.
  `betanorm::bmiq_calibration()` returns `$success` and `$failures` (`sample_index`, `sample_name`,
  `stage`, `message`); `bmiq_panel()` keeps only `$calibrated`. A `samples_coverage()` reason column
  is the intended fix.
- **Not a gap:** BMIQ failing on `sim_DNAm()` output is expected. `random_betas()` is U(0,1), which
  has no 3-state structure. Tests needing a normalizing clock must build methylation-shaped input --
  `test-normalize.R` jitters the gold standard.

## 10. Tests

Always-on, no new dependency -- the in-memory source is what makes this possible and is the whole
reason it exists. Every numeric comparison is `expect_equal`; `expect_identical` is banned
package-wide.

Landed in `test-chunk-invariance.R`. One helper builds a cohort carrying all three missingness
shapes at once (ordinary partial NA, partial-but-all-NA-in-block-1, fully absent) and every test
reuses it:

- **Chunk invariance** (sec 4.2): one cohort scored whole vs in three blocks, through the seam
  directly. Both paths end in the same `finalize_cross_sample()`, over **every** clock in the
  sequence with no exclusion. The request is asserted to contain a cohort-reducing clock, derived
  rather than named.
- **The loop defers exactly the declared set**: `names(pending)` equals `spec$cross_sample`, those
  clocks leave the loop with a wider-than-one-column per-sample intermediate and no score, and the
  finalize turns each into an ordinary `n x 1` column.
- **The sample-axis split is derived**: `split_cross_sample()` partitions `resolve_clocks("all")`
  totally, and each sex-routed alias matches its members.
- **The block-dependency case**: a column partial cohort-wide but all-NA within one block lands in
  `partial_fill` and counts as `imputed_partial`, not `imputed_full`.
- **Coverage assembles by concatenate, never by sum.**
- **The clocks gate throws out of `mc_cohort()`**, before anything is scored.

Landed in `test-value-gates.R`: `+/-Inf` stops with and without an unrelated NA present, the abort
names the sample and CpG, each range flag warns on its own, an M-value matrix trips both while still
returning scores, ordinary betas and all-NA columns pass in silence, and the kernel bails on `Inf`
with `stats` unset.

Landed in `test-bind.R` (Phase 4, sec 8): one test per gate (colliding ids, differing columns,
differing `pheno_id`, differing `normalized` -- all throw; disjoint records bind); the folds agree
(`rbind(r1, r2, r3)`, `Reduce` and `do.call` produce the same record); the label is a function of
the sample ids and nothing else, argument names are refused, and re-association returns an
`identical()` record with each sample still under the batch it was scored in; a run given no
`pheno` still carries the id column and binds against one that was given a pheno;
`clocks_coverage()` is one row per (clock, batch) and `samples_coverage()` carries each
sample once under its own batch; a probe all-NA in one batch is recorded there rather than merged
away; and an opt-in `refinalize_clocks()` over retained `pending` reproduces the single-pass
cross-sample column. That last one is the only numeric test in the file -- it derives its golden
from a whole-cohort run of the same matrix, so it needs no fixture and no seed.

Parked with Phase 6, listed so the park is reversible:

- **The headline test**: `calc_clocks_chunked(matrix, budget)` equals `calc_clocks(matrix)` over the
  whole record. Needs no Suggests.
- **Gates fire before pass 2**: via a source that records how many times each block was read, and
  that each block is read exactly twice.
- **Warnings fire once, not per block**: the three range/Inf flags over a multi-block run.
- **Cross-block diagnostics**: dead rows, duplicate ids and missing pheno ids name the cohort's set,
  not one block's.
- **Reader canonicalization**: the null-spelling set resolves every spelling in sec 5.1's table,
  `Inf` survives, a stray token errors naming column and token, and orientation is detected from the
  header on both a tall and a sample-major file.

Skipped when Suggests are absent, per existing practice: the duckdb reader tests, and the
`DelayedArray` source including one transposed input and one hostile chunk geometry.

For Phase 7: `calc_clocks(DNAm, ...)` and `calc_clocks(prep(DNAm, ...))` agree; a lazy `prep` and a
`materialize()`d one agree; a clock whose normalization is declined scores what it scores today.

Parity is untouched -- it scores complete cohorts and is not part of "run the tests".

## 11. Non-goals

- No random-access source contract beyond what sec 5.3 uses. Two scans is the whole requirement.
- No pre-imputation of a block outside the scorer (sec 1.1).
- No single prepared `DNAm` out of `prep()` (sec 9.2).
- No eager materialization by default, and no whole-array normalization inside `calc_clocks()`.
- No chunk-size autotuning. Take the budget from the caller and say what was used.
- No parallelism in the first cut. Chunking trades speed for RSS by design, and the run is
  deflate-bound and serial anyway (sec 5).
- No float32 anywhere on the input path.
- No `data.table` dependency -- duckdb owns file reading.
- No hand-rolled on-disk format, and no second one: the store is a duckdb file (sec 5.2).
- No `rbind` across cohorts *reconciling* anything -- no re-imputation, no merged coverage
  denominators (they nest, sec 8.3). It binds and labels. Re-z-scoring is the one exception and it
  is not silent: opt-in only, exact rather than approximate, and keyed on `spec$cross_sample`
  (sec 8.4).
- No new identity key. Ids are the identity; a batch **label** is not a `batch_set_id`.
