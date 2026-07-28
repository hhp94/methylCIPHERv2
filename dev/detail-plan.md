# methylCIPHERv2 Rewrite -- Detail Plan

Canonical long-form design for the rewrite. Overview: [`migration-plan.md`](migration-plan.md).
Decisions log: [`DECISIONS.md`](DECISIONS.md).

This document folds the architecture discussion (engines, packs, memory, coverage) and the
operational contracts that used to live only in the overview plan. Prefer updating **here**
when behavior is specified; keep the overview short.

**Current truth only.** Superseded designs are removed from these sections, not annotated inline;
their history lives solely in [`DECISIONS.md`](DECISIONS.md).

---

## 1. Product and public API

### 1.1 Exported surface

| Function | Purpose |
|---|---|
| `list_clocks()` | Discover / filter (format, covariates, batch_dependent, external, ...) |
| `get_clock(id)` | Lightweight metadata for one id |
| `get_clock_probes(id)` | Required probes + imputation refs |
| `calc_clocks(DNAm, clocks, pheno = NULL, ...)` | Main scorer |
| `clocks_coverage(x)` | Per-clock coverage data.frame (wide) from the result record |
| `samples_coverage(x)` | Per-(sample, clock, panel) coverage data.frame (long) from the result record |
| `report(DNAm = NULL, result = NULL, ...)` | One QC report over the DNAm input, an `mc_result`, or both; writes a PDF (sec 4.3) |
| `augment(scores, data, ...)` | Join scores to analysis / pheno tables |
| `clear_mc_cache()` | Report cached external packs; consent-gated removal (never auto-deletes) |
| `mc_data_download()` | Explicit pre-fetch (cache dir overridden by `options(mc.cache_dir=)`) |

### 1.2 Removed or non-primary surface

- Per-clock exported `calc*` as the **implementation** (optional thin legacy wrappers OK).
- Exported `MiAge_*` helpers, `calcUserClocks`.
- Return types that flip between vector and pheno-appended data.frame as the default.
- Row-order phenotype alignment.
- `library()` / `require()` inside scoring paths.

### 1.3 Result contract

`calc_clocks()` returns an **S3 record over `list`**, class `"mc_result"`:

```r
structure(
  list(
    scores     = <n x k double matrix>,   # dimnames = samples x clock ids
    pheno      = <aligned pheno, or NULL>,# id column + required covariates only (sec 5.1)
    coverage   = <per-column coverage>,   # see sec 4
    provenance = <clocks, requested, dependencies, covariates_used, pheno_id>
  ),
  class = "mc_result"
)
```

Inherits `list`, **not** `matrix`: a matrix subclass drops class + attributes on the first
`x[, "Horvath"]` / `t()` / arithmetic, silently discarding coverage and provenance. Verbs are
methods so no operation loses data:

| Method | Behavior |
|---|---|
| `print` | dims, then a labelled preview per component (`<$pheno>`, `<$scores>`, ...) so the record reads as the list it is; never a full matrix dump |
| `as.matrix` | `$scores` -- the naked-numbers escape hatch |
| `as.data.frame` | scores as a data.frame (sample id column + score columns). **Never** `$pheno` -- see sec 5.1 |
| `[` | subset rows/cols of `$scores` **and** the matching pheno/coverage/provenance -> `mc_result` |
| `cbind` | bind score columns; require equal `sample_id` sets (read straight off `$provenance$sample_id`) |
| `rbind` | **refuses** -- stacking samples leaves any cohort-dependent score computed against the wrong cohort. Score per cohort, bind the `as.data.frame()` outputs |
| `augment` | analysis-ready table: `as.data.frame()` + the record's aligned covariates (`$pheno`) + an optional user `data` joined by sample id. Plain exported function, not a `broom` generic (avoids the clash). Built. |

Built so far: `print`, `as.matrix`, `as.data.frame`, `augment`. Still unbuilt: `[`, `cbind`.

Not methods: `clocks_coverage()` / `samples_coverage()` format `$coverage` (never re-touch beta) and
replace the `summary` that earlier drafts listed; `codebook()` and `bibliography()` are plain
functions (they take a result / clock tokens / `"all"`), since `utils::citation` is not an S3 generic.
`bibliography()` reads the **sync-vendored** `inst/bibliography/clocks.bib` (written by `sync.R`'s
`vendor_bibliography()`, keyed by `bib_key`) and emits full BibTeX (`format = "bibtex"`) or an enriched
data.frame (citation / title / journal / year / doi / pmid / url), falling back to a key+PMID stub for
any catalog `bib_key` absent from the .bib. `codebook()`'s training-population fields still await the
`master_source_of_truth.csv` sync. See DECISIONS 2026-07-23, 2026-07-24, 2026-07-27.

Rules:

- One clock -> still `n x 1` scores, never a bare vector.
- Scores only; no automatic pheno columns. This governs `$scores` and the `as.data.frame()`
  shape. `$pheno` is a separate component holding the aligned, narrowed covariates (sec 5.1) --
  it is what `augment()` and `cbind` gate 2 read, and keeping it off `as.data.frame()` is what
  keeps subject data from riding the normal export path.
- Multi-output / pack requests -> multiple columns.
- **Dependency clocks are returned as columns**, after the requested ones (§2.1). They were
  genuinely computed and each carries its own coverage row, so dropping them would hide the state
  of a composite's inputs — a `DNAmFitAge` resting on a `GrimAgeV1` scored from few CpGs must be
  inspectable. Column order is requested-first (request order), then auto-added deps (compute
  order); `$provenance$requested` / `$dependencies` partition `$provenance$clocks`. Callers wanting
  only what they asked for subset downstream.
- `augment()` is the join step; optional `augment = TRUE` sugar may exist later (non-default).

### 1.4 Legacy `calc*` (collaborators)

If collaborators need old names:

- Generate thin wrappers -> `calc_clocks()` + optional legacy return shape (append pheno).
- Do **not** reimplement old bodies from `dev/legacy/R-pre-rewrite.R`.
- Deprecate for one release, or keep as permanent soft compatibility -- product choice in
  `DECISIONS.md` when locked.
- Fixtures and tracker always target `calc_clocks` / catalog ids.

---

## 2. Runtime architecture

### 2.1 Mental model

```text
calc_clocks
  |- normalize / expand aliases (packs)
  |- prepare_inputs (once)
  |    check_dnam, align_pheno by ID, covariates_required
  |    full-panel warn if needed (localized)
  |- compute_coverage (once, keyed by clock id) -- needs no score
  |- for each work unit -> score_*  (engine or pack)
  |    returns just the score matrix
  \- construct_mc_result -> record (scores + coverage + provenance)
clocks_coverage(x) / samples_coverage(x)  # pure reads of x$coverage
```

**Resolve + prepare-once front end** (in `R/resolve_inputs.R`). Two phases, both *before* any
scoring, both run exactly once — never inside the per-clock loop:

- `resolve_clocks(clocks)` — **pure namespace resolution**: user tokens -> catalog `clock_id`s.
  Accepts (any mix) the alias `"all"`, a **keyword** from `MC_TAGS` (`gestational`, `mitotic`,
  `mortality`), a `group_id` (expanded to its members), or a bare
  `clock_id`, resolved by a fixed precedence, **maximal clock-id set first**: `"all"` > keyword >
  `group_id` (whole group) > bare `clock_id`. A keyword expands to tokens that resolve by the same
  group/clock rules — one level only, so a keyword never names another keyword and there is no
  cycle to chase; a keyword naming a token the catalog no longer has is a **hard stop**, not a
  silently shorter set. Keyword names are lowercase and clock/group ids are not, so the namespaces
  cannot collide. When a token is **both** a `group_id` and a `clock_id` (e.g.
  `"DNAmFitAge"` names both the 7-member group and one composite member), the **group wins** so a
  group request returns as many member ids as possible — one rule covering a non-clock group_id
  (`GrimAge`, `SystemsAge`), a group_id that is also a member id (`DNAmFitAge`), and a singleton
  (`group_id == sole clock_id`, `MiAge`). Returns unique ids in first-seen order. Unknown tokens
  are an error naming them. It does *no* math and *no* ordering, and does **not** resolve
  dependencies — that is `resolve_clocks_sequence`'s job (next bullet). **Input contract
  (enforced):** `clocks` must be a non-empty character vector; `NULL` / missing / `character(0)` /
  `NA` / empty-string tokens error at the front door (callers that want everything pass `"all"`), so
  an empty request can never silently resolve to nothing.
- `resolve_clocks_sequence(clock_ids)` — **the dependency plan step** (settles the mechanism §2.1
  used to leave open; see DECISIONS 2026-07-17). Takes the resolved set and returns its **transitive
  closure** over `clock_inputs` (the recipe's `inputs`, or a sex-routed alias's two members; read
  via `clock_depends_on()`), **topologically sorted** so a
  clock always comes after everything it depends on (cycle-guarded; stable first-seen otherwise).
  Cross-clock deps DO exist and are declared in the catalog — only 3 clocks: `GrimAgeV1/V2` -> their
  CpG surrogates; `DNAmFitAge` -> gait/grip/VO2max + `GrimAgeV1`, which **crosses group boundaries**.
  Like the covariate union, deps are pulled in **on demand** (a request that never touches them adds
  nothing). **Auto-added ids** (in the plan, not in the requested set) are scored *and returned* as
  columns (§1.3), tracked apart in `$provenance$dependencies`. This is why the prepare-once unions
  below run over the **plan**, not the raw request.
- prepare-once, **parametrized by the resolved set, not per clock**:
  - sample identity is `rownames(DNAm)`, mandatory: `check_DNAm()` hard-errors identity-less DNAm
    and hands the caller the one-liner to name anonymous rows (§5.1). The package never manufactures
    ids.
  - `check_DNAm(DNAm)` — *universal* invariants only: double matrix; unique CpG colnames; unique
    non-NULL sample-id rownames; a `^cg`-prefix orientation
    warning (the `nrow > ncol` dimensional guess was dropped -- it false-positives on small
    panels). Carries no clock knowledge.
  - covariate check is a **union over the compute plan** (requested + auto-added deps):
    `extra_columns <- unique(unlist(lapply(plan, clock_covariates_required)))`, fed once to
    `check_pheno(pheno, ID, extra_columns)`. Unioning over `plan` (not the raw `clock_ids`) is
    load-bearing: a dep can require a covariate the request does not — the single `DNAmFitAge` clock
    declares only `Female`, but its `GrimAgeV1` / `DNAmVO2max` deps also need `Age`. Requirements are
    a per-clock catalog field (`covariates_required`), so requesting a single component vs. a whole
    group is automatically correct — **no per-clock/per-group check registry, nothing hand-coded per
    clock.**
  - **missing covariates** — an `NA` in a required covariate is legal and propagates: that sample's
    score is `NA` for every clock consuming the covariate, the row is never dropped, and no error is
    raised. `check_pheno()` warns **once**, naming each affected column and its NA count, scoped to
    rows that survive the id-join (an `NA` on an unscored cohort row does not warn). Clocks that do
    not declare the covariate still score finite for that sample.
  - **coverage floors** — `calc_clocks()` carries two independent floors, named for the axis each
    reads, and they differ in *severity*, not just in axis:

    | Floor | Axis | When | Action |
    |---|---|---|---|
    | `min_clocks_coverage` (0.75) | columns -- which of a clock's CpGs exist in `DNAm` at all | before scoring | **stop** (graded, see below) |
    | `min_samples_coverage` (0.75) | rows -- how much of that panel a given sample observes | after scoring | **warn**, always |

    A missing column is missing for *every* sample alike -- it is a property of the input matrix,
    so no sample escapes it and the call stops. A sparse row is one bad sample among good ones:
    NA-ing it would silently discard data the caller may still want, so the score stands and the
    warning names it. **No branch NA-s a row for coverage.**
  - **row floor** — `check_row_coverage(coverage, min_samples_coverage)`, run once over the hoisted
    coverage structure (`compute_coverage()` already masked routed members): a sample's observed
    fraction is `(score_present - sample_miss) / score_needed`, or the `norm_*` pair for a clock
    that normalizes (DunedinPACE), since that is the panel `count_sample_miss()` counted over. No
    branch takes the floor as an argument and no branch re-implements the check -- it is one
    function (`row_coverage()`) reading the per-clock record + its `sample_miss` vector. It runs
    over **every clock computed**, not just the returned columns: clocks with a `NULL` record
    (sex-routed aliases) are skipped, and the routed members that carry their panels are checked
    under their own ids even though they have no column. `calc_clocks()`
    validates the floor up front so a bad argument does not cost a full scoring pass.
  - **column floor** — `check_coverage(cpg_list, min_clocks_coverage = 0.75)`, run once on the
    `resolve_cpgs()` skeleton (§2.3a). **Graded**, on a clock's **scoring** panel (`score_present /
    score_needed`) -- the only panel that can stop the call:

    | Band | Action |
    |---|---|
    | 0 observed CpGs, or under `min_clocks_coverage` | **stop**, naming each clock + counts + percentage |
    | within 10% above `min_clocks_coverage` (`< 1.1 x min_clocks_coverage`, capped at 1) | **warn**, same listing |
    | otherwise | silent |

    A panel at **0 observed CpGs always stops**, regardless of `min_clocks_coverage`. So
    `min_clocks_coverage = 0` relaxes the quality gate but never buys a score computed from nothing --
    `Horvath1` on a foreign panel used to return `anti_trafo(0)`, a plausible ~40 years, for every
    sample. Full coverage never warns, whatever the threshold. The warn band is *relative* so it
    tracks a caller who moves the floor, rather than pinning an absolute second number. Clocks with
    no panel at all are not gated here; they fail later on the unimplemented-scorer branch, which
    is the clearer error.
    A clock's **normalization background** (only DunedinPACE carries one) is a *separate,
    warn-only* check in the same function: under `min_clocks_coverage` it warns and never stops,
    because absent gold CpGs fill to their reference mean before quantile-normalizing, so the
    clock still scores. "Cannot score" and "normalized against a thin background" are different
    failures and do not share a band.
    **Sample-invariant** — it reads set sizes, so it fires on panel/array mismatch (wrong array, a
    group's weights off-manifest, a simulated panel built from too narrow a clock set), never on
    per-sample NA, which is tier 2's job (§4.2). Running it over the **plan** is what makes a
    collapsed *input* visible: a `GrimAgeV1` at 0 present CpGs is named even when only
    `DNAmFitAge` was requested. The listing caps at 10 clocks -- mass failure is near-always one
    root cause, so a 50-line dump would bury it. It lives in `calc_clocks()`, not inside
    `resolve_cpgs()`, which stays pure set math (§2.3a), and before `build_partial_cache()`, so
    nothing expensive runs ahead of it.

Deliberately **not** in this front end: the full-panel warning (it is the `sample_scale`
transform's own precondition — §2.4), imputation checks (they live with the impute step — §2.3),
and array-normalization (never executed — §2.4a). Same principle throughout: a clock-specific
precondition lives with the step that consumes it, not hoisted into the universal prepare path.

### 2.2 Dispatch: linear engine + a finite branch set

There is no recipe interpreter -- not even as a reference/fallback. Every work unit routes on
the catalog pair `(weights_format, computation_type)` to one of a small, closed set of branches.
Meta `weights_format` / recipes still **define** the math (and fixtures check it); the branches
are pragmatic, hand-optimized implementations of that contract.

| Bucket | Implementation | Typical content |
|---|---|---|
| Linear engine | Shared `linear_score()` | Most `cpg_coefficient` clocks |
| Pre-transforms | Small modules feeding linear | Zhang `sample_scale` (stats on full panel -> apply to coef subset) |
| Family orchestrators | Dedicated internal functions | **GrimAge** (`score_GrimAge`); **Dunedin** (`score_Dunedin`); **EpiTOC2** (`score_EpiTOC2`); **DNAmFitAge** composite (`score_DNAmFitAge`) |
| Sex-routed aliases | `score_sex_routed()` -- selects a member per sample | The 7 DNAmFitAge stems |
| Batched pack scorers | One shared subset + one matmul per group | **PCClocks**, **SystemsAge** packs |
| Custom | Dedicated helper | **MiAge** (`score_MiAge`) |

Bundled packs may own their orchestration (shared intermediates, multi-column assembly) but call
the shared `linear_score()` / impute helper for every linear sub-step, so imputation lives in
exactly one place (sec 2.3).

Every branch, whatever its shape, ends in the same three-field record, so the boilerplate that
builds it lives once in [`R/utils.R`](../R/utils.R): `cached_cols()` (which present CpGs the
cohort-mean cache covers), `count_sample_miss()` (per-sample fill counts feeding
`coverage$score_imputed_partial`), and `score_matrix()` (the n x 1 `dimnames = list(sample_id, id)`
matrix). A new branch should reuse these rather than re-inline them -- they are the shape contract,
not a convenience.

**Batched pack scorers** (`score_pack_group()` in `R/score_pack.R`) are the exception, for the two
external groups whose members share one large CpG panel (PCClocks ~78k, SystemsAge ~125k). Scoring
those per clock repeats the panel-wide `DNAm[, present]` subset once per member -- the dominating
cost (13.7x measured for 14 PC clocks; the matmul itself is negligible). So a group is scored in a
single `pack_design()` (one subset, reused) + one `pack_linpred()` matmul over all requested
columns, with per-clock intercept, covariate, and output transform applied to the resulting
columns. This is a *batched* linear kernel, not `linear_score()` -- but it reproduces the identical
imputation contract (partial-NA -> cohort cache; absent -> vendor-mean offset from `pack$impute`),
so imputation semantics still live in one described place. It is **only** worth it where the panel
is huge and shared: the bundled clocks (union ~20.7k CpGs, 1.6x overlap) score all 86 in ~60ms, so
they stay on the per-clock engine -- no general "group linear clocks" layer (see DECISIONS
2026-07-21). Requesting a single pack member returns just that member (no family expansion); the
batch simply collapses however many members are in the plan into one matmul.

### 2.3 Linear engine

```r
linear_score(DNAm, coefs, intercept, impute_spec, transform = NULL)
```

Imputation / missingness -- **two missing-kinds, two sources, never crossed** (implement in one
place):

| Case | Source | Behavior |
|---|---|---|
| **Partial NA** on a probe **present** in `DNAm` | current cohort | Fill from the cohort (mean over observed samples for that probe). Cohort-dependent by design; a **fallback** for users who did not pre-handle NAs |
| **Completely absent** probe (not in `DNAm`), policy fill | vendor | Fill from the clock's vendored ref (mean / median / sex-wise, per meta), then weight |
| Completely absent, policy `drop` (default) | -- | Drop the term (zero contribution) |
| Sex-keyed vendor ref | vendor | Require `Female` before fill; pick female/male ref vector |

The split is load-bearing: vendored-filling a partial NA injects a foreign cohort's data into the
user's samples; cohort-filling a fully-absent probe is impossible (no data to borrow). `mean` is
the current fallback statistic, not a fixed rule. Never apply one package-wide impute policy to
all clocks -- the vendor spec comes from meta (`imputation`, recipe steps).

### 2.3a Partial-NA is a shared cohort cache, not per-scorer work

The partial-NA (cohort) fill is **precomputed once per call and shared**, not redone inside each
scorer. Component families (DNAmFitAge's 7 members, GrimAge surrogates) reuse the same CpGs, so a
shared missing probe would otherwise get its cohort mean recomputed once per clock; hoisting it also
makes every scorer read the **identical** filled column for a shared CpG (a consistency guarantee,
not just a speed win). Prepare-side pipeline, gated so clean betas pay nothing:

- **Global gate.** `anyNA(DNAm)` is monotone: no NA anywhere ⇒ no clock can have a partial NA ⇒ skip
  the whole cache + per-sample layer. This is the common case (pre-handled betas).
- **One classifying pass.** `slideimp::mat_miss(DNAm, col = TRUE)` (per-column NA counts, no logical
  mask) splits every column three ways: `0` clean, `1..nrow-1` partial, `nrow` **all-NA**. `mat_miss`
  keeps the pass allocation-free.
- **All-NA present column → reclassify ABSENT.** Zero observed values ⇒ no cohort mean exists ⇒ it is
  informationally identical to a probe off the panel. Drop it from the **usable** column universe
  (`setdiff(colnames(DNAm), all_na_cols)` — a name-set op, **not** a `DNAm[, keep]` copy) so it flows
  to the vendor/drop path and is counted in `n_cpg_score_miss` (§4). This is why `mean_imp_col`'s
  own "0 observed → unchanged" branch never has to fire in our path.
- **Empty sample (all-NA row) → HARD ERROR.** `mat_miss(DNAm, col = FALSE)`; `row_miss == ncol`
  cannot be scored. Load-bearing reason it errors rather than warns: cohort-mean imputation would
  otherwise fill its every cached cell with the column mean and emit a plausible-looking
  **fabricated** score instead of an honest `NA`. (An all-NA row is `isnan`-skipped in every column
  mean, so it never biases other samples; the damage is confined to its own fake score, but that is
  reason enough to refuse it. Like the mandatory-rownames rule in §5.1, this is a hard error with no
  permissive mode; the honest output is `NA`, not a fabricated number.)
- **Cache build.** `cache_cpgs = intersect(present_needed_union, partial_na_cols)` — present (a mean
  exists), needed by some clock (don't cache dead columns), and actually partial. Subset **first**,
  then impute: `slideimp::mean_imp_col(DNAm[, cache_cpgs])`. `mean_imp_col` returns a matrix the same
  width as its input (untouched columns are `memcpy`'d), so calling it on the full panel would
  allocate a second `n × p` copy (the §3 spike); narrowing columns first bounds the cache to the
  handful of partial-NA needed probes (empty when clean). The cache is an `n × k` matrix scorers read
  from; raw `DNAm` is never mutated.

Front end split (which function owns what): `scan_missing_cpgs()` (numeric — the gate, the classify,
the empty-row throw, `usable_cols`) and `build_partial_cache()` live with the impute machinery
(`R/impute_DNAm.R`); `resolve_cpgs()` (pure set math over `usable_cols`, §4 aggregate skeleton +
`present_needed_union`) lives with the other resolvers (`R/resolve_inputs.R`). Order in
`calc_clocks()`: `scan_missing_cpgs → resolve_cpgs → build_partial_cache` (the cache needs
`resolve_cpgs`'s present union; `resolve_cpgs` needs `scan`'s `usable_cols`).

### 2.4 Zhang / `sample_scale` transform

Do **not** scale and keep a full `n x p` matrix.

```text
For each sample i:
  mu_i, sigma_i <- moments over ALL probes in DNAm[i, ]
Keep only EN / coef probes j:
  z_{ij} <- (DNAm[i,j] - mu_i) / sigma_i
linear_score(Z, coefs, policy = drop)
```

- Moments need full width -> user should pass full beta (450k/EPIC).
- Scoring only needs the small coef submatrix after stats.
- This is a **transform module** in front of the linear engine, not a separate scientific
  family.

**Full-panel notice (Zhang2019-only, hard-coded):** do **not** fire a generic
`needs_full_probe_panel` / `ncol(DNAm) < 1e5` warning across flagged clocks. Only `Zhang2019`
needs this today, so hard-code the exception: when the worklist includes `Zhang2019`, emit a
`message()` (not `warning()`) — Zhang2019's original code computes per-sample moments over **all**
CpGs, but a large-enough subset is usually sufficient. Fired here in the transform, not in the
prepare front end.

There is no `batch_ops` catalog marker any more, and `sample_scale` was never a cross-sample op:
it is a **within**-sample z-score over one sample's own row, so Zhang2019 is `per_sample` and
chunk-safe (`cross_sample_at` is `NA`; the retired `extract_batch_ops()` wrongly lumped it with
`cohort_zscore`). Nothing keys the transform case off the catalog: `clock_needs_full_panel()` and
`check_DNAm_extra()` are **gone**, the generic machinery having collapsed into
`resolve_DNAm_extra()` with one `"Zhang2019" %in% clock_ids` special case, as no other clock
exercises it.

### 2.4a Normalization policy — annotate, never execute

`sample_scale` above is a *scoring-recipe* transform, not array normalization. Array
normalization proper is the per-clock `normalization` field: `none` ×104, `BMIQ` ×7,
`quantile` ×1 (DunedinPACE), `noob` ×1 (Horvath2). **The package executes none of it, except
DunedinPACE's `quantile`** (see below).

- BMIQ / noob are squarely **upstream** (sesame / minfi) — the user's responsibility.
- Horvath's BMIQ-to-golden-mean is **deliberately skipped**: a correctness bug in RPMM. Parity
  fixtures show ~0.9999 correlation of no-BMIQ vs the Horvath server, so re-implementing it buys
  nothing and inherits the bug. Users who want BMIQ run the RPMM pipeline themselves first
  (most won't).
- DunedinPACE's `quantile` is the **one executed** normalization: the `score_Dunedin` branch
  quantile-normalizes the gold panel to `gold_standard_means` via `betanorm::quantile_norm`
  (bit-exact with the author's `preprocessCore`) before the linear score. It is intrinsic to the
  clock (the score is defined on normalized betas), not array prep, so it cannot be pushed upstream.

So `normalization` is a **coverage / provenance annotation** — surfaced via a (future)
`clock_norm_scheme()` accessor feeding `norm_needed` (§4) — **not a compute step and not a
check**. There is no shared "norm intermediate" layer: nothing is jointly normalized by us
(`sample_scale` is per-sample and Zhang-only today; QN is external), so building one would be
speculative. The only norm-adjacent thing the package runs is the `sample_scale` transform.

### 2.5 GrimAge pack (orchestrator)

**User-facing:**

- Prefer `"GrimAge"` and/or **component** ids (`DNAmADM`, `DNAmPACKYRS`, ...).
- Do not push `GrimAgeV1` / `GrimAgeV2` as the primary UX (catalog may still retain those ids
  for fixtures and provenance).

**Under the hood (product policy):**

- Always compute enough shared work for both V1 and V2 pipelines (fast enough).
- Component columns come from the **V2** path (not duplicate V1 surrogates).
- Full pack return includes V2 components + GrimAgeV1 + GrimAgeV2 columns (exact set TBD in
  implementation; fixtures pin catalog ids).

Too special to force through a generic path -> `score_GrimAge()`.

### 2.6 SystemsAge pack

External asset bundle carrying `$organs`/`$systems`/`$age`/`$impute` matrices + a small
systems_PCA tensor tree; loaded **once** upfront and threaded to every member (no per-member
reload). The 11 organ sub-clocks are plain `cpg_coefficient` linear (coef from `$organs`, shared
engine). The two component-matrices composites (`Age_prediction`, `SystemsAge`) are the
`score_systemsage_group()` family orchestrator: age-linear front -> quadratic, and for the overall index
11 raw system predictors + poly-scaled age -> center/scale -> systems_PCA project -> linear head.
Pipeline *shape* is hard-coded (no recipe walker); constants come from the catalog recipe via
`systemsage_*` accessors. Each member is scored independently on the shared linear kernel -- the
composite recomputes its raw system vectors (cheap; no shared-intermediate cache).

### 2.7 Worked call

```r
calc_clocks(DNAm, c("Zhang2019", "GrimAge"), pheno = pheno)
```

1. Resolve: Zhang2019 -> linear+sample_scale transform; `"GrimAge"` -> pack orchestrator.
2. Prepare: check DNAm; if Zhang in list and `p < 1e5` -> warn; align pheno; require Age +
   Female for GrimAge pack.
3. Zhang: full-matrix row moments -> EN subset -> scale -> linear(drop) -> 1 column + coverage.
4. GrimAge: orchestrator -> multi-column matrix + per-column coverage.
5. Assemble -> `mc_result` record; `clocks_coverage(out)` / `samples_coverage(out)` without
   re-touching DNAm.

---

## 3. Memory and copies

| Do | Don't |
|---|---|
| Pass **raw** `DNAm` into scorers (no upfront global intersect) | Build union-of-all-clocks submatrix in `calc_clocks` before dispatch |
| Let each scorer `DNAm[, probes]` -> small copy | Scale full DNAm in place for Zhang |
| Pack orchestrators may extract **one** union for the pack | Micro-optimize column views / sparse before profiling |

R: passing `DNAm` does not copy; `DNAm[, j]` does allocate the subset -- expected and fine for
clock-sized `j`.

---

## 4. Coverage: `clocks_coverage()` and `samples_coverage()`

Scorers **must** return coverage; assembly stores it on `x$coverage`. Coverage has **two tiers** at
two granularities, both recorded at score time (never re-intersecting beta):

### 4.1 Tier 1 -- per-clock aggregate (`clocks_coverage()`)

Recorded **per role** (a clock's scoring CpG set and its normalization/background set are different
sizes) and splitting the two impute sources. **Sample-invariant** — these are set sizes against the
usable column universe, the same for every sample:

| Field | Meaning |
|---|---|
| `clock_id` | Score column / catalog id |
| `norm_needed` / `norm_present` | Normalization / background panel (Zhang, QN clocks); `NA` when none |
| `score_needed` | Scoring CpGs (materialized `probe_sets[role=scoring]`; external groups: the pack's `$cpgs`) |
| `score_present` | Scoring CpGs found in `usable_cols` (colnames minus all-NA, §2.3a) |
| `score_used` | Terms that entered the sum |
| `score_imputed_partial` | Present-but-NA cells filled from the **cohort** cache |
| `score_imputed_full` | Absent probes filled from the **vendor** ref |
| `score_dropped` | Absent probes dropped by policy |
| `missing_cpgs` | Character vector of absent probes (incl. all-NA columns reclassified per §2.3a) |
| `policy` | Impute policy used |

```r
clocks_coverage(x)  # data.frame, one row per clock computed; NA where a stage does not apply
```

One row per clock **computed** (not per score column): a `role` column tags each row `returned`
(any score column, including a sex-routed alias with an all-NA panel row, and dependency columns) or
`routing_target` (a member kept for coverage, never a column -- it carries the per-sex denominator).
The counts are just `length()` of the present/absent sets `resolve_cpgs()` (§2.3a) already computed,
so they ride the skeleton to assembly and `clocks_coverage()` only **formats**; it never recomputes.

- **Free** if recorded at score time; `clocks_coverage()` must not re-intersect beta.
- Zhang: `score_needed` = EN coef count (not full array); `norm_needed` = panel size.
- Messaging that implies confidence must not ship without coverage context.

The `score_present / score_needed` ratio is also read at **prepare** time by the coverage floor
(§2.1), which stops before any scoring happens rather than waiting for the user to call
`clocks_coverage()`. Same numbers, two moments: the floor is the gate, `clocks_coverage()` the full
account of what got through it.

**Who records coverage.** A clock assembled from other clocks' scores records coverage **iff every
component contributes to every sample**. `GrimAgeV1` (8 surrogates) and `DNAmFitAge_{Sex}` (3
same-sex members + `GrimAgeV1`) qualify -- every sample they score consumes all of them, so the
declared panel describes each sample truthfully.

**A composite's panel is what the catalog declares, not the closure over its dependencies.**
`DNAmFitAge_{Sex}` declares 172 (F) / 190 (M) CpGs -- exactly the union of its three *fitness*
members, sharing a single CpG with `GrimAgeV1`'s 1030. So its coverage figures, on both axes,
describe the fitness panel only; the `GrimAgeV1` contribution to the KDM mix is accounted on
`GrimAgeV1`'s own column (and its 8 surrogates'), all of which are returned as dependency columns.
Both floors read the same declared panel, so the two agree with each other -- but "`DNAmFitAge_Male`
is at 100%" is a statement about the fitness CpGs, not about everything feeding the score. A
**sex-routed alias** does not qualify at all: its two members apply to
disjoint halves of the cohort, so no aggregate over their union is true of any sample. Its
`per_clock` entry is `NULL` -- the panel counts live on the members, which keep their `per_clock`
rows even though they are not columns. This is a scoring-contract rule (`score_sex_routed()` emits
no coverage record), not a catalog fact -- the alias also has no panel, so the sec 2.1 floor skips
it too.

**What crosses a routing split is decided per-sample vs per-panel, not score vs coverage.** A
routed sample was scored by exactly one member, so anything indexed *by sample* routes with it and
stays true of its row: the score, and the tier-2 `sample_miss`, both stitched in
`score_sex_routed()`. Anything indexed *by panel* does not, because the panels differ (172 F / 190
M) and a count only means something against the one it was counted over -- "4 imputed" is 4 of 172
for a woman and 4 of 190 for a man. Merging those numerators into one figure while the
denominators stay on two separate rows yields a number nobody can read, so `per_clock[[alias]]`
stays `NULL` and the caller reads the member rows for the denominators.

The same asymmetry is why `check_row_coverage()` runs over every clock computed rather than the
returned columns (sec 2.1): the row floor **is** a per-panel division, so it can only be evaluated
on the member records, which have no columns of their own.

### 4.2 Tier 2 — per-clock x per-sample missingness (QC matrix, `samples_coverage()`)

Tier 1 collapses the sample axis; tier 2 keeps it. Absence is sample-invariant, but **partial NA is
sample-specific**: a sample 10% missing *globally* can be ~100% missing *for one clock* if its gaps
concentrate in that clock's scoring set — invisible to the global empty-row check (§2.3a), and **not
an error** (one clock being locally empty for one sample must not fail the batch — record it, don't
throw). `compute_coverage()` builds a length-`n` vector = row-wise NA count over a panel's present
CpGs, **once per distinct panel** and **per panel role** (score always; norm when the clock
normalizes -- only DunedinPACE), fanned out to clocks via `resolve_cpgs()`'s panel index. Assembly
stacks these into `$coverage$sample_miss = list(score = <n x k>, norm = <n x k' over just the
normalizing columns>)` on the result record, for a sample x clock coverage heatmap per panel.

This is what the **row floor** consumes: `check_row_coverage()` (§2.1) reads the per-panel miss
alongside `$coverage`, derives each sample's observed fraction over the clock's declared row-gate
panel (`normalizes` -> norm, else score) as `(present - sample_miss) / needed`, and warns. That is
the whole reason no scoring branch takes `min_samples_coverage` as an argument -- coverage already
carries the per-sample signal, so the gate is one function over the finished records rather than a
check repeated in ten branches.

It is literally the row-wise decomposition of tier 1's per-panel `score_imputed_partial` /
`norm_imputed_partial` (sum the sample axis of that panel's vector and you get the scalar back), so
no new mechanism — same `slideimp::mat_miss` primitive, `col = FALSE`, coerced to integer.
Efficient by the same gate cascade: skip entirely when `!scan$has_na`; per clock,
`intersect(score_present, partial_na_cols)` is a **free** set op (reuses the prepare-side
`mat_miss(col = TRUE)` result — no re-scan), and `mat_miss(col = FALSE)` runs only over that clock's
NA-bearing columns (clean columns add 0 to every sample). Store only the partial vector; combine with
tier 1's absent scalar at render (`observed_fraction[i,c] = (score_needed − absent − sample_miss[i,c])
/ score_needed`). Optional future hook: a soft policy that flags/`NA`s a cell below a coverage
threshold reads this matrix — default stays record-and-report.

### 4.3 `report()` -- one QC entry point (DNAm input and/or scores)

`report()` is the paper's `qc(DNAm)` and `report(result)` collapsed into a **single verb that routes
on what it is handed**: a DNAm matrix (input QC), an `mc_result` (score QC), or both (both). A bare
`mc_result` in the first position is auto-routed to `result`, so `report(res)` works. At least one of
the two is required; neither is an error. It builds a structured `mc_report` (S3 over `list`), prints
headline problems to the console immediately (format, out-of-range betas, low-coverage clocks), and
writes a **PDF** (base `grDevices::pdf()` + base graphics -- no rmarkdown/LaTeX/pandoc, no new hard
deps, CRAN-safe and offline). Computation and rendering are split: `report_dnam.R` / `report_score()`
return data; `report_render.R` draws it. Tests assert the `mc_report`, never the PDF pixels.

- **DNAm input QC** (`report_dnam.R`). Reuses the scoring machinery rather than re-deriving it:
  `coerce_dnam()` (an all-numeric data.frame -- how some cleaned datasets ship DNAm -- is coerced to
  a matrix with a note; a non-numeric column is refused; shared with `calc_clocks`);
  `report_check_dnam()` (structural gate -- aborts on a non-matrix / missing dim names, records
  softer issues like a transposed matrix); `detect_array()` (nearest-size guess over cg-probe count
  vs `MC_ARRAY_SIZES`, with an EPICv2 replicate-suffix override -- a heuristic, never a hard identity,
  since array is **not** a catalog field); `check_beta_range()` (whole-matrix range, count of values
  `<0` / `>1`, NaN/Inf, offending probes via `matrixStats::colRanges`, a scale note for
  M-values/percentages); missingness split partial-NA vs fully-missing vs absent-from-matrix, straight
  off `scan_missing_cpgs()`; a **per-sample** block (`report_samples()`: per-sample NA fraction,
  mean/median beta, middle-band fraction, MAD-outlier flags, near-identical-pair detection via
  correlation, and compact per-sample beta density curves); and a per-clock coverage table from
  `clock_panels()` + `resolve_cpgs()` over `resolve_clocks(clocks)` (default `"all"`). Bimodality is
  judged from each sample's density **shape** (center-density / flanking-mode ratio; ~0.1 for real
  methylation, ~1 for flat/mis-scaled) -- **not** a middle-band fraction, which real data (lots of
  intermediate CpGs) drives to ~0.5 even when clearly bimodal. Sample-level flags are all **relative
  to the cohort** (a sample unlike its peers); a whole cohort that is not bimodal is a separate
  one-off signal (`cohort_bimodal_ok`), not a flag on every sample. External-pack
  clocks are assessed **only when their pack is already cached or `assets=` is passed** -- `report()`
  never triggers a surprise download; skipped groups are named in the report.
- **Score QC** (`report_score()`). Per-clock distribution summary (mean/sd/min/median/max, NA count),
  the samples that failed to score (NA/non-finite), a clock-clock correlation matrix (internal
  consistency), and `clocks_coverage(result)`. The score-distribution **reference is not in yet**:
  `mc_score_reference()` returns `NULL`; the comparison is wired to consume a `data.frame(clock_id,
  ref_mean, ref_sd, expect_lo, expect_hi)` -- `ref_mean`/`ref_sd` drive `|z|>2` mean flags,
  `expect_lo`/`expect_hi` drive out-of-(training-)range flags -- lit up by `mc_score_reference()` or
  the `score_reference=` arg, and otherwise a "no reference yet" note.
- **Age/sex association check** (`score_associations()`, needs `pheno` with `Age`). Recomputes each
  clock's age correlation in the user's data and compares it to a **shipped, meta-analytic
  expectation** (`inst/extdata/clock_reference.csv`, read by `mc_clock_reference()`): per clock, the
  pooled age correlation + a 95% prediction interval across ~136 public datasets (built by
  `data-raw/build_clock_reference.R`; see DECISIONS 2026-07-27). A clock outside its interval, or with
  the wrong sign, is surfaced -- but **advisory only**: age correlation shrinks with a narrow cohort
  age range, so it never enters the verdict, and the section states the cohort's age span as a caveat.
- **Verdict.** `report_verdict()` grades each section PASS/WARN/FAIL and takes the worst as the
  overall verdict (stored in `$meta$verdict`, shown in `print()` and atop the PDF).
- **Rendering** (`report_render.R`). Base `pdf()` text/table pages plus base-graphics plot pages:
  per-clock coverage histogram, overlaid per-sample beta densities (flagged samples in red), score
  histograms (small multiples), and a clock-clock correlation heatmap. Default output is
  `methylCIPHER-report.pdf` in the **working directory** (findable, not a temp dir); pass `file=`
  for elsewhere. The record is returned invisibly; the important flags (verdict + actionable issues)
  print to the console via `cli` on every run.

See DECISIONS 2026-07-27.

---

## 5. Phenotype and covariates

| Name | Encoding |
|---|---|
| `Age` | Numeric |
| `Female` | `1` female, `0` male |

- Only canonical names at scoring (aliases may be documented later, applied in prepare).
- Align by sample id only, never row order; full identity/alignment contract in §5.1.
- Covariate requirements are **one flattened catalog field** (`covariates_required`, computed by
  sync's `extract_covariates`). It unions every source that can imply a covariate: top-level
  `covariates` and any recipe step key ending in `covariates`. R reads the flattened field once at
  runtime; it does **not** re-derive from the meta keys in the scoring path. (The source set is
  empirical -- it grew as clocks were implemented -- so it is the accessor/fixture surface, not a
  fixed grammar; a wrong extraction surfaces as a covariate error or score mismatch in that
  clock's golden fixture.)
- Covariate *weights* have two homes and `clock_covariate_coefs()` reads both: top-level
  `covariates` on single-tensor clocks, and the recipe step producing `score` on recipe-borne ones
  (the DNAmFitAge members). Missing the second reads as "no covariates" and silently drops the
  term, so it is asserted directly in `test-score-fitage.R`, not just through a score golden.
- No clock branches on `Female` to pick its own weights any more. Sex is resolved in the
  `clock_id` (`{stem}_Female` / `{stem}_Male`); `Female` reaches the FitAge family only as
  GrimAgeV1's `stack` covariate.
- `calc_clocks` fails clearly when required covariates are missing.

### 5.1 Sample identity and pheno alignment

`rownames(DNAm)` is the canonical `sample_id`: unique and **mandatory**. Identity-less DNAm is a
hard error -- `check_DNAm()` refuses it and hands the caller the one-liner to name anonymous rows
themselves (`rownames(DNAm) <- paste0("sample", seq_len(nrow(DNAm)))`). The package never
manufactures ids: a nameless matrix is a degenerate input, and the honest response is an error, not
a silent guess whose meaninglessness then has to be tracked through the whole record as a
`positional`/`synthetic` flag (DECISIONS 2026-07-24, reversing the earlier `allow_positional_ids`
stamp-and-refuse-at-cbind design). Relocating the `sample1..N` step to a caller-written line moves
the "these rows are positional" admission to where the knowledge is. A caller with genuine sample
duplicates must give each its own id; a caller wanting a clock's score twice duplicates the score
column, not the input rows. `sample_id` is derived **once**, from `rownames(DNAm)`, stamped on
`$provenance$sample_id`, and used as `rownames($scores)`. pheno is never the identity.

`pheno_id` (the id column name) always has a default; only `pheno` itself may be `NULL`. Alignment
(`resolve_pheno()`) is **always an id join** -- there is no row-order mode, because there are always
real ids to match on:

| Case | Rule |
|---|---|
| `pheno = NULL` | `sample_id = rownames(DNAm)`; no covariate side-table |
| `pheno` given | `pheno[[pheno_id]]` present, **unique**, non-missing, and `rownames(DNAm)` a subset of `pheno[[pheno_id]]`; else error naming missing ids. Subset + reorder pheno to `rownames(DNAm)` order |

- **Id join.** `rownames(DNAm)` is a subset of `pheno[[pheno_id]]`, so we filter first; no
  `unique(pheno)` dedup -- duplicate ids are the caller's error, never silently collapsed. pheno may
  carry **extra** rows (other cohorts); they are ignored after the subset. This subset rule is what
  makes streaming ergonomic: pass the full-cohort pheno to every DNAm row-chunk and each chunk
  narrows to its own ids (`dev/id-streaming-plan.md`).
- Post-condition: the returned pheno has `nrow(DNAm)` rows in `sample_id` order, with
  `pheno[[pheno_id]] == sample_id` (reached for free since it matched on that column). Row names are
  **dropped** -- identity lives in the id column and
  nowhere else, so there is no second copy to drift, and data.frame row names do not survive a
  tibble/dplyr round-trip anyway. The frame is **narrowed** to `c(pheno_id, <covariate
  union>)` -- for the current catalog at most `ID`, `Age`, `Female`, since those are the only
  covariates any clock requires. Columns the caller supplied but no clock needed are dropped, so an
  arbitrary clinical table cannot ride into the record.
- That narrowed frame both feeds the covariate scorers and is stamped on the record as `$pheno`
  (sec 1.3), which is what lets `augment()` derive age acceleration and `cbind` gate 2 compare shared
  covariates per sample. `$scores` stays scores-only, and `as.data.frame()` surfaces `sample_id`
  under the `pheno_id` name beside the scores -- **covariates are not on that path**.

---

## 6. Batch-dependent clocks

| Clock | Op | Scope |
|---|---|---|
| `DNAmPhysAge` | `cohort_zscore` | Samples in current call |
| `DNAmPhysAge_years` | `cohort_zscore` | Samples in current call |
| `Zhang2019` | `sample_scale` | All probes within each sample (moments) |

Rules:

- Subsetting samples before vs after scoring can change batch-dependent outputs. Different
  cohorts are analyzed separately (dedupe samples); there is no "same sample across cohorts."
- Partial-NA imputation (sec 2.3) is also cohort-dependent by design -- imputation borrows
  within-cohort information -- but it is a fallback, not a scoring op, so it is not listed here.
- Do not global-union-probe-optimize away Zhang's full panel.
- **No stored `batch_set_id`.** It was `digest(sort(sample_id))` -- a fingerprint of the record's
  own ids, carrying nothing not already in `$provenance$sample_id`, and hashing the *current* ids
  meant it could not even survive a `[` subset to catch a cohort mismatch. The compatibility check
  reads the id set directly: `cbind` requires equal sets, a streaming `rbind` (future) requires
  disjoint sets. Removed 2026-07-24; the future direction (`dev/id-streaming-plan.md` Phase 3) moves
  every cross-sample op out of the scoring loop, so no cohort is baked into a score to fingerprint.

---

## 7. Provenance (result record)

`x$provenance` carries:

- `sample_id`
- `pheno_id` -- the id column name, so `as.data.frame()` and `augment()` agree on the join key
  without either taking an extra argument
- `clocks` — every scored column's catalog id, in `$scores` column order
- `requested` / `dependencies` — the partition of `clocks` into what the caller asked for and what
  the plan pulled in (§1.3, §2.1). The only place that distinction survives; both are real columns
- Covariates actually used — unioned over **all** returned columns, deps included, so a request for
  `DNAmFitAge` alone reports `Age` (its deps' requirement) and not just `Female`
- Coverage (see sec 4)

No commit SHA / pin is stamped on the result: correctness is proven by fixtures. The one sync-side
identity key that remains -- the external packs' content-address `payload_hash` -- stays in
`data-raw/` and never reaches a result record.

`augment()` compares covariates actually used; mismatch = error by default (warn/ignore
overrides optional).

### 7.1 `cbind` compatibility

`cbind.mc_result` ([`R/generics.R`](../R/generics.R)) binds score columns of two or more records
after **two** independent gates. (No positional-id guard is needed: rownames are mandatory, so
every `sample_id` is real identity -- DECISIONS 2026-07-24.)

1. **`sample_id` set.** Sets must be equal. Equal set + same order binds directly; equal set +
   different order reorders the later records to the first record's `sample_id` order and
   re-verifies `identical(sample_id)` before binding; unequal sets throw (different samples).
2. **Shared-covariate consistency** -- evaluated **per covariate**, only for covariates that appear
   in `covariates_used` of **more than one** record. Their per-sample values (keyed by `sample_id`)
   must agree, so the same `Age` / `Female` is bound to the same sample across records. Moot when
   only one side carries covariates, or both do but they are disjoint (nothing shared to compare); a
   covariate used by exactly one record is never checked. Reads the values from each record's
   `$pheno` (sec 5.1) keyed by `sample_id`, so a mismatch can name the offending samples. That frame
   is already narrowed to the id column plus required covariates, so this never touches the whole
   pheno table (sec 12 non-goal). An earlier draft stored a per-covariate hash instead; a hash cannot
   be recomputed after a row subset, which would silently void this gate for any record narrowed
   by `[` (DECISIONS 2026-07-23).

There is no third batch-cohort gate. Gate 1 already compares the full `sample_id` set, so a
`digest(sample_id)` gate would be redundant; the "same samples scored in two cohorts" case it
nominally guarded is dissolved by moving cross-sample ops to `augment` (§6, `dev/id-streaming-plan.md`
Phase 3). `batch_set_id` was removed 2026-07-24.

---

## 8. Source of truth and catalog discovery

Canonical remote: `https://github.com/hhp94/methylCIPHER-meta.git`.

Sync-time crawl (`data-raw/sync.R`):

```r
metas <- list.files("weights", pattern = "\\.meta\\.json$", recursive = TRUE, full.names = TRUE)
clock_metas <- metas[basename(metas) != "_group.meta.json"]
group_metas <- metas[basename(metas) == "_group.meta.json"]
```

- Discriminate clock vs group by basename, not depth.
- **The crawl locates metas; `manifest.json` decides which are vendored.** `manifest_clocks()` takes
  the clock set from the manifest's `clocks[]` (exactly the `weights_status=done` set) and errors on
  a manifest row with no meta on disk. Upstream explicitly permits staging a clock's meta before it
  is stamped done, and the crawl alone would vendor that unstamped clock.
- **No path is discovered by pattern.** `declared_tensors()` reads the named pointers
  (`imputation.ref`, `components/shared/probe_sets[].file`, `code_ref`, `code_deps`) and applies the
  `coef_path()` rule -- `weights/{group_id}/{clock_id}.csv.gz`, derived from the two keys -- for
  `cpg_coefficient`. Required there, so absence is an error, never a signal. Tensor shape comes from
  the declaring `row_key` / `col_key` (`col_key` optional: a rotation's PC columns are generated),
  asserted against the file's own header rather than inferred from column count.
- Singleton: `group_id == clock_id`; multi-member groups have `_group.meta.json`.
- At **sync time** R reads `weights/`, `manifest.json`, and `bibliography/{papers.csv,clocks.bib}`
  (the last joins `pmid -> bib_key` in memory and vendors `clocks.bib` to `inst/`; `papers.csv` is
  never shipped). Nothing under `control/`, `papers/`, or `scripts/` is read. At **runtime** R reads
  only the built `sysdata` objects, via the accessor layer -- never the meta tree.
- Family membership from `group_id` / sidecar path, not free text.

**Schema = accessors, not a doc.** The sysdata shape is defined by `build_index()` /
`build_catalog()`. The executable schema is the accessor layer -- `get_clock()`,
`clock_scoring_cpgs()`, `clock_norm_cpgs()`, `clock_impute()`, `clock_coefs()`,
`clock_group_bundle()` -- documented via short `#` comments for now (roxygen deferred to alpha)
and whose structural `testthat` test asserts `names(mc_index)` + a sample `mc_catalog` entry. A
shape change breaks the test and forces the comment update. `calc_clocks` code consumes only accessors, never raw nested lists. No
hand-written `schema.md` (it would drift, the class the meta repo's `meta_schema.py` exists to
eliminate); generate one from the objects if ever needed.

Meta remains the scientific contract; Python `ops.py` / `transforms.py` / `covariates.py` in
upstream are reference implementations for parity thinking.

---

## 9. Sync, lockfile, packaging

### 9.1 Sync-internal identity (not product provenance)

The full `sysdata.rda` rebuild is cheap (~2s), so there is **no** catalog build-skip cache and no
`manifest_key` (removed 2026-07-20). The identity keys that remain are scoped to the three external
qs2 packs (SystemsAge, PCClocks, PCBrainAge):

- `payload_hash` -- the **only** identity key. A content-address over a version-pinned serialization
  (`serialize(version = 2L, xdr = TRUE)` + `digest` sha256); sets the pack filename
  `<group>-<payload_hash>.qs2`, its GitHub release tag, and the reupload skip key. It is
  maintainer-side only: it is not a field on the shipped registry row (the filename already carries
  it) and R never recomputes it at runtime. Transfer integrity and bit rot are qs2's own
  `validate_checksum` on read (see 9.4).

A gitignored `data-raw/assets/lockfile.rds` (keyed on the external clocks' `bundle_hash` vector from
`manifest.json` + presence of the staged packs) skips only the external-pack rebuild; it is **not**
a product pin, not stamped on results, and not a catalog-change detector. No SHA is treated as
scientific identity. `bundle_hash` replaced `source_git_sha` as the key on 2026-07-24: the sha moved
on every upstream commit including docs-only ones, so the cache busted constantly and still could
not say *which* clock changed, whereas a `bundle_hash` moves iff that clock's meta or one of its
declared artifacts moved.

### 9.2 Sync workflow (`data-raw/sync.R`)

1. Resolve meta commit; checkout (`source_git_sha`).
2. **Always** build catalog + accessors' backing objects + small bundles -> `R/sysdata.rda` (no
   build-skip cache; the rebuild is ~2s).
3. External packs: if `force = FALSE` and the asset `lockfile.rds` hits (every external clock's
   `bundle_hash` unchanged and every staged pack still on disk), reuse them and restore their
   resolved probe sets from the lockfile; else rebuild the three content-addressed packs and
   rewrite the lockfile.
4. `upload = TRUE` publishes packs to GitHub Releases; the content-address (`payload_hash` ->
   filename -> tag) plus the remote "asset already present" skip make reupload idempotent.

Rules: no tensor rehash in R; sync does not score clocks; bundle inputs only whole `weights/`
paths (never `papers/` or `scripts/`); missing referenced weights path = error.

**Sex-routed aliases.** A group whose `_group.meta.json` carries `routing.sex` gets one alias clock
per stem, minted by `attach_sex_routed_aliases()` after probe-set resolution (so an alias never
acquires a panel). The alias is `kind = "sex_routed_alias"` -- a package-minted classification in
its **own** column, leaving `weights_format` / `computation_type` `NA`, because neither `"routed"`
nor `"sex_routed"` exists in upstream's `meta_schema.py` enums and a user filtering on
`weights_format` must never see a locally-invented member. `score_type()` and the coverage alias
split check `kind` *before* reading either enum. The alias
owns no weights, declares `covariates_required = "Female"`, and depends on its two members. The
members stay in the catalog but leave the **callable pool**: `resolve_clocks()` derives the
not-callable set from the same routing tables (`sex_routed_members()`), so the pool, the refusal
and its suggested alias can never disagree. Requesting one is a hard error naming the alias.
They also leave the **output**: `drop_routed_members()` (same source) keeps them out of
`output_ids`, so a routed member is scored, feeds its alias, keeps a coverage row -- and never
becomes a score column.

**Scoring-panel resolution.** `resolve_group_scoring_probe_sets()` is a memoized DAG walk, not a
fixpoint with fallbacks. The rule is total:

```
scoring_cpgs(clock) = union of the clock's own cpg-keyed tensors
                      (coef_path for cpg_coefficient; row_key == "cpg" components otherwise)
                    = if it owns none (score-assembled):
                      union(scoring_cpgs(c) for c in every recipe step's `inputs`)
```

Upstream guarantees `inputs` holds only `clock_id`s (the operand namespaces `inputs` / `internal` /
`covariates` are disjoint and validated), so this resolves in dependency order with a cycle guard.
Out-of-group inputs are excluded from the panel on purpose: `DNAmFitAge_{Sex}` reads `GrimAgeV1`,
but GrimAgeV1 is its own column with its own coverage row, and folding its 1030 probes in would
double-count and swamp the family's own 172/190. The DAG **edge** is still kept (`clock_inputs`) --
this is a coverage-reporting decision, not a dependency one.

`assert_declared_n_cpgs()` then checks the derived panel against upstream's declared `n_cpgs` for
every clock, hard error on mismatch, no exemption list. That is what proves the exclusion above is
right, and it is the cheapest guard in the sync: the retired three-tier version (own tensors ->
`depends_on_clocks` falling back to `covers` -> the group's shared bare CpG list) handed the
DNAmFitAge composites the family-wide 627-probe `data_prep2` prep panel against a declared 172.
`covers` is the **output** set (every clock the blob or recipe tree serves) and is never an input
list.

**Custom-group payloads.** A `weights_format = "custom"` group declares frozen author code as its
definition -- but "the code is the definition" was never a licence to leave the code's own inputs
undeclared. Upstream closed that: a `custom` meta declares `components` (frozen parameter tensors,
`row_key: cpg` like any other) and `code_deps` (the code closure `code_ref` reads or sources), and
an undeclared file under a `custom` group dir fails upstream validation. So MiAge's site-specific
`(b, c, d)` are one declared `weights/MiAge/site_parameters.csv.gz` (`cpg,b,c,d`) that loads like
any other cpg-keyed tensor, and its panel resolves by the ordinary rule above -- **no special case
at all**. The old `CUSTOM_GROUPS` registry (a loader plus a grafted component declaration, needed
when the blob was an undeclared `.Rdata` aligned to the panel by position) is deleted, along with
the positional-alignment guard that was the only thing catching a silent mis-key. `code_deps` files
ride into the pack as `type = "r_source"` alongside `code_ref`; before this, `function_library.r`
was bundled by nothing and the vendored scorer could not have run.

### 9.3 Distribution tiers

| Tier | Contents | Delivery |
|---|---|---|
| Bundled | Small groups (~0.8 MB class) | `R/sysdata.rda` |
| External | SystemsAge, PCClocks, PCBrainAge | Release assets |

- `data-raw/` not in tarball; no large tensors in `inst/`.
- No download at install / build / check / CRAN test.
- Double precision coefs only (no float32).

### 9.4 External asset resolution

`load_mc_assets(groups, assets = NULL, ask = TRUE)` in [`R/mc_data.R`](../R/mc_data.R)
is the single runtime entry (deliberately small -- see the 2026-07-21 DECISIONS entries). It returns
a **named list of packs keyed by `group_id`** (even for one group). `calc_clocks()` calls the
identical function internally, so a pre-loaded object and an auto-loaded one cannot drift. Flow:

1. **Cache dir precedence:** an `assets` **path** > session option `mc.cache_dir`
   (`options(mc.cache_dir = )`) > `MC_CACHE_DIR` (.Renviron) > `tools::R_user_dir(.., "cache")`
   (`mc_default_cache_dir()`). `mc_cache_dir()` takes a path or `NULL` **only** -- a loaded pack
   names no directory, so it is an error there rather than a silent fall-back to the default cache
   (the cache verbs below pass their raw `assets` in; see DECISIONS 2026-07-23).
2. **Open vs closed set (from `assets`).** `assets = NULL` -> **open**: resolve each group from the
   cache dir; missing packs are consent-downloaded. `assets` **explicitly provided** -> **closed**:
   resolve only from what is given, **never download**; a needed group not covered is a hard error.
   `assets` accepts a cache-dir path **or** loaded object(s) -- a bare pack (a list with `$group_id`),
   a list of packs, or a path all canonicalize to the named-list registry (objects key by their own
   `$group_id`); an asset the plan does not need warns and is ignored.
3. **Expected file** = `mc_provenance$external_assets[[group_id]]$file` (`<group>-<payload_hash>.qs2`).
   Present -> read it. Open-set download is consent-gated: `ask = TRUE` prompts interactively (one
   **batched** prompt for the union of missing packs) and **refuses** non-interactively; `ask = FALSE`
   is the explicit-consent signal. `ask` is a **strict flag** -- anything that is not a single
   non-NA `TRUE`/`FALSE` is an error, never silent consent. Staged to `<file>.part`, validated by a
   qs2 read, then atomically renamed (`fs::file_move()`, which errors rather than returning FALSE)
   -- a truncated transfer never lands as cached. `download.file()`'s status code is checked as well
   as its conditions, since it reports some failures only through the return value.
4. **Read** with `qs2::qs_read(validate_checksum = TRUE)`; qs2's own checksum is the sole integrity
   guard. There is no second hash: no sha256 recompute, and no re-hash of the loaded payload (the
   content-addressed filename already asserts which pack it is). A just-fetched pack is **not**
   re-read: `mc_fetch()` hands back the payload from its validating read.

**No memoise.** A cold `qs_read` is ~0.1-0.2s for all three packs (benchmarked 2026-07-21), so
`calc_clocks()` resolves the needed groups **once per call** (in the prepare phase, before the pure
scoring loop) and threads the registry down; there is no session-cache/global-env tier. No silent
first-use download.

Pack encoding (`canonical_matrices`) is per-group: PCClocks/PCBrainAge carry one
`coefficient_matrix` + `impute` vector; SystemsAge carries `organs`/`systems`/`age` matrices +
`impute`. Within a group the probe set is shared once as matrix columns.

`groups = "all"` (the `mc_data_download()` / `clear_mc_cache()` default) resolves to every external
group; an **empty** request stays empty -- `character(0)` from `pack_groups_needed()` means "no
external pack needed", never "all of them".

Verbs: `mc_data_download(groups, assets, ask)` pre-fetches; `clear_mc_cache(groups, assets, ask)`
removes cached packs under the same consent gate as download -- `ask = TRUE` prompts interactively
and **refuses** non-interactively, `ask = FALSE` is the explicit-consent signal, and declining is a
quiet no-op (not an error). It only ever deletes registry-named files inside the resolved cache dir,
and stops if a file survives the delete. Both prompts print an aligned group + size table (sizes via
`fs::fs_bytes()`, totalled when more than one pack); the delete prompt labels rows by **group id**,
not the content-addressed file name, which is why `mc_cached_files()` names its result by group.

### 9.5 Verification status

From manifest: `"" | pending | verified | disputed | skipped`. Separate from format / fixture
parity. It is a **maintainer/dev signal only** -- it records why a clock is (e.g.) `skipped`
(Horvath online variants we deliberately omit) for the human running sync. It is **not** part of
the user surface: `list_clocks` does not expose it and it is not canonicalized into the shipped
catalog. Sync still validates the value against the allowed set at build time (a warning on an
unexpected value), but the field itself stays on the maintainer side (manifest), not in
`mc_index` / `mc_catalog`. R does not read upstream ledger files.

---

## 10. Fixtures and CRAN

### 10.1 Two test tiers, different jobs

- **Engine / machinery unit tests -- always run, no meta dependency.** Hand-authored toy inputs
  (a few probes x a few samples, coefficients chosen so the expected dot product is written by
  hand in the test). Cover `linear_score` arithmetic, the impute accounting (partial vs full),
  accessors, coverage math, and result-class methods. Golden values tie to nothing upstream ->
  zero drift.
- **`sim_DNAm` smoke -- always run, no meta/cohort dependency, no golden.** `sim_DNAm(n, clocks)`
  ([`R/sim_DNAm.R`](../R/sim_DNAm.R)) reads scoring probe sets via `clock_scoring_cpgs()` -- the
  shipped `mc_catalog` for bundled clocks, the pack for external groups, so it takes `assets` /
  `ask` and requires the pack for an external request -- and throws
  random uniforms at them; scoring is asserted with `expect_no_error` (the data need not be
  scientifically meaningful). This exercises dispatch, coef loading, and assembly over every
  bundled clock -- and it fails loudly if a clock resolved to an empty scoring set -- without any
  `expect_equal` and so without any drift risk. It is a *machinery* check, not a correctness one.
  It resolves its panel through **`resolve_clocks_sequence(resolve_clocks(clocks))`**, the same plan
  `calc_clocks()` computes — mirroring namespace resolution alone would emit a panel missing every
  cross-clock dependency's CpGs (`"DNAmFitAge"` -> the group's 627, where the plan needs 1643), and
  since the GrimAge surrogates are policy `omit` they would score intercept-only and the composite
  would come back a plausible-looking number rather than an error. `remove = k` drops `k` CpGs to
  exercise the missingness path, and is the intended way to trip the §2.1 coverage floor in tests.
- **Parity fixtures -- the single clock-golden source, cohort-gated.** Upstream golden fixtures vs
  **every registry cohort**: `fixtures/{cohort}/beta.duckdb` for `cohort_EPICv1` and `cohort_450K`
  (gitignored; regenerable via `fixtures/build_cohort.R`; paths under the meta clone
  `data-raw/methylCIPHER-meta/`). Upstream ships one `fixtures[]` block per cohort, so each
  (clock, cohort) pair is its own test, with one duckdb connection per staged cohort. Skipped unless
  `MC_PARITY=1` is set and that cohort is staged (`file.exists()`); run locally via the dev-only,
  build-ignored `test_parity()` (`R/dev-utils.R`).

  **Tolerance is ours and is keyed on a declared field, not a clock list.** Upstream retired its
  `fixtures[].parity` policy enum. A fixture whose `server_normalization` is non-empty was produced
  by the Horvath calculator on betas the **server** normalized (BMIQ / noob); we score raw betas and
  deliberately vendor no BMIQ gold standard, so those are graded by **correlation** (> 0.99). Every
  other fixture's oracle saw the same raw betas we do and is held to **exact** (max_abs_diff <
  1e-6). That splits the 244 fixture blocks 186 exact / 58 correlation, and it is the real numeric
  gate. A hand-kept list of "clocks that cannot be exact" would rot the moment a clock's oracle
  changed. `KNOWN_PARITY_GAPS` holds only genuine skips (today: `Zhang2019`).

No slice of the golden cohort is committed into `tests/` -- that is a second copy of upstream
golden values and it drifts; a small slice also cannot validate a `correlation`-policy clock.
CI may stage/regenerate the cohort and run the full parity gate; CRAN skips it. DuckDB is for
fixture lookups only, not a runtime tensor store.

### 10.2 CRAN checklist

- Drop or clearly deprecate legacy exports before submission as required.
- Valid License; no `library`/`require` in scoring; Suggests soft-fail for optionals.
- Teaching deps must not block install; tests without external assets; no network at check.

---

## 11. Implementation sequence

1. Accessor layer over `sysdata` (executable schema) + its structural test; catalog readers.
2. `prepare_inputs`, result record, coverage, `clocks_coverage()` / `samples_coverage()`.
3. Linear engine + impute policies (partial/full split); fixture-drive `cpg_coefficient` batch.
4. `sample_scale` / Zhang transform -> linear; localized full-panel warn.
5. GrimAge orchestrator (alias + components + V1/V2 columns).
6. SystemsAge (and PC\* as needed) + asset resolver.
7. Remaining component_matrices as more branches/packs; external + MiAge; DNAmFitAge sex-resolved
   members + routed aliases (done).
8. Full fixture suite; optional legacy wrappers; docs / CRAN cleanup. **(current)**

Track clocks in `dev/clock_tracker.csv` (`uv run python dev/build_clock_tracker.py`).

---

## 12. Non-goals

- Expanding per-clock exported calculators as the real engine.
- A recipe-walker runtime (rejected outright, not kept as fallback).
- Returning full phenotype tables from scorers.
- Hashing entire phenotype tables.
- Silent network downloads.
- Deep OOP for weighted sums.
- Rewriting published clock mathematics.
- Overriding metadata `drop` / impute policies package-wide.
- Crossing the impute sources (partial=cohort, absent=vendor).
- float32 coefficients.
- Bioconductor.
- SHA / pin as result provenance; R-side tensor rehashing.
- A `matrix`-subclass result; a hand-written `schema.md`.
- Micro-optimizing matrix copies before evidence.
- `$` on catalog / pack / tensor structures (partial matching returns a wrong value silently).
- Locating a payload by `grep`/regex instead of resolving its declared pointer.
- A coverage figure aggregated over components that no single sample consumes.
- Exposing every catalog `clock_id` as callable when some are routing targets.

---

## 13. Follow-ups

- A discovery helper listing the callable pool, with Levenshtein "did you mean" on a typo. Must
  read `sex_routed_members()` so a rejected routing target answers with its alias rather than a
  distance guess (sec 8).
- Whether `mc_groups[[g]]$members` should mean "callable" (aliases) or "real clocks" (today's 14).
  `test-fixtures-parity.R` reads it expecting scoreable ids.
- Make the exact-tolerance parity grade **scale-aware** (numpy-allclose:
  `|got - exp| <= PARITY_ATOL + PARITY_RTOL*|exp|`, atol `1e-8`, rtol `1e-6`), so a `1e-8` diff on a
  `2.7e6`-scale PC biomarker is not judged like `1e-8` on an age. **Not implemented** -- the tier
  currently grades `max_abs_diff < 1e-6` (or correlation, sec 10.1). It matters most for the PC
  clocks, which are external-pack-only and therefore skipped until the packs are cached, so it is
  not currently blocking anything. The old motivating case is gone: the DNAmGrip / GrimAge-surrogate
  skip class was dissolved on 2026-07-24 when grading moved to the `server_normalization` split, and
  `KNOWN_PARITY_GAPS` now holds only `Zhang2019`.
- Whether SystemsAge members can be derived from shared components only.
- Whether fixture cohort should be published for cross-language consumers.
- Permanent vs one-release legacy `calc*` wrappers.
- Exact GrimAge pack column set and alias errors (`GrimAgeV1` typed by user -> ?).
- **Finalize license + author/copyright credits before the alpha / public release.** Copyright
  holder is currently Yale University (DESCRIPTION `cph`, LICENSE / LICENSE.md); the BSD-3 choice
  vs bundled third-party weight terms is still open (see DECISIONS 2026-07-17 / 2026-07-21).

---

## 14. Pipeline diagram

```text
methylCIPHER-meta
  weights/ + manifest.json + fixtures/
           |
    data-raw/sync.R  (lockfile.rds = external-pack rebuild skip)
           |
    R/sysdata.rda  +  external family assets
           |
    accessor layer (executable schema)
           |
    calc_clocks()
      prepare | linear engine | pack orchestrators | external/custom
           |
    mc_result record  --clocks_coverage()/samples_coverage()--> coverage data.frame
           |
        augment() --> analysis join
```
