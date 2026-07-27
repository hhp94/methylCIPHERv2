# methylCIPHERv2 Rewrite Plan (overview)

Compressed product plan for the catalog-driven CRAN rewrite. **Detailed design, contracts,
and worked examples live in [`detail-plan.md`](detail-plan.md).** Design *why* / reversals:
[`DECISIONS.md`](DECISIONS.md). Per-clock agent tracking: `clock_tracker.csv` (regen via
`build_clock_tracker.py`).

These plans state **current truth only**. Superseded designs are not kept inline; their history
(the "we tried X, reversed to Y") lives solely in [`DECISIONS.md`](DECISIONS.md).

---

## Objective

> Pick clocks -> pure scores -> optionally `augment()` with phenotype.

- One public scorer: `calc_clocks()`.
- Catalog + tensors from [methylCIPHER-meta](https://github.com/hhp94/methylCIPHER-meta.git);
  fixtures are the scientific gate.
- CRAN, not Bioconductor.
- Legacy per-clock `calc*` are not the implementation model. Optional thin wrappers that call
  `calc_clocks` may exist for one deprecation window (or longer if collaborators insist) --
  never reimplemented old bodies.

Bocklandt / Garagnani stay out (no published coefficients).

---

## Source of truth (one screen)

| Item | Rule |
|---|---|
| Remote | `https://github.com/hhp94/methylCIPHER-meta.git` |
| Inputs R may read | `manifest.json`, `weights/**`, `bibliography/{clock_citations.csv,clocks.bib}` (never `control/`, `papers/`, `scripts/`, `bibliography/papers.csv`) |
| Clock meta | `weights/{group_id}/{clock_id}.meta.json` |
| Group sidecar | `weights/{group_id}/_group.meta.json` (multi-member only) |
| Discovery | Recursive `*.meta.json`; clock vs group by basename |

The package carries no commit pin as product identity: correctness is proven by fixtures, not a
SHA. The only sync-internal identity is the external packs' content-address `payload_hash`; a
gitignored `lockfile.rds` skips just their rebuild -> [`detail-plan.md`](detail-plan.md) sec 9.

---

## Runtime architecture (current intent)

Runtime is **one linear engine + a small, closed set of branches** -- not a clock-function
factory and not a recipe interpreter. Each work unit is dispatched on the catalog pair
`(weights_format, computation_type)`. Meta/fixtures remain the scientific contract; the branches
are optimized implementations of that contract.

```text
calc_clocks(DNAm, clocks, pheno = NULL, ...)
  prepare once          # DNAm check, pheno align-by-ID, covariate requirements
  resolve requests      # callable pool only; group ids expand; routing targets refused
  expand deps           # transitive clock_inputs (recipe `inputs`), deps first
  route each unit on (weights_format, computation_type):
    linear engine       # most cpg_coefficient (+ optional pre-transforms)
    family packs        # GrimAge, SystemsAge, Dunedin (shared intermediates; not generic)
    sex-routed alias    # picks a sex-resolved member per sample
    external / custom   # MiAge
  mask routing targets  # blank rows a one-sex model does not apply to
  assemble              # mc_result record: scores + coverage + provenance
clocks_coverage(result) / samples_coverage(result)  # data.frames from coverage
augment(result, data)   # join scores to analysis tables
```

| Bucket | Role |
|---|---|
| **Linear engine** | Impute policy -> subset probes -> sum(w*x) + intercept -> optional transform |
| **Transform modules** | e.g. Zhang `sample_scale`: row moments on **full** matrix, apply to coef subset only |
| **Family orchestrators** | GrimAge / SystemsAge / DNAmFitAge composite: pack UX, shared work, multi-column return |
| **Sex-routed aliases** | The 7 DNAmFitAge stems: no weights of their own; select a member per sample |
| **One-offs** | `custom` (MiAge) -- `external_package` was retired upstream 2026-07-24 and left zero metas |

Packs may own their orchestration (shared intermediates, column assembly) but call the shared
linear/impute helpers for every linear sub-step -- imputation lives in exactly one place.

**Memory:** pass the raw beta into scorers; no global intersect/copy in `calc_clocks`. Each
scorer subsets to the probes it needs (small copy). Do not micro-optimize further until
profiling says so.

**Full-panel notice:** hard-coded to `Zhang2019` (the only `sample_scale` clock today) -- a
`message()` (not a warning) that its moments are computed over all CpGs but a big-enough subset
usually suffices. No generic flagged-clock warning; not on every call.

Details, imputation table, coverage/`clocks_coverage()`+`samples_coverage()`, GrimAge pack policy,
batch rules ->
[`detail-plan.md`](detail-plan.md).

---

## Public API (summary)

| Function | Purpose |
|---|---|
| `list_clocks()` / `get_clock()` / `get_clock_probes()` | Discover and inspect |
| `calc_clocks(DNAm, clocks, pheno = NULL, ...)` | Score -> `mc_result` record |
| `clocks_coverage(x)` | Per-clock coverage table (per-role needed / used / imputed / missing) |
| `samples_coverage(x)` | Per-(sample, clock, panel) coverage table (long; carries per-row denominators) |
| `augment()` | Join scores to phenotype / analysis data |
| `cite_clocks(x)` | Citations for ids/groups or an `mc_result` -> `mc_citation` (print / `as.data.frame` / `toBibtex`) |
| `get_mc_assets_dir()` / `set_mc_assets_dir()` / `list_mc_assets()` / `download_mc_assets()` / `load_mc_assets()` / `clear_mc_assets()` | Heavy external assets |
| Optional legacy `calc*()` | Thin -> `calc_clocks`; not the engine |

`calc_clocks()` returns an S3 record over `list` (class `"mc_result"`): `$scores` (n x k
double), `$pheno`, `$coverage`, `$provenance`. Verbs are methods (`print`, `as.matrix`,
`as.data.frame`, `[`, `cbind`, `augment`) -- so subsetting never silently drops
coverage/provenance; `rbind` refuses, and `clocks_coverage()` / `samples_coverage()` / `codebook()`
are plain functions, while citations dispatch on the package's own `cite_clocks()` generic
(detail-plan sec 1.3, sec 8.1). `$scores` and `as.data.frame()` are **scores only** (no
auto-appended pheno); `$pheno` separately retains the aligned id column plus required
covariates, which is what `augment()` reads. Align
pheno by sample id, never row order. `rownames(DNAm)` is the canonical sample id and is
**mandatory** -- rowname-less DNAm is a hard error that hands the caller the one-liner to name
anonymous rows themselves; the package never manufactures ids (DECISIONS 2026-07-24). Canonical
covariates: `Age`, `Female` (0/1).

The sysdata schema is the **accessor layer** (`get_clock`, `clock_scoring_cpgs`, ...) plus a
structural test -- not a hand-written schema doc. Accessors read the catalog with `[[` (never
`$`, which partial-matches) and resolve declared pointers rather than searching for them.
Covariate requirements are one flattened catalog field, read once, never re-derived in the
scoring path.

**Callable pool != catalog, and neither is the output.** Some clock_ids exist only as routing
targets: the 14 sex-resolved DNAmFitAge members are scored but never returned as columns, and
requesting one by name is a hard error naming its alias (`DNAmGrip_wAge_Female` -> use
`DNAmGrip_wAge`), because a one-sex model returns a plausible number for the other sex rather than
failing. A sex-routed family is one column per alias, filled for every sample. `"all"` and group
ids expand to callables only. A planned discovery helper (available options + Levenshtein "did you
mean") reads the same routing tables, so the pool, the refusal, the suggestion and the output
filter cannot drift apart.

**Coverage never describes a sample it is not true of.** A clock assembled from other clocks'
scores reports coverage only when every component contributes to every sample; a sex-routed alias
therefore reports no panel counts (its members cover disjoint halves of the cohort, with
different-sized panels) and those stay on the member rows, which survive without columns. **What
crosses the routing split is per-sample, not per-panel**: the score and the per-sample
`sample_miss` route to the member that scored each row, while a count that is only readable
against its panel stays on that panel's row -> [`detail-plan.md`](detail-plan.md) sec 4.

---

## Packaging (summary)

| Tier | What | Where |
|---|---|---|
| Bundled | Catalog + small tensors | `R/sysdata.rda` |
| External | SystemsAge, PCClocks, PCBrainAge | Release assets + user cache |

No install/check-time network. Sync via `data-raw/sync.R` (a gitignored `lockfile.rds` skips only
external-pack rebuild; not a product pin).

Where the upstream contract does not match what a caller needs, sync adapts it in a **small closed
registry** rather than adding a runtime code path. One remains: `attach_sex_routed_aliases()` (one
alias per `routing.sex` stem), emitting ordinary catalog entries so nothing downstream is
special-cased. `CUSTOM_GROUPS` (MiAge's frozen payload) is gone -- upstream now declares those
tensors, which is the better fix. See [`sync-boundary-migration.md`](sync-boundary-migration.md).

---

## Testing (summary)

| Tier | Runs where | Job |
|---|---|---|
| Engine units + `sim_DNAm` smoke | Always (no meta dependency) | `linear_score` arithmetic, impute accounting, accessors, coverage math, result methods (golden values hand-authored in-test); plus `sim_DNAm` `expect_no_error` over every shipped clock |
| Parity fixtures | `MC_PARITY=1` + cohort staged (`file.exists` gate); dev `test_parity()` | Upstream golden fixtures vs **every registry cohort** (`cohort_EPICv1`, `cohort_450K`) -- the single clock-golden source. Tolerance is ours and is **harsh on both axes**: `max_abs_diff` AND `max_rel_diff`, both reduced with `max`. `PARITY_REL_TOL` is `1e-10` everywhere with no exception; only the scale-sensitive `PARITY_ABS_TOL` varies by block (`core`/`fitage` `1e-10`, `packs` `1e-6`, measured). Correlation is banned as a gate package-wide, as are `median`/`mean` reducers. Four blocks from `parity_block()`, each selected by a declared field: `core`, `fitage`, `packs`, `horvath`. **214 pass / 32 skip / 0 fail.** The `horvath` block (30) is skipped on evidence -- the online calculator filled absent probes with an unpublished constant, and pairs with zero absent probes already match to ~1e-8, so the gap is the oracle's input, not our arithmetic |

No shipped slice of the golden cohort (it would drift). CI may stage the cohort and run parity;
CRAN skips it. Details -> [`detail-plan.md`](detail-plan.md) sec 10.

---

## Implementation sequence (high level)

1. Accessor layer over `sysdata` (executable schema) + catalog readers.
2. Prepare path + result record + `clocks_coverage()` / `samples_coverage()`.
3. Linear engine + impute policies + fixture batch for `cpg_coefficient`.
4. Zhang-style pre-transform into linear engine.
5. GrimAge / SystemsAge orchestrators; then other packs as needed.
6. External + MiAge; asset resolver. **(done)**
7. Full fixture suite; optional legacy wrappers; CRAN cleanup. **(current)**

Per-clock status: `dev/clock_tracker.csv`.

---

## Non-goals (short)

No new per-clock exported calculators as the real API; no recipe-walker runtime; no pheno tables
from scorers; no silent downloads; no float32 coefs; no Bioconductor; no SHA/pin as result
provenance; no rewriting published clock math; no crossing the impute sources (partial=cohort,
absent=vendor); no `$` on catalog structures; no searching for a payload an accessor should have
had declared; no coverage figure that is true of no sample.

---

## Doc map

| File | Contents |
|---|---|
| **This file** | Overview and pointers |
| [`detail-plan.md`](detail-plan.md) | Full design: API, engines, memory, packs, sync, fixtures, sequence |
| [`DECISIONS.md`](DECISIONS.md) | Dated why/reversals -- the only home for superseded design |
| `clock_tracker.csv` | Temporary per-clock agent board |
| `legacy/` | Frozen pre-rewrite `R/` sources |
