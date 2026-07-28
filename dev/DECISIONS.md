# Decisions log (package rewrite)

Append-only, date-stamped. Records **why** we chose a design during the rewrite under
`dev/migration-plan.md` — the “we tried X / other maintainers will ask this” history that
should not bloat the plan’s operational sections.

**Scope:** R package rewrite, packaging, API, and local maintainer workflow. Upstream
metadata contract decisions live in `data-raw/methylCIPHER-meta/control/DECISIONS.md`.

Newest first. Add an entry when a decision reverses a prior approach or is likely to be
second-guessed; do not restate rules already stated in the migration / detail plans.

---

## 2026-07-28 -- `from` -> `ext_data`: an argument name has to be greppable

**Decision.** The external-data argument on `calc_clocks()` / `sim_DNAm()` / `load_mc_assets()`
is now **`ext_data`**. Internals follow: `mc_spec()` and `mc_resolve_assets_dir()` take `ext_data`,
and `mc_canonicalize_from()` is `mc_canonicalize_ext_data()`. Semantics are untouched -- `NULL` is
still the open set, a path or loaded pack(s) still the closed set, precedence still
`ext_data` > `mc.assets_dir` > `MC_ASSETS_DIR` > `R_user_dir`.

**Why this reverses 2026-07-24.** That entry split the polymorphic `assets` argument into
"`assets` = the things, `from` = the source", and `from` was the right *semantic* choice: it names
where the packs come from rather than what they are. The problem is operational, not semantic.
`from` is a preposition that appears 383 times in this repo and dozens of times in `R/` alone, in
comments, cli strings and ordinary prose. There is no search that separates the argument from the
English word, so the one question you actually ask about an argument -- where is it used -- has no
answer. A name you cannot grep is a name you cannot refactor safely; that outweighs reading well at
the call site.

**Why not `assets_from`,** which would have kept both halves of the vocabulary and needed no
reversal: rejected as a compromise that keeps the preposition's reading without committing to
either noun. `ext_data` says plainly what the argument supplies and collides with nothing.

**Cost, paid in full here.** 4 formals, 3 pass-throughs, 1 renamed helper, 8 `{.arg from}` cli
strings (an unrenamed one names an argument that no longer exists), 15 named call sites across 3
test files, and the plan/CLAUDE.md text. Nothing serializes the argument -- `mc_spec()` stores the
resolved `packs`, not the source -- so `$provenance` and saved `mc_result` objects are unaffected
and there is no on-disk compatibility break. Pre-alpha, so user code calling `from =` breaks with
no deprecation shim, which is intended.

## 2026-07-28 (latest) -- cli is front-door only; internals use plain stop()

**Reverses** the earlier "every user-facing error goes through cli" rule. cli stays only on
the top interactive surface; catalog/score-branch/package-bug paths use `stop()` /
`warning()` / `message()` with `call. = FALSE`.

**Keep:** assets consent/download/clear (`mc_data.R`); discovery (`list_tags`,
`print.mc_citation`, `list_clocks` unknown group); public refusals (`rbind.mc_result`,
`cite_clocks.default`); calc_clocks front door -- `resolve_clocks` token errors +
did-you-mean, DNAm/pheno structure, coverage gates, value gates / all-NA samples,
`mc_cohort` missing pheno, `sim_DNAm` unresolved panels.

**Strip:** accessors, score branches, pack dispatch/finalize, `unroutable()`, normalize
arg validation, dependency-cycle, citation internals, `check_mc_result`, soft-dep
(`require_betanorm` plain stop with install hint), BMIQ sample-failure `warning()`,
`note_full_panel_clocks` -> `message()`.

Why: cli pluralization/`{.val}` machinery is high cost for package-bug and mid-score
paths the user cannot act on beyond "report it" or "fix the matrix"; the keep set is
exactly where multi-line guidance and consent UX earn their keep. `cli` remains in
Imports.

## 2026-07-28 -- comment hygiene: short `#` notes only; kernel *why* lives here

Package sources (`R/`, `src/`, `tests/`, `data-raw/*.R`) keep **short lowercase `#`
comments** that say *what* only when the code does not. Rationale essays, measured
tradeoffs, and "do not reverse this" notes do **not** live in-source -- they belong
here or in the plans. C++ keeps only `// [[Rcpp::export]]` (an Rcpp attribute, not
prose); every other `//` block is gone.

**C++ kernels (moved out of `src/`).** `col_stats()` is one column-major sweep over
the panel slice. Three cheap-inner-loop choices, measured on 100 x 39,025 at -O2:
raw column pointer (not `obj(i, j)`), `std::isnan` (not out-of-line `ISNAN` /
`R_finite`), and **one** predicate in the hot loop -- Inf is allowed into `sum`,
which poisons to non-finite, and one `isfinite(sum)` per column catches it (nr * nc
tests become nc). The cold rescan of that column separates Inf (report 1-based row)
from overflow (row NA; `check_col_values()` words it differently). Classification
the R side depends on still holds: NA/NaN skipped, Inf poisons, else summed.
`fill_imp_col()` is a branchless select (`isnan ? mean : v`) so the loop can
vectorize; it fills **missing** only (`std::isnan`, not `!isfinite`), so an Inf
that somehow arrived would stay visible in the score rather than be swapped for a
cohort mean. Measured 5.5 -> 1.2 ms vs the sugar-indexed original on that slice.

**BMIQ / notes channel.** Fully absent probes are dropped under BMIQ (not filled
from the gold target), because target-drawn values pull each sample's mixture fit
toward the gold and shrink the correction. A sample whose mixture will not fit is
NA in the score while coverage (computed before scoring) still reads full -- so the
failure is recorded in `notes`, not invented into coverage.

No behaviour change. Verify with `devtools::test()`; parity and `R CMD check` were
not run.

---

## 2026-07-27 -- a stack's column labels are declared upstream, so both copies of the `^raw_` strip are gone

Closes housekeeping item 11, which was parked as "needs an upstream declaration, not a downstream
patch". Upstream landed it (`methylCIPHER-meta` `ca6833c`, `STACK_COLUMNS_FIELD` in `scripts/ops.py`,
`weights_extraction.md` sec 4): a `stack` step's output columns carry **labels**, and the tensors a
later step names (`center`, `scale`, `rotation`, `coef`) are keyed by those labels **by name**. The
default label is the operand's own name; an optional `columns` list, parallel to the
`inputs`-then-`internal`-then-`covariates` concatenation, overrides it elementwise. SystemsAge is
the only stack in the repo that declares one -- GrimAgeV1/V2 and PhysAge hold under the default.

**The strip was never a second valid implementation, only an undeclared one.** `member_coef_files()`
radix-sorts, so the packed SystemsAge matrix is ordered `Blood, Brain, Heart, Inflammation, ...`
while the recipe declares `Blood, Brain, Inflammation, Heart, ...`: Heart/Inflammation and
Lung/Metabolic were **already transposed** between the pack and the recipe, and name-indexing was
load-bearing. Measured: permuting two entries of the declared `columns` in memory moves SystemsAge
by **93.19 years**, which is the local proof the labels reach the tensors.

**The guard that existed covered 11 of 12 and missed the only one that mattered.** It compared
`names(member_coef_files())` with the stripped recipe outs -- the same strip on both sides, so it
could only catch a spelling drift, and `Age_prediction` (a real member `clock_id` with no
`coef_path`, hence absent from `member_coef_files()`) entered no comparison at all. That was
precisely the label the hardcoded `if (operand == "ap_scaled") "Age_prediction"` produced. The new
guard compares the 12 declared labels against the group's 12 **member clock ids** (via
`group_member_ids()`), with no strip on either side, and keeps the narrower 11-way coefficient
check underneath it.

**Deleted, with no fallback path:** `SYSTEMSAGE_RAW_PREFIX`, `systemsage_label()`, the `ap_scaled`
literal, and `sub("^raw_", ...)` in `data-raw/sync.R`. A retained strip would be a second copy of
the convention that the upstream gate cannot see, which is the whole defect; and `sync()`
regenerates the catalog wholesale, so there is no old-catalog/new-code window to bridge. `"ap_scaled"`
and `"sysscores"` survive as recipe `out` names -- those are declared identifiers, not spelling
arithmetic. What replaces it is `stack_label_map()` / `stack_labels()` in `R/accessors.R`, alongside
the existing operand-space `stack_roles()` / `stack_operands()`. **Both spaces stay**, because they
are genuinely different: PhysAge indexes its `linear_mean` steps by recipe `out` (an operand) while
SystemsAge indexes the PCA tensors by label. Collapsing them would have been correct only while the
two coincide.

**Acceptance, all measured.** `payload_hash` **unchanged** (`2c180afc...`) across a maintainer
`sync(force = TRUE, upload = TRUE)` -- the packs were rebuilt from scratch and landed on the same
content address, which is the sharp proof that labelling moved and nothing else. (It was predicted
first: `systemsage_system_components()` is `identical()` old vs new, names, values and order, and
the packed matrix order still comes from `member_coef_files()`, deliberately untouched.) The 12
computed labels are identical to the pre-change baseline, set and pairing. **Parity passes** on the
post-change tree. `devtools::test()`: **731 pass / 0 fail / 248 skip**, +8 expectations from three
new tests exercising the default branch on GrimAgeV1/V2/PhysAge and the override on SystemsAge.
`R CMD check` was **not** run, per the standing invariant.

---

## 2026-07-27 -- housekeeping tail: MiAge memoized bit-identically, a scoring-failure channel, and the unbuilt verbs cut from the invariant

Three items off `dev/housekeeping.md`, chosen by the maintainer. `devtools::test()` after:
**723 pass / 0 fail / 248 skip**, from a re-verified 718/0/248 baseline (+5 expectations, both new
tests). Parity was **not** run and `R CMD check` was **not** run, per the standing invariants.

**Item 18 (MiAge) does not move a number after all, and is 3x not 1.65x.** The worklist proposed
hoisting `log(b)` *and* binding `b^(n-1)` as `d * log(b)`, measured 1.65x, and flagged it as
`max_abs = 2.76e-08` -- close enough to `PARITY_REL_TOL` to need a parity run, which is why it had
sat there. Folding `d` into the hoist is what cost the bit-identity, and it buys nothing: `log(b)`
is the expensive term and `d` is already a vector multiply. Hoisting **only** `log(b)` and leaving
the multiply order exactly as written (`res * bn * logb * d`) reproduces the original doubles
operation for operation. The bigger win was elsewhere: L-BFGS-B asks for the objective and the
gradient at the same `n`, and both need `b^(n-1)` *and* the residual, which the old closures
recomputed four times per point. Binding both once per `n` gives **2.85-3.04x** over four
independent unseeded draws with `identical()` TRUE and `max_abs = 0` every time. MiAge is ~40% of a
default run. No parity requirement, because there is nothing to verify.

**A scoring-time failure is now on the record, because coverage structurally cannot carry it.**
BMIQ calibrates a sample whose mixture will not fit to NA, so its score is NA while its coverage
reads a whole panel -- and coverage is computed *before* scoring precisely because it depends on no
score, so it can never learn this. On a saved record that NA was indistinguishable from a missing
covariate. `$provenance$scoring_failures` is now `clock id -> sample ids`, empty list in the normal
case, plus a `cli_warn` naming the clock and the samples.

The channel is deliberately **not** a second return value from a branch: the invariant that a
branch returns only its score is what keeps the dispatch `switch` one line per tag, and it is not
worth spending on a diagnostic. `mc_block()` carries a write-only collector (`new_notes()`),
branches write to it via `note_scoring_failure()`, and `score_cohort()` harvests it once per block
into a fourth return element. It composes under chunking for free -- each block harvests its own
rows -- and nothing routes or scores on it.

**The cross-package mock for this test is a trap: it works under `test_file()` and silently does not
under `test_local()`.** `local_mocked_bindings(bmiq_calibration = ..., .package = "betanorm")` hit
the mock and passed when the same file was run from outside `tests/testthat/`, and did not fire at
all in the real suite -- the test then failed in the one way that matters while passing every way it
was checked. Verbatim copy inside vs outside the test dir was the decisive experiment. The test
provokes a real failure instead: a constant sample has no three-component mixture, so
`DNAm[2, ] <- 0.5` fails BMIQ deterministically with no mocking anywhere.

**The five verbs CLAUDE.md asserted do not exist and are now cut.** `as.data.frame.mc_result`,
`[.mc_result`, `cbind.mc_result`, `augment` and `codebook` have zero occurrences in `R/`; the built
surface is `print`, `as.matrix`, `cite_clocks` and `rbind` (which exists to refuse). The invariant
read as a description of shipped behaviour and was a description of an intention, including the live
behavioural claim that `as.data.frame()` returns "scores plus the id column". CLAUDE.md and
`migration-plan.md` now name only what is built; `detail-plan.md` keeps the designs -- it is the
design doc -- behind a **Built** column, so a reader cannot mistake a spec for a contract. The
scores-only rule survives as a constraint on any future frame conversion rather than as a claim
about one that exists. Not touched: `migration-plan.md` still lists `get_clock()` /
`get_clock_probes()` in its API table and those do not exist either (already noted in the
2026-07-23 entry) -- same class of staleness, outside what was asked for here.

---

## 2026-07-27 -- housekeeping Waves 2-5: the `block` record, the narrow-once hoist, and four declarations the code was guessing at

Applied from `dev/housekeeping.md` items 8, 9, 11, 12, 13, 14, 15, 16, 17, 19, all of Wave 4, and
the tractable half of 10. `devtools::test()` after: **718 pass / 0 fail / 248 skip** (parity
auto-skipped, 1 network). `R CMD check` was **not** run, per the standing invariant. Item 18
(MiAge) was left alone: it is the only remaining item that moves a number, and that is the
maintainer's call.

**The maintainer then ran the parity tier: 0 fail / 30 skip / 657 expectations, 71.1s.** The 30
skips are exactly the known `horvath_online` block and the census test passed, so no fixture went
missing. That upgrades everything below from "argued bit-identical" to measured: every clock is
oracle-verified against both registry cohorts with the `block` record, the narrow-once hoist, the
policy-driven `absent_fill()`, the Dunedin delegation and the shared `sex_rows()` all in place.
Worth recording that the run *looks* like a 3x regression against the standing state and is not --
CLAUDE.md's "217" counts `test_that` blocks and testthat's `PASS` counts expectations, of which
`expect_parity()` carries three. That line now says so.

**`block` is now the one channel between the front end and every branch (8).** The
`(DNAm, partial_cache, pheno, packs)` quadruple threaded through 16 functions, each taking a
different subset in a different order, which is why the dispatch `switch` was argument re-typing
rather than dispatch. `mc_block()` builds the record once per block and every branch takes
`(id, cpgs, block, results)`. The switch is one line per tag. The one deliberate exception is
`linear_score(cpgs, block, observed)`: it is the shared engine, not only a branch -- `score_Dunedin`
and `score_normalized` call it too -- so it keeps the engine signature rather than carrying a
`results` argument it would never read.

**The narrow-once hoist is bit-identical, and that is measurable (17).** `mc_block()` narrows to
the panel union up front instead of letting each of ~100 clocks gather its own panel out of the
caller's full-width array. Measured over 86 bundled clocks, scores are **max_abs = 0** against an
unnarrowed run, and `clocks_coverage()` / `samples_coverage()` compare equal. The scaling is the
point: per-clock gathering costs 0.22s at 39k columns and 0.67s at 399k, while narrow-once is flat
at 0.22 -> 0.25s (2.7x at EPIC width, widening with the array and again under chunking).

**Narrowing has exactly one exception, and it is declared.** A `sample_scale` recipe op z-scores
each sample over every probe in the input matrix, so feeding Zhang2019 the panel union would move
its moments and its score. `block[["DNAm_full"]]` carries the caller's matrix untouched for that
one branch. The predicate is now `clock_needs_full_panel()` in the package (13) rather than
`"Zhang2019" %in% clock_ids`; the parity tier's own copy was deleted and aliased to it, so the
scorer and the fixture loader cannot disagree about which matrix a clock needs.

**A routing gap must fail at routing, not fall through (12).** `score_type()` returned `normalized`
for any scheme in `NORM_SCHEMES`, but `score_normalized()` implements only `bmiq`. The obvious fix
-- narrow the tag test to the implemented schemes -- is **wrong**, and worse than the bug: a
`quantile` clock outside the Dunedin group would then fall through to the `linear` arm and be
scored raw, silently dropping a normalization step the catalog declared. So the test still fires on
the full `NORM_SCHEMES` and `stop()`s when the scheme is not in `NORM_SCHEMES_ROUTED`. `unroutable()`
now reads all four routing facts (including the scheme) off the entry itself, so no call site can
report a subset of them.

**One `sex_rows()` for the four copies of the sex split, and the NULL case resolved (15).** Two
sites built index vectors (`which(female == 1)`), two built masks with a separate NA sweep, and the
four disagreed on a NULL/absent `Female`: the mask sites left member rows **unmasked**, reporting a
full per-sample miss count for a sample the member never scored. `sex_rows(female, n)` returns two
masks and treats unknown sex and no sex column alike as "neither" -- the same answer the alias
already gives, which is NA. That is the behaviour kept, because the alternative reports coverage
for a sample it is not true of.

**Four sites hard-coded an absent-CpG policy that `coverage_record()` reads from the catalog (14).**
`score_Dunedin` and `physage_raws` always vendor-filled, `score_GrimAge`'s internal surrogates
always dropped. All three happen to match the catalog today (Dunedin/PhysAge `vendor_mean`, GrimAge
`omit`), so a sync that flipped a policy would have moved the reported `score_used` /
`score_imputed_full` / `score_dropped` while the score silently did not move. `absent_fill()` now
owns the read for all four including `linear_score()`, and returns the term count as well as the
offset -- a dropped CpG has to leave a `mean` reduction's denominator too, which the old PhysAge
`/ length(coef)` did not do.

**`packs` was dead through three layers (9).** `score_type()` routes every external clock to a
`pack_*` tag before the `linear` / `normalized` checks, so `clock_coefs()`'s and
`clock_impute_ref()`'s external branches were unreachable from `calc_clocks()`. Rather than leave
an unreachable branch with a passing test, the two were settled differently: `clock_impute_ref()`
got a real caller (`pack_design()` was re-implementing its body inline, minus the missing-vector
diagnostic), and `clock_coefs()`'s external branch was **deleted** -- the batched pack scorers
multiply the whole coefficient matrix and never want one column. Deleting it also retired the
`identical(group_id, "SystemsAge")` payload dispatch inside it, which was the third bullet of item
13. Four tests that pinned the deleted branch went with it (722 -> 718 passes; no coverage lost,
the branch is gone).

**Item 11 is half a real finding.** `systemsage_pca()` guessed the rotation's key column
(`if ("system" %in% names(rot_df)) ... else names(rot_df)[1]`) where the component **declares**
`row_key: "system"` -- that is now `component_row_key()`, and an absent or mismatched declaration
stops. The second bullet does not survive contact: `raw_Blood` -> `Blood` is not a package-side
guess but an upstream naming convention that `data-raw/sync.R` applies itself to label the pack
matrices, and which it hard-fails on if the recipe steps and the group's member clocks disagree.
The component `name` the item points at (`sys_raw_Blood`) carries the same prefix, so driving off
it moves the strip rather than removing it. Nothing in the catalog declares the bare organ label.
The strip is now in one place (`systemsage_label()` over `SYSTEMSAGE_RAW_PREFIX`) instead of four,
and a real fix needs an upstream `label` field, not a downstream patch.

**GrimAge classifies stack operands by namespace, not by name (13).** The stack step declares
`inputs` / `internal` / `covariates` separately and `stack_operands()` concatenated them, discarding
the boundary; the branch then re-derived it from the operand's spelling (`nm %in% c("Age",
"Female")`, `startsWith("_internal_")`), so a covariate other than Age or Female would have fallen
through to the upstream-scores arm. `stack_roles()` keeps the boundary and `grimage_stack_roles()`
returns it in Cox-coefficient order, stopping on an operand the stack step never declared.

**Item 16 needed no new golden.** `score_Dunedin`'s non-QN path now delegates to `linear_score()`
as its first line. `test-score-dunedin.R` already derives the PoAm38 score independently from the
real coefficients with 2 CpGs absent, so it *is* the equivalence proof; it passes unchanged. What
remains in the branch is only what QN actually needs.

**Item 19 measured 11.3x, not the 32x the worklist guessed.** `resolve_cpgs()` called `match()`
once per distinct panel, rebuilding the hash of `usable` each time -- 26.65 ms over 89 panels and
39,025 usable columns, against 2.35 ms for one `match()` plus a `split()`. The `split()` factor
must carry all levels: most norm panels are `character(0)`, and dropping the empty groups would
shift every later panel's index. Output is `all.equal`-identical to the per-panel version.

**Item 10 landed for `accessors.R` and `resolve_inputs.R`, and is refused for `mc_data.R`.** The 18
clock-family accessors were contiguous at the end of `accessors.R` and moved next to their branches
(720 -> 533 lines); `resolve_inputs.R` split cleanly on its three phases into `validate_inputs.R` /
`resolve_inputs.R` / `coverage_gates.R`. `mc_data.R`'s claimed "15/13 exactly along the reachability
line" is **not** there: `mc_fetch`, `mc_consent`, `mc_ask_yes_no`, `mc_pack_paths`,
`mc_resolve_assets_dir` and `mc_resolve_groups` are all shared by `load_mc_assets()`,
`download_mc_assets()` and `clear_mc_assets()`, so the split would need a third shared file and
would not be the split the item describes. Left alone. Every move was checked by diffing the
top-level definition name set across all of `R/` before and after (210 both ways, zero difference).

**Wave 4 kept the `invisible()` sentinels but killed the computed ones.** The six side-effect-only
checkers now uniformly return `invisible(NULL)`. Dropping the `invisible()` entirely, as the
worklist suggested, would have made each function's return value whatever its last `if` block
happened to evaluate to -- a worse outcome than the unread constant. `clock_panels()`'s
`normalize = NULL` default was **kept**: the worklist called it unreachable from package code, but
three test files (including the parity tier) rely on it, so it is a used API, not a phantom one.

---

## 2026-07-27 -- housekeeping Wave 1: six correctness fixes, and the row gate moves onto the panel slice

Applied from `dev/housekeeping.md` Wave 1 (the only wave with correctness stakes). Most items are
"the fix the file already implied"; the three below changed observable behavior and are worth the
record. `devtools::test()` after: **717 pass / 0 fail / 248 skip** (parity auto-skipped, 1 network).
`R CMD check` was **not** run, per the standing invariant. The parity tier was **not** run either --
nothing here moves a score, but items 2 and 3 change which inputs are *refused*, so the maintainer
may want it before merging.

**1. The all-NA-sample gate now reads the panel slice, and `col_stats()` grew a `row_obs` field.**
This reverses the 2026-07-27 entry below, which kept the row half behind `anyNA()` on the grounds
that it scans the full matrix width. That framing was the bug: requiring `row_miss == ncol(DNAm)`
means **every** column NA, including the ~430k outside any scoring panel, so a sample with zero
observed CpGs on every panel sailed through and was imputed end to end from other people's cohort
means. The kernel now accumulates a per-row observed count over the submatrix it is already
walking, so the correct gate is also the cheaper one -- the `anyNA()` + `mat_miss()` full-width
passes are gone. The `length(present_needed)` guard keeps the diagnostic honest when the matrix
shares no CpG with any panel: that is the coverage gate's message to give, not this one's.

**2. `min_clocks_coverage = 0` is off only where a fully absent panel is scoreable.**
`ratio == 0 || ratio < threshold` short-circuited, so 0 was never actually "no gate" -- which
silently defeated the parity tier, since it scores at exactly 0. The `ratio == 0` arm is now
conditioned on the **declared** policy: `vendor_mean` fills every absent CpG from the vendored ref,
so an empty panel is scoreable by design and 0 must let it through; under `omit` every term goes
with its CpG and the score is undefined at any threshold. The existing "zero observed CpGs stops
even at 0" test kept passing untouched because all three clocks it names declare `omit` -- it was
testing the right thing for the wrong reason, and now says so.

**3. The gates label a routed member with its alias, and keep the member's sex.**
`check_coverage()` / `check_row_coverage()` walk the whole compute sequence, so a thin matrix on
`DNAmFEV1_wAge` was telling the caller to drop `DNAmFEV1_wAge_Female` from `clocks` -- a clock
`resolve_clocks()` hard-errors on. `gate_label()` puts the alias first. It does **not** collapse the
two members into one row, which is what the backlog item proposed: their panels differ (77 vs 73
here, 172 vs 190 for FitAge) and a count only means something against the panel it was counted over,
so merging them would be exactly the fabricated figure the routed-members invariant forbids. The
line reads `DNAmFEV1_wAge (female model): 23/77 scoring CpGs`, which is both actionable and true.

**4. Both kernels rewritten: raw column pointers, one predicate in the hot loop, no early exit.**
Measured at `-O2` on a 100 x 39,025 slice (an EPIC panel union, 2% partial NA), ms per call, best
of 5 x 300:

| `col_stats()` | | `fill_imp_col()` | |
|---|---|---|---|
| `obj(i, j)` + `R_finite` (was) | 7.13 | `obj(i, j)` + `R_finite` (was) | 5.50 |
| column pointer + `std::isfinite` | 3.73 | column pointer + branch | 1.60 |
| one predicate, Inf caught per column | **2.70** | branchless select on `isnan` | **1.23** |

Net **2.6x** and **4.5x**. Three separate changes, and the order they mattered in is not the order
one would guess:

- **Raw column pointers** (`double *col = x + j * nr`) instead of `obj(i, j)`, which recomputes the
  offset per element -- and which `fill_imp_col()` was paying twice per missing entry, once to test
  and once to assign. Worth ~1.3x on its own; at `-O2` gcc already hoists most of the index
  arithmetic, so this is the smallest of the three.
- **`std::isnan` / `std::isfinite` instead of `ISNAN` / `R_finite`.** These are the big one for
  `col_stats`: the libR predicates are **out-of-line calls**, so nothing across them vectorizes.
  Isolated on a bare sum loop: 6.0 -> 4.5 ms from pointers, 4.5 -> 1.5 ms from the predicate.
- **One predicate, and the `Inf` check hoisted out of the loop** (below).

**Only `fill_imp_col()` wanted a branchless form.** The test per element is unavoidable; the
*branch* is not, and a branch is what stops a loop vectorizing. Writing the fill as
`col[i] = isnan(v) ? m : v` is 1.6 -> 1.2 ms. The same rewrite on `col_stats` measured 2.70 -> 2.60
ms, i.e. nothing -- its loop carries five dependent accumulators, so the compiler was already doing
what it could. It is not worth the readability, and was not taken. (An earlier draft of this entry
claimed 5x for that rewrite and 9x for the fill; both were timer quantization. `system.time` on a
0.5 ms call is noise, and the fill benchmark was re-running on an already-filled matrix, so the
branch never fired. Numbers above use a 300-iteration inner loop and re-poke the NAs each rep.)

**The loop now asks one question, not two, and that is the interesting part.** There are two
different questions to ask of a value -- "is it missing?" (`isnan`: `NA_real_`, `NaN`) and "is it
usable?" (`isfinite`: also excludes `+/-Inf`) -- and the old kernel asked both per element to get a
three-way split. It only needs the first. An `Inf` is let through into `sum`, poisoning it, and one
`isfinite(sum)` per **column** catches that: `nr * nc` tests become `nc`, and the early `return` --
which also blocked vectorization -- leaves the inner loop entirely. Locating the row is a cold
rescan of the single bad column, at most once per call and on the way to an abort. The
classification the caller depends on is unchanged; it is established per column rather than per
element.

**That opened one new case, and it is handled rather than ignored.** `!isfinite(sum)` also fires on
a finite column whose sum **overflows**. That needs values around 1e306, which no beta matrix has --
but the value gate only *warns* on out-of-range input, so nothing stops one arriving. The rescan
tells them apart: no `Inf` found means overflow, the kernel reports row `NA`, and
`check_col_values()` words it as a separate abort naming the column instead of inventing a sample.
Both are hard stops. Tested at the kernel and through `calc_clocks()`.

**None of this changes semantics**, and the reason is that R requires IEEE-754, on which the
predicates agree with the macros exactly. The one thing that would break it is `-ffast-math`, which
CRAN forbids and neither `Makevars` sets. Loop counters also moved to `R_xlen_t`, which stops
`j * nr` overflowing `int` past 2^31 elements -- reachable at ~5,000 samples x 450k probes.

**Worth knowing, but it does not apply where you would expect: `pkgbuild::compile_dll()` builds at
`-O0`.** It appends `-UNDEBUG -Wall -pedantic -g -O0` after R's `-O2` and the last `-O` wins, so
every `devtools::load_all()` session runs an unoptimized `src/` -- ~4x slower on these kernels.
**This does not taint the Wave 5 profile**, contrary to a first draft of this entry: MiAge is pure
R (`stats::optim` over R closures, `R/score_MiAge.R`) and `qnorm_target_rows_cpp` belongs to
`betanorm`, an installed package built at `-O2` by its own install. `compile_dll(".")` rebuilds only
this package's two kernels, and neither appears as a line in that profile. It matters only when
benchmarking `col_stats` / `fill_imp_col` themselves -- as here, which is why every number above
was taken through `sourceCpp` at an explicit `-O2`.

**Also in this pass**, without behavior change worth arguing about: `finalize_PhysAge()` now refuses
a constant surrogate instead of returning `NaN` for the entire cohort (`rowSums(scale(.))` propagates
one dead column to every sample; the `n < 2` guard folded in beside it, since both are "scale()
cannot z-score this"), `recipe_step_out()` stops on >1 matching step rather than folding ambiguity
into the `NULL` that means "no such step" -- eight clocks take their `Age` term off that fallback, so
the old behavior would have scored them with no covariate term -- two `cli` pluralization crashes on
the closed-set assets path are bound with `cli::qty()` (the handler was throwing **in place of** the
diagnostic, a regression against 2026-07-23), and `digest` is declared in `Suggests` with a
`skip_if_not_installed()` at the top of `test-mc_data.R`, where 12 of 16 always-on tests were relying
on it being incidentally installed.

---

## 2026-07-27 -- value gates ride the `col_stats()` sweep; the `anyNA()` short-circuit was hiding an `Inf`

**The defect, found in an audit of the Phase 5 diff.** `scan_missing_cpgs()` returned early on
`!anyNA(DNAm)`. `anyNA()` does not see `+/-Inf`, but the new `col_stats()` kernel treats every
non-finite value as missing -- so whether an `Inf` got filled depended on whether an **unrelated** NA
existed somewhere else in the matrix. Measured on one Hannum request: the same `Inf` scored `-Inf`
with no NA present, and a plausible `60.00852` once an NA was added in a different column. The
pre-Phase-5 `slideimp::mat_miss()` counted only NA, so the old code was consistently the former;
Phase 5 introduced the *inconsistency*, and it was silent and data-dependent, which is the bad kind.

**Decision.** The panel scan runs **unconditionally**, and `col_stats()` returns a list rather than a
statistics matrix: `stats` (2 x ncol, named rows `sum` / `n_obs`), the global flags `any_lt0` /
`any_gt1`, and `inf_at`. Only the row half (the all-NA-sample abort) stays behind `anyNA()`, because
it scans the full matrix width rather than the needed panel and is genuinely skippable when there is
nothing missing.

**Fail fast on `Inf`, and make the caller's obligation explicit.** The scan stops at the first
`+/-Inf` and returns its 1-based `c(row, col)`; `stats` is then `NULL` and the range flags cover only
the prefix scanned. Nothing in the shape of the list enforces that, so it is a stated contract:
`check_col_values()` tests `inf_at` first and reads nothing else in that branch. Cheaper than
finishing a traversal whose result is about to be thrown away, and it is what lets the abort name the
exact sample and CpG.

**Asymmetry in what the gates report, on purpose.** The abort carries a position; the warnings carry
none -- not a column list, not the observed range. An `Inf` is a point defect worth locating and the
scan is already stopping there. Out-of-range values are a property of the whole matrix, where naming
ten columns out of 866k is noise rather than information, so two booleans suffice and the per-column
`min`/`max` an earlier cut carried was dropped. Two flags, not one: a matrix can be out of range on
both sides and neither implies the other, so they warn separately.

**Missing and bad are counted apart.** NA and NaN are "not observed" and fill from the cohort mean.
`+/-Inf` is bad data that no fill can repair, so it gets its own counter and its own outcome: an
abort. Folding it into the NA count -- which is what the kernel was effectively doing -- means
silently imputing over a number that is wrong rather than absent.

**The kernel reports; R decides.** `col_stats()` raises nothing. `check_col_values()` (R) issues both
messages, so they go through cli with `call = NULL` per the CLI rule, and can name the offending
columns. `Rcpp::stop()` would have been fewer lines and the wrong shape.

**Range: warn, do not stop, and let 0 and 1 through.** A value outside [0, 1] is the beta-vs-M-value
tell -- `check_DNAm()` checks `mode = "double"` but not the range, so an M-value matrix scores
plausible garbage today. It is strong evidence rather than a definition, so it warns and still
returns scores. Exactly 0 and exactly 1 pass: those are ordinary saturated betas, and rejecting them
would be a false positive on real data.

**`n_obs == 0` is answered by classification, not by a new error.** The instruction to guard the
all-missing case explicitly was taken as "make it impossible and say so" rather than "abort": an
all-NA column is an ordinary expected case that lands in `all_na_cols`, leaves `usable_cols`, and
takes the clock's vendored ref or the drop policy -- the standing imputation invariant. Only
`partial` columns are divided and `n_obs > 0` holds for every one by construction, so `0/0` never
arises. Aborting instead would refuse any matrix missing a single probe and would contradict
CLAUDE.md's imputation rule.

**`fill_imp_col()`'s Inf branch is now unreachable and stays anyway.** `col_stats()` scans the same
columns earlier and aborts on the first `Inf`, so only NA/NaN reach the kernel through the supported
path. The `!R_finite` test is kept: filling is the right response to NaN, and narrowing it would make
the kernel's correctness depend on a gate in another file. Its comment says so; the in-place mutation
contract and the fresh-allocation requirement are untouched.

**Also dropped: `scan_missing_cpgs()$has_na`.** Set on every call and read by nothing.

**Not verified against parity.** The always-on tiers are green with **zero** warnings, so no
synthetic input trips the range gate. Real cohort betas should not either (0 and 1 are allowed), but
if a staged fixture carries an out-of-range value the parity run will now emit a warning where it
previously ran clean.

---

## 2026-07-27 -- Phase 3: the cohort reduction leaves the scoring loop, and `rbind` still must not redo it

**Decision.** `score_PhysAge()` split at its `cohort_zscore` step (recipe position 11) into
`physage_raws()` -- the per-sample surrogate matrix, computed in the scoring loop -- and
`finalize_PhysAge()` -- the reduction, run once after assembly. `score_cohort()` gained a third
return element, `pending`, and both front ends end in `finalize_cross_sample()`.

**The loop routes on the catalog, not on a tag or a clock id.** A branch's output goes to `pending`
instead of `scores` iff the clock is in `spec$cross_sample`, which `mc_spec()` derives from the
declared `cross_sample_at`. So the `switch` over `score_type()` stays total and unchanged in shape,
no branch decides its own deferral, and a clock that starts declaring `cohort_zscore` upstream defers
with no downstream edit. This is the amendment to CLAUDE.md's "a branch returns only its score": a
declared cohort-reducing clock returns its per-sample intermediate instead.

**`finalize_cross_sample()` is called unconditionally, which is the whole point.** `pending` is empty
for 127 of 129 clocks, so it is a no-op for a typical request. That is what keeps the public surface
magical: `calc_clocks()` and (later) `calc_clocks_chunked()` return the same number, and nothing
anywhere -- not a branch, not a front end, not the record -- tests whether the run was chunked.

**Carry the raws; do not stream sufficient statistics.** Per-surrogate `(sum, sumsq, n)` accumulated
per block would let the first z-score stream, but the raws are 8 doubles per sample (~640 KB at
n = 10,000), so carrying them costs nothing and avoids a second accumulator that would have to keep
bit-agreeing with the first. Decisive point: it would not have worked anyway for
`DNAmPhysAge_years`, whose **second** `cohort_zscore` reduces over `phys`, itself a function of the
first reduction -- not streamable, trivial with the raws in hand.

**The `n >= 2` guard moved to the finalize.** It was in the branch, where it would have rejected any
1-row block under chunking. It is a statement about the *cohort*, so it now sits where the cohort is
whole: a 1-row block scores its raws fine, a 1-row cohort still errors.

**Rejected: re-finalizing at `rbind`.** The natural-sounding "after rbinding, PhysAge is recalculated
over the combined statistics" is right for chunk assembly and wrong for `rbind.mc_result`, and they
are different layers. Chunk blocks are one cohort whose reduction has not run yet -- finalizing after
assembly is the *only* calculation, not a recalculation. Bound records are several cohorts that each
already finalized and whose numbers the user has already seen; re-reducing over the union would
silently rewrite them and would standardize across batches that were deliberately imputed
differently. So assembly (Phase 6, fragments, sums coverage) and `rbind` (Phase 4, records, keeps
coverage per record) do not share an implementation, and `rbind` binds and labels exactly as sec 8
said. Corrected on the way: `dev/id-streaming-plan.md` sec 8 claimed "after Phase 3 no score column
is cohort-dependent", which was never true and contradicted its own closing paragraph -- Phase 3
relocates the reduction within a run, it does not make the column batch-independent.

**Phase 6 does not depend on Phase 4.** The chunked front end never calls `rbind.mc_result`, so
`rbind` can go on refusing while chunking ships.

**Measured, both bit-exact.** A 3-block split now finalizes to the single-pass answer for both
PhysAge clocks -- previously off by 2.8 and **11 years** respectively, purely from re-standardizing
against 10 rows instead of 30. And single-pass output is `identical()` to the pre-Phase-3
implementation, so the standing parity run is undisturbed and does not need repeating for this
change. The chunk-invariance test now compares **every** clock with no exclusion, and asserts
`length(spec$cross_sample) > 0` so it cannot quietly stop exercising the deferred path.

---

## 2026-07-27 -- Phase 2: a sex-routed alias derives its sample axis instead of assuming it

**The defect.** `attach_sex_routed_aliases()` (`data-raw/sync.R`) minted every alias with a literal
`cross_sample_at = NA_integer_`. Every *other* alias field is derived -- `normalization`,
`imputation_policy`, `pmid`, `license` from the donor, `clock_inputs` from routing -- so this one
constant was the only assumption in the block, and it asserted the strongest possible claim about
the members without consulting them.

**Why it matters more than it looks.** `batch_dependent` is `!is.na(cross_sample_at)`, so a
cross-sample member under a per-sample alias would tell a chunked engine the alias may be scored
block by block. That is a **wrong score, not a missing one**, and it cannot be caught downstream:
the alias's rows would look complete and plausible. Phase 2 is the point where the engine starts
trusting the field, which converts a latent hole into a live one.

**Fix:** `alias_cross_sample_at(members)` -- the min over the two members, NA when neither reduces.
Positions are not comparable across two different recipes, so the min is read as "the first index
past which the alias cannot be assumed per-sample" rather than as a step pointer.

**No data moved.** All 14 routed members are `DNAmFitAge` and per-sample, so the derived value equals
the old constant for all 7 aliases and `R/sysdata.rda` needs no regeneration. Verified against the
built catalog, not assumed. The rest of the sync-side classification was audited at the same time and
is correct: `cohort_zscore` is the only cross-sample op in the whole upstream corpus (3 occurrences,
2 clocks); `sample_scale` (Zhang2019) is a within-sample z-score; and `center_scale` (SystemsAge)
reads vendored `systems_pca_center` / `systems_pca_scale` tensors, so it is per-sample too.

**Package side, no clock lists anywhere.** `clock_cross_sample_at()` / `clock_is_cross_sample()`
read the declared field, `split_cross_sample()` partitions a compute sequence, and `mc_spec()` stores
the reducing half as `spec$cross_sample`. The chunk-invariance test excludes exactly that vector
instead of a literal, so it empties itself when Phase 3 lands; a census test walks
`resolve_clocks("all")` and asserts the partition is total and that each alias matches its members.
That last assertion is the standing guard against this decision being quietly reversed -- it compares
an alias to its members rather than to a constant, so re-hardcoding NA fails the always-on tier.

---

## 2026-07-27 -- Phase 5: the seam is cut, `partial_fill` is the channel, and `fill_imp_col()` stays in place

**Decision.** `calc_clocks()` is now `mc_spec()` + `mc_cohort()` + `score_cohort()` +
`construct_mc_result()` (`R/score_cohort.R`), split by what each depends on, per
`dev/id-streaming-plan.md` sec 4. Nothing about the public surface moved.

**`facts` carries `partial_fill` (cohort means), not a prebuilt cache.** The alternative -- having
`mc_cohort()` build the filled sub-matrix once and hand it to `score_cohort()` -- is a smaller diff
and carries zero numeric risk, but it makes `score_cohort()` take a whole-cohort object, which is
exactly the thing Phase 6 would then have to re-cut. A seam that a later phase must re-cut is not a
seam. So `score_cohort()` builds its own cache per call from the means, which is what lets Phase 6
add a front end and two adapters instead of reopening this. The names of `partial_fill` are the
column classification, not just a lookup key: a block must never re-derive which columns are
cohort-partial, because a column that is partial cohort-wide can be entirely NA inside one block.

**Consequence: `slideimp::mean_imp_col()` is out of this path.** Splitting mean-computation from
filling means the mean has to be available on its own, so `scan_missing_cpgs()` now runs the
package's own `col_stats()` -- one traversal returning per-column `(sum, n_obs)` -- and derives
`all_na_cols`, `partial_na_cols`, `usable_cols` and the fill values from that single sweep. Two
slideimp calls collapse into one, and the classification and the fill can no longer disagree because
they come off the same pass. `slideimp::mat_miss(col = FALSE)` keeps the row half (the all-NA-sample
abort), which is per row and correct against any block. Semantic change worth knowing: `col_stats()`
skips non-finite values, so NaN and +/-Inf count as missing where `mat_miss()` counted only NA. A
beta matrix has none of the three, and an Inf in one would be bad data worth filling anyway.
**Parity was not run** -- this moves the numeric path for every clock with partial NAs, so it should
be run before merge.

**`fill_imp_col()` stays in-place and void, reversing sec 5.2.** The plan said the kernel "must slice
and fill in one traversal, returning a new matrix", on the grounds that an in-place kernel that ever
sees the caller's `DNAm` corrupts it invisibly. That risk is real but it is bounded by one caller:
`build_partial_cache()` always passes a fresh `DNAm[, cols, drop = FALSE]`, and R matrix subsetting
always allocates, so the caller's matrix is never reachable from the kernel. Making the kernel copy
would allocate the slice twice for a guarantee the single call site already provides. Accepted cost:
the safety property is now a **caller contract rather than a kernel property**, so it is written down
in `src/fill_imp_col.cpp` and at the call site, and a second caller must re-establish it. Revisit if
one ever appears.

**Two orderings moved, both deliberate.** `check_row_coverage()` runs after the scoring loop now
(sec 4.1: warn-only, reads assembled counts). And `load_mc_assets()` is in the data-independent tier,
so it precedes `check_DNAm()` and the pheno checks -- a caller asking for an external pack with a
malformed matrix gets the download prompt first and the error second. That is a genuine small
regression, accepted because hoisting resolution out of the per-block path is the point of the tier
split; the fix would need a pre-check that does not re-warn when `mc_cohort()` repeats it.

**The Rcpp wiring lives in a roxygen block, not in NAMESPACE.** Calling a compiled kernel from R
needs `useDynLib(methylCIPHERv2, .registration = TRUE)`; without it the kernels build fine and
`.Call()` fails with "not available for .Call() for package". Adding it by hand to NAMESPACE does not
survive -- NAMESPACE now carries roxygen2's "do not edit by hand" header and `document()` regenerates
it from tags, silently dropping anything hand-added (observed: a hand-added pair vanished on the next
`document()` and took the whole compiled path down with it). So the two lines come from
`R/methylCIPHERv2-package.R` (`@useDynLib` + `@importFrom Rcpp sourceCpp`), the standard
`usethis::use_rcpp()` block. This partially reverses CLAUDE.md's "no roxygen yet" and "never edit or
regenerate NAMESPACE": roxygen is on for whatever Rcpp requires, and NAMESPACE is now a generated
file, so the maintainer owns the **tags**, not the file. CLAUDE.md still states the old rule.

**Measured, and better than the doc promised: bit-exact chunk invariance.** Sec 4.2 only claimed
"equal, not identical", expecting the mean's association order to move the last bit. Over a 3-block
split of a cohort carrying all three missingness shapes, every chunk-safe clock agrees exactly --
because Phase 5 computes each mean once and every block reads the same number, so there is no
re-summation to reassociate. This does **not** carry to Phase 6, where pass 1 accumulates block by
block; the doc's warning stands for that. `DNAmPhysAge` is the sole non-invariant clock (~1.1 off),
exactly as sec 6 predicts, and Phase 3 is what fixes it.

---

## 2026-07-27 -- `normalize=` is an argument, not `prep()`; BMIQ defaults off and drops where QN fills

**Decision.** Per-clock normalization is controlled by `normalize=`, a named logical argument on
`calc_clocks()` and `sim_DNAm()` -- `normalize = c(Horvath1 = TRUE)`, or a bare scalar policy.
`resolve_normalize()` (`R/resolve_inputs.R`) turns it into one decision per clock in the compute
sequence, before any DNAm is read. Constitutive schemes are on and cannot be declined; everything
else defaults **off**.

**This decouples normalization from `prep()`, and retracts the framing above.** The entry below has
`prep()` as "where per-clock normalization is decided". That conflated two features: a *policy*
(a small named vector) and a *data structure* (the per-clock prepared matrices, potentially 24.7 GB).
The policy needs no record -- it is an argument that flows to where normalization already happens,
inside the clock's branch. `prep()` remains worth building as an **audit** feature ("what number did
Horvath1 actually multiply?"), but nothing about normalization control waits on it, and `prep()`
takes the same `normalize=` argument when it lands. Building them together would have blocked a
2-day change on a multi-week one.

**Constitutive vs declinable is read off the scheme, not a new upstream field.** The entry below
asks upstream for that declaration. Upstream shipped the narrowing (below) but **no such field** --
Horvath1 and Knight declare `normalization: ["BMIQ"]` exactly as DunedinPACE declares
`["quantile"]`. Resolving it from the scheme is acceptable where a clock-id exception would not be:
`NORM_CONSTITUTIVE` is a value in a closed enum that upstream owns, so a new BMIQ clock inherits the
answer with no downstream edit, and adding a third scheme is a one-line change next to the other
three `NORM_*` constants. If upstream ever needs a per-clock override of the scheme-level answer,
that is when the field becomes necessary.

**Upstream state, replacing the census below.** meta `a7148b3` narrowed `normalization` from 7 BMIQ
clocks to **2**: EpiTOC, EpiTOC2, HypoClock, Mayne and PedBE went to `["none"]` (generic type-2
corrections and training-side sentences, not scoring contracts), leaving Horvath1 and Knight. It
also **vendored the calibration target**, reversing 2026-07-26's "the BMIQ gold standard is
deliberately NOT vendored": `goldstandard2.csv.gz`, 21,368 probes, new `bmiq_gold_standard`
probe_set role, committed once per group. Verified downstream: both scoring panels are strict
subsets of that panel (353/353, 148/148), and the two vendored copies are `identical()` as R
objects -- same probes, same order, same values. So the panel union across both clocks is exactly
21,368 and there is nothing to union; `dedup_panels()` already collapses them.

**Default `FALSE` for Horvath1 and Knight, deliberately, with a known cost.** `TRUE` would match
Horvath's calculator; upstream measured a **7.7 year** gap without calibration. `FALSE` wins anyway
on three counts: it is a no-op on existing behavior, so no silent numeric change; BMIQ needs 21,368
probes present to fit a meaningful mixture and most callers on targeted panels have ~353, where
defaulting on would silently produce a bad calibration; and it is what the ecosystem does (biolearn
defines these with no preprocess, methylclock's `normalize` argument is dead code). It also keeps
**Knight's parity untouched** -- its fixture is a `frozen_reference` on raw betas, so a default of
`TRUE` would have turned two green `core` rows red. Reversing the default later is one constant.

**BMIQ drops absent probes; QN fills them. This asymmetry is the point.** `score_Dunedin()` fills
fully-absent norm CpGs with the gold mean itself and then quantile-normalizes the full panel.
Copying that for BMIQ would be wrong, not merely different: BMIQ estimates the **sample's** own
3-state mixture from the panel, so injecting values drawn from the **target** distribution pulls the
sample's fit toward the gold standard and shrinks the correction being computed -- silently, and
more the more probes are missing. Dropping keeps the fit on real observations. Absent scoring CpGs
are then handled by the clock's declared `policy: "omit"` exactly as today, and partial NAs are
already filled from the cohort-mean cache upstream, so `bmiq_calibration()` never sees an NA.
Consequence: the thin-background warning's "Missing background CpGs are filled from the reference
mean" is false for a BMIQ clock and must become scheme-aware.

**No tuning surface.** `normalize=` takes a logical, never a per-clock options list. Every BMIQ
setting is fixed at Horvath's published choice (`nL = 3`, `doH = TRUE`, `th1.v = c(0.2, 0.75)`,
`niter = 5`, `tol = 0.001`), with `nfit = ncol(panel)` the single departure -- fit on every present
probe rather than a 20,000 subsample, since the panel is only 21,368 and a caller who wants
different preprocessing can do it to their own matrix before scoring. Determinism is not at risk:
`betanorm::draw_fit_indices()` saves and restores `.Random.seed` around a fixed seed, so scoring
neither consumes nor perturbs the caller's RNG stream.

**It adds a branch rather than collapsing one, and that is the drop-vs-fill decision showing up in
the routing.** The plan predicted `score_Dunedin()` would generalize into one normalize-then-linear
branch serving both schemes. It cannot: the two schemes disagree about what an absent probe *is*.
Rerouting DunedinPACE through the shared linear engine would change its numbers whenever a scoring
CpG is absent -- today it fills that CpG with the gold mean and normalizes it, where the linear
engine would drop it or take a vendored offset -- so it would put a currently-passing parity clock
at risk to serve a tidiness goal. `score_Dunedin()` therefore stays as-is and `score_normalized()`
is a second pre-transform branch, tagged `normalized` and routed by declaration
(`clock_norm_scheme(p) %in% NORM_SCHEMES`, catalog-only so the smoke tier still proves totality).
What *is* shared is everything after the transform: `linear_predictor()` / `linear_score()` gained
one optional `observed` argument, so the normalized branch reuses the policy handling, covariate
terms, mean reduction and output transform instead of reimplementing them.

**Row-level failure is now a third kind, and `samples_coverage()` does not yet report it.** BMIQ runs
`on.sample.error = "continue"` with `failed.sample = "NA"`, so a sample whose mixture fit fails gets
an NA score while the cohort completes -- deliberately, since aborting a run over one sample is
harsh and an NA is visible. But a reader then cannot tell *why* a row is NA: missing pheno, excessive
NA, and now a failed calibration all land the same way. `betanorm::bmiq_calibration()` returns a
per-sample `$success` vector, which is the signal a future `samples_coverage()` row-reason column
should carry.

**The `horvath` parity block stays skipped, and BMIQ does not change that.** `Horvath1@cohort_450K`
looked like the one fixture where calibration is the *sole* difference from the oracle -- it
declares no `missing` block, so zero absent scoring CpGs -- and therefore the one candidate for
un-skipping now that `normalize = c(Horvath1 = TRUE)` exists. The measurement already on record
(detail-plan sec 12) closes it: `betanorm::bmiq_calibration()` improves that pair **140x, to
1.9e-3**, which is seven orders of magnitude outside `PARITY_ABS_TOL["horvath"]` (`1e-10`). Moving
that bound would be an argument about algorithmic agreement, and the only licence this table has is
about units. So BMIQ's numeric gate is the always-on golden in `test-normalize.R`, not parity --
which is a thinner gate than a fixture and should be understood as such.

## 2026-07-27 (later) -- `prep()` is a lazy per-clock plan; normalization is decided per clock, from declarations

**Decision.** Add `prep(DNAm, clocks, ...)` (plan Phase 7): an S3 record, list-like and keyed by
clock id, whose element for a clock is the matrix that clock's coefficients multiply -- its panel,
imputed clock-wise and panel-wise, with its normalization applied or not. `calc_clocks()` accepts a
`DNAm` or a `prep`. This invents no boundary: `score_cohort()` already factors as `score(prep(...))`,
so Phase 7 makes an existing seam public.

**Lazy, and that is forced by measurement.** `Sigma |panel_i|` against the union is 53,127 / 39,025
(**1.36x**) over the 101 bundled clocks but 3,083,623 / 378,363 (**8.15x**) over the 28 pack clocks --
PCBrainAge alone declares 357,852 CpGs, the 13 SystemsAge organ clocks 125,175 each, the 14 PCClocks
78,464 each. Eager materialization at n=1000 is 425 MB bundled and **24.7 GB** external, against a
~6.9 GB full-EPIC input: 3.6x the matrix being prepared. So the record holds a plan (DNAm reference,
the shared cohort-mean cache, per-clock panel spec, per-clock normalization decision) and
materializes one clock at a time; peak RSS is unchanged from today. `materialize()` is the opt-in
eager form, which the bundled set can always afford.

**Why `prep()` does not return one prepared matrix.** Neither form is well defined. Imputation is
per clock -- a fully absent CpG takes the clock's own vendored ref as a scoring CpG and the
gold-standard mean as a normalization CpG (`score_Dunedin.R`'s `fill_ref` split), so the same probe
gets different numbers depending on which clock asks. Normalization is per clock -- different panels,
different targets, and normalizing the shared input would silently change what every non-normalizing
clock is scored on. The only line that exists is cohort-global vs clock-specific: cohort means are
hoisted, everything else stays in the branch. `prep()` is a list because of that, not in spite of it.
This also disposes of a `prep()`-as-pipeline-stage reading, in which imputation would have to be
hoisted out too -- there is nothing hoistable left to hoist.

**Normalization state, which is worse than it looked.** `normalization` is non-`none` on 9 clocks
(`BMIQ` x7: EpiTOC, EpiTOC2, Horvath1, HypoClock, Knight, Mayne, PedBE; `noob` x1: Horvath2;
`quantile` x1: DunedinPACE), but **DunedinPACE is the only clock in the catalog with a normalization
probe_set and a target file pointer**. BMIQ is declared on seven clocks and implemented nowhere;
`normalization: "BMIQ"` is a bare string with no artifacts behind it. `noob` is moot -- an
IDAT-intensity background correction, always done upstream, unreachable from a beta matrix.

**Horvath1's BMIQ is the DunedinPACE shape, not a whole-array transform.** Per the maintainer it
picks ~21k probes (a superset of its 353) and normalizes each sample against a vendored golden-mean
vector -- structurally identical to PACE's declared 20,000-CpG background plus gold target, and
matching `betanorm::bmiq_calibration(datM, goldstandard.beta, ...)`, which already takes exactly that
`(matrix, target)` shape. So it generalizes the existing branch rather than needing a new category:
`dunedin_gold_means()` becomes a scheme-agnostic target accessor and the branch becomes
normalize-then-linear dispatching on `clock_norm_scheme()` -- the existing **pre-transform** kind, no
new entry in the closed branch set. An earlier framing of BMIQ as an unavoidably whole-array
operation was wrong for Horvath1 and is retracted.

**Two facts upstream must declare, per clock.** (1) Is the normalization expressible as a declared
panel plus a vendored target? If not (genuine whole-array BMIQ needing type I/II annotation) it is
not expressible in the catalog at all and stays the caller's preprocessing. (2) Is it **constitutive
or published preprocessing**? PACE's QN is part of the clock's definition and cannot be disabled;
Horvath1's BMIQ is what the published pipeline did and a caller may decline it. Without that field,
"PACE cannot be turned off" becomes a hardcoded clock-id exception, which the declared-routing
invariants exist to prevent.

**No chunking cost.** `quantile_norm()` is per-sample against a fixed target and measured
bit-identical on a row subset (max abs diff 0); BMIQ fits its mixture within a sample. Neither adds
to `cross_sample_at`, so Phases 2/3 are untouched. When a BMIQ clock lands, Horvath1 goes from
needing 353 CpGs to 353 to score plus ~21,000 to normalize -- the shape the norm-panel coverage axis
and the 2026-07-23 gate split already handle, at which point CLAUDE.md's "(only DunedinPACE)"
parenthetical about norm panels goes stale.

---

## 2026-07-27 -- chunking passes a fill vector inward; pre-imputation is rejected, and `expect_identical` is banned

**Reverses 2026-07-26**, below, on its central mechanism. That entry had the chunked wrapper
pre-impute each block from cohort means and hand a clean matrix to the scorer, and concluded that
coverage therefore had to be hoisted out of `calc_clocks()` into the wrapper's first pass. Both
halves are wrong, and the second follows from the first.

**Decision.** The shared scoring internal takes `partial_fill`: a named numeric vector whose names
are the cohort-partial columns and whose values are the cohort means. Blocks are handed over
**raw, NAs intact**, and the fill happens inside, in the cache that `mean_imp_col()` builds today.
Public `calc_clocks()` calls the internal with `partial_fill = NULL` -- derive it locally -- so the
complete-cohort precondition remains literally true of the public surface and no chunk-aware
parameter appears on it.

**Why pre-imputation was wrong.** It is score-equivalent and that is exactly the trap: it erases the
NA pattern, and `count_sample_miss()` reads the NA pattern, not the fill. A pre-filled block reports
`score_imputed_partial = 0` for every clock, an all-zero `sample_miss`, and a `check_row_coverage()`
that reduces to a per-clock constant -- a silent breach of "coverage is never reported for a sample
it is not true of". Passing the vector inward instead keeps the NAs visible exactly where coverage
counts them, so **the coverage machinery does not move at all** and the hoist that the prior entry
built its case on is not needed. Assembly is a concatenate of the per-sample vectors plus a sum of
the two `imputed_partial` scalars; every other `per_clock` field derives from the shared `cpg_list`
and is already identical in every block.

**The vector's names are the classification, and that is the load-bearing half.** A column that is
partial cohort-wide can be entirely NA inside one block; `col_miss == nr` would call it structurally
absent and route it to vendor-ref-or-drop. So when `partial_fill` is supplied the column half of
`scan_missing_cpgs()` is **bypassed, not overridden** -- one source of truth, never a local
computation a caller corrects afterwards. The row half (the all-NA-sample abort) is per-row, correct
against a block, and stays.

**What the seam is now.** Still a split of `calc_clocks()` into `mc_spec()` / `mc_cohort()` /
`score_cohort()`, but justified by sharing and by cohort-set predicates (id uniqueness, the pheno
union, column classification, means, the clocks gate) -- not by coverage, which no longer has a
problem to solve.

**Phase 0 (`partial_cache` becomes a mean vector) is deleted, not deferred.** It was justified as a
memory win and measurement does not support it: against the current `mean_imp_col()` cache, a mean
vector with per-clock refill ran 2-5x slower with peak RSS within 1-2% (n=2000-8000, union
8k-30k probes, 40 panels). The persistent `n x k` cache is real but is not the binding constraint
next to `DNAm` itself, and the per-clock slice is a fresh allocation in either design, so nothing is
saved by not materializing the cache. Measured overlap across the 101 bundled panels is only 1.36x
(union 39,025 CpGs, `Sigma |panel_i|` 53,127), so the amortization the cache buys is modest -- but
modest and free beats a regression.

**Two kernels are ours, in this package's C++.** The package is going Rcpp so the `betanorm`
normalization functions can be vendored in and the dependency dropped; `fill_imp_col(obj, values,
subset)` and `col_stats(obj)` (per-column sum / n_non_na / n_na / min / max in one traversal) land
there rather than in `slideimp`. `fill_imp_col()` must **slice and fill in one traversal returning a
new matrix**, never mutate in place: the slice is already a fresh allocation, so in-place buys
nothing, and a kernel that can reach the caller's `DNAm` corrupts it invisibly. A pure-R fill is the
2-5x regression measured above, which is what makes the kernel load-bearing rather than a
convenience.

**`expect_identical()` is banned package-wide; use `expect_equal()`.** `identical()` is bit-exact on
doubles and also trips on differences that carry no meaning here (integer vs double storage, a
dropped attribute, a name reordering), so a result correct to every digit anyone can act on fails
with a diagnostic that reads like a numeric regression. Chunking makes this concrete: a single pass
sums a column in one traversal, a chunked run sums block by block, and the association order alone
moves the mean in the last bit. Where a bound is the actual claim, state it
(`expect_equal(x, y, tolerance = ...)`). Recorded in CLAUDE.md under "Test altitude".

---

## 2026-07-26 -- chunking is a wrapper that pre-imputes; the seam is a hoist, not an argument

**Decision.** `calc_clocks()` always assumes its `DNAm` is a **complete cohort** -- a precondition,
never a parameter. Chunked scoring becomes a separate front end that pre-imputes each chunk from
cohort-wide means *before* handing it to the shared scoring internal, so no cohort-statistics
argument is added to the public function. Design in `dev/id-streaming-plan.md` (local-only); phases
1-4 keep their existing numbers because four citations in the tracked docs refer to them.

**Why a precondition and not a flag.** Partial-vs-complete missingness is only definable against a
column's NA count over every sample, so `col_miss == nr` in `scan_missing_cpgs()` is wrong the moment
`nr` is a chunk's row count rather than the cohort's. Making completeness a precondition means the
existing machinery is correct by assumption instead of correct by inspection, and no scorer ever has
to ask whether it is seeing everything.

**Why pre-imputation beats an injected `cohort_stats`.** An earlier proposal in the same session was
to pass cohort means plus the column classification into `calc_clocks()`. Pre-imputing outside is
strictly better: it makes the scoring path a pure function of the matrix it was handed, and it is
provably score-equivalent because `observed_panel()` only *orders* cached-before-raw columns and no
scorer branches on fill provenance. A cohort-partial column arrives filled and reads as clean; a
cohort-all-NA column cannot be filled and still routes to vendor-ref-or-drop; and the pathological
case (partial cohort-wide, all-NA within one chunk) is exactly what the fill erases. Rejecting the
argument keeps a chunk-aware parameter off the public one-cohort surface.

**What pre-imputation does *not* fix, which is the real finding.** It erases the record of what was
imputed. `panel_sample_miss()` counts NAs in the raw matrix over cached columns, so with
`partial_cache = NULL` a chunked run would report `score_imputed_partial = 0` for every clock, an
all-zero `sample_miss`, and `check_row_coverage()` would never fire -- breaking "coverage is never
reported for a sample it is not true of" **silently**. The clocks gate survives (`usable_cols` is
cohort-invariant); only the samples gate breaks. So coverage must be built by the wrapper's first
pass from the pre-imputation NA pattern. That is a hoist rather than a duplication -- coverage is
already computed once upstream of the scoring loop -- and it delivers the wanted behavior that
`min_samples_coverage` throws before any scoring.

**Why the fill policy was never negotiable.** Falling back to a vendored constant for partial NA
would have made streaming trivial and was rejected outright: a fixed external constant is a
*different estimand*, not a cheaper estimate of the same one, so it converts a technical missingness
problem into biological error inside a value that then flows into a clock. Two passes is therefore
not a performance preference but the only exact way to honor the existing invariant.

**Chunking and binding are separate problems.** Chunking caps resident set size with **one** set of
fill values; binding is wider and admits records with **deliberately different** fill values (three
time points, imputed per time point). Both yield disjoint ids and identical columns, so no id gate can
distinguish them -- only a recorded per-sample batch label plus a per-batch imputation summary can.
`rbind` therefore records and never refuses on differing fill regimes. A batch label is not
`batch_set_id` returning: it is provenance about origin, not a fingerprint used as a gate.

**Input: an iterator contract with two adapters.** The source must only report dimnames without
loading and yield sample-blocks -- both passes are sequential scans, so there is no random-access
requirement. An **in-memory block iterator** (no dependency) is the always-on test vehicle and the
supported manual path: chunk invariance is a property of the algorithm, not of I/O, so splitting a
matrix in memory proves it as well as reading from disk, and the correctness tier never needs a
Suggests installed. **`DelayedArray`** (Suggests) is the user-facing adapter -- users bring their own
HDF5 and the friction is one line, not an ETL. Bioconductor stays in Suggests only; its twice-yearly
breaking cycle must never be able to archive this package.

**We own the accumulator; DelayedArray only slices.** Its whole used surface is `dim()`, `dimnames()`,
`[` + `as.matrix()`, and `chunkdim()`, because breakage exposure scales with surface area.
`DelayedMatrixStats` is deliberately absent for three reasons that have nothing to do with dependency
hygiene: (1) `rowMeans2()` traverses in its own order and lands ~5.6e-16 from a plain-matrix
`colMeans()`, so one shared accumulator is what keeps the chunked and single-cohort paths
`identical` by construction -- the same failure mode as duckdb's threaded `avg`, from a different
direction; (2) a fused single traversal measured 1.31s against 2.33s for two DelayedMatrixStats
sweeps, and yields the value range and id vector free, where pass 1 needs ~6 quantities; (3)
`setAutoBlockSize()` is a session global a package must not set. `matrixStats` is already in Imports,
so the per-block kernels cost nothing.

**Assume neither orientation nor geometry.** Bioconductor is probes x samples and this package is
samples x probes -- *opposite* -- so orientation is inferred from dimnames, verified against the union
panel, and errors rather than guesses when ambiguous (extending `check_DNAm()`'s existing transposed
warning). Chunk geometry is read from `chunkdim()` and turned into an advisory: measured spread is 12x
on a 100-sample block, `writeHDF5Array`'s auto default is one of the *worse* choices for this access
pattern, and the two passes pull in opposite directions so a square-ish chunk is near-best on both.
Slicing along samples is fixed by fiat -- a choice about our traversal, not an assumption about input.

**npy rejected despite winning every benchmark.** 1.00x storage, 0.018s per 100-sample block (6x
HDF5), 0.1s for pass 1 (19x HDF5), and readable with base R `seek()` + `readBin()` with no package at
all (`RcppCNPy` cannot do partial reads, so it would not have been the route). It loses because npy
carries no dimnames, so it needs sidecars, so *we* would own a format spec plus its validator,
versioning, and edge cases (v1/v2 headers, `fortran_order`, endianness, rejecting float32 whose ~1e-7
precision blows `PARITY_REL_TOL`) permanently. Hand-rolled correctness is a liability a speed win does
not buy off.

**duckdb deferred, not rejected.** It is CRAN, already in Suggests, and the parity fixtures already
use the exact wide schema with a working reader. Cohort validation as SQL is its real strength. It
loses on storage -- measured 1.65x in-memory against HDF5's 0.53x, i.e. 3.1x HDF5, with `DOUBLE`
mandatory so it cannot be tuned down -- and on requiring an ETL, since "duckdb reads whatever you
have" is false and the freeze is the expensive full read+write chunking exists to avoid. If it lands,
the rule is **SQL owns identity and integers, R owns every float**: duckdb's threaded hash aggregate
is not reproducible across thread counts *or* across repeat runs (~1e-15), while `SET threads=1`,
counts, and distincts are. `bigmemory` is rejected outright -- weak dimnames support is disqualifying
when dimnames *are* the identity.

---

## 2026-07-26 -- a fixture census guards the generated parity blocks, behind the parity gate

**Decision.** `test-fixtures-parity.R` gains one census test asserting that every catalog clock
declares a fixture for every `PARITY_COHORTS` cohort, with the sex-routed aliases as the only
clocks declaring none. It is gated on `MC_PARITY` alone -- **not** on a staged cohort, because it
reads only the committed catalog.

**Why it exists.** The four blocks are generated by `parity_targets()` looping over
`clock_fixtures()`. That makes the tier self-maintaining when fixtures are added, but it also means
a fixture upstream *drops* produces no test at all: not a skip, not a failure, just 216 passes
becoming 214 in a run that is still green. Skips are loud -- testthat prints every one with its
reason -- so the skipped case was never the risk. The risk was the test that is never generated.
Before this, the only tripwire was the hand-written standing-state count in `CLAUDE.md`, i.e. a
human remembering a number. This is the same silent-shrinkage failure the routing invariant already
forbids on the other axis ("a sync that adds a routing pair no branch claims must fail the
always-on tier, never silently shrink it"); the fixture axis had no equivalent.

**Why a rule and not a snapshot.** The declaration surface is exactly regular: 129 catalog clocks,
122 declaring both cohorts (244 targets = core 130 + fitage 28 + packs 56 + horvath 30), and the 7
with none are `setequal` to the sex-routed aliases. So the census states a rule -- "both cohorts, or
you are an alias" -- rather than a committed list of ids or per-block counts. A count ledger would
catch the same drops but churn on every legitimate upstream addition, which trains people to bump
it without reading it. Both halves of the rule derive from existing declarations
(`clock_fixtures()`, `sex_routed_members()`), so a new clock or a new sex-routed family needs no
edit here. The residual hole is a clock vanishing from `mc_catalog` entirely -- the census iterates
the catalog, so it cannot see that; closing it needs a real id ledger, and a missing score column
is a louder signal anyway.

**Why gated, given it costs nothing to run.** Ungated it would run on every `devtools::test()` and
on CRAN, which is precisely when the tier it guards is skipped -- the argument for always-on. It
lost on scope: a fixture census asserts an upstream contract's completeness, which is maintainer
and CI concern in the same way parity itself is, and CRAN can neither see the meta repo nor act on
a dropped fixture. CI stages the cohorts and sets the flag, so the gate is open where a bad sync
actually lands. Accepted consequence: a plain local `devtools::test()` will not catch a dropped
fixture.

---

## 2026-07-25 -- Zhang2019 parity lands by fixing the fixture's *input*, not its tolerance

**Decision.** `Zhang2019` leaves `KNOWN_PARITY_GAPS`, which is now empty. The parity tier gained
`needs_full_panel()` + `cohort_betas_full()`: a target whose worklist contains a clock declaring a
`sample_scale` recipe op is scored on the cohort's **whole array**, not on the union of scoring
panels. Both cohorts pass at the unchanged `core` bounds -- `max_abs_diff` 9.8e-13 @ EPICv1 and
3.5e-12 @ 450K, `max_rel_diff` 6.2e-14 and 8.5e-14.

**Why the old skip reason was wrong.** It read "sample_scale moments over needed-CpG subset, not
full panel; exact parity unreachable" -- but "not full panel" described the *test harness*, not the
clock. `score_Zhang2019()` already takes its moments over whatever `calc_clocks()` is handed, and
real users hand it their array; only the fixture loader was narrowing the input to the 514-CpG
scoring panel. So the unreachable-ness was self-inflicted. Measured on the subset the miss is 1.8e1
absolute / 82% relative -- which is exactly the size you get from shifting a z-score's centre, and
is why no tolerance was ever the answer here.

**Why the rule is read off the recipe and not a `FULL_PANEL_CLOCKS` list.** The 2026-07-23 entry
below proposed a `FULL_PANEL_CLOCKS` constant; it was forward-looking text and never landed (that
entry's "Zhang leaves `KNOWN_PARITY_GAPS`" paragraph described work that stopped at the scorer).
A hand-kept list is the thing every other block in this file was refactored away from: the fact
"this clock's moments span the input matrix" is already declared upstream as `recipe[[i]]$op ==
"sample_scale"`, so the predicate reads that and a future sample-scale clock gets the full array
without anyone remembering to add it.

**Cost.** The full array is 866k x 71 (EPICv1) / 473k x 80 (450K) -- ~550 MB and ~350 MB, about 1s
to read and 1s to score, and only for the two Zhang targets. That is affordable in a tier that is
already cohort-gated and CRAN-skipped, and it is the only honest input.

---

## 2026-07-25 -- the horvath_online gap is the oracle's input, not our arithmetic; skip it

**Decision.** A fourth parity block, `horvath`, selected by the declared `fixtures[].oracle ==
"horvath_online"`, skipped wholesale. Parity is now 214 pass / 32 skip / 0 fail. This entry also
**corrects two claims made earlier the same day** (see "Corrections" below) -- both were wrong, and
the earlier entry is left standing per append-only.

**What the oracle actually did.** The staged submission (`fixtures/_horvath_submit/{cohort}/
server_matrix.csv`, 32,136 probes) has **2,567 (EPICv1) / 3,208 (450K) rows that are entirely NA
and zero partial NAs**. The server filled those server-side with a per-probe constant, then -- for
the `DNAmAge` column only -- BMIQ'd the 21k panel. Neither the fill values nor the BMIQ'd matrix
are published. We hold the raw betas, so we drop what it filled.

**Why this is provably not our bug, which is why a tolerance would be a lie.** Every one of these
clocks is a matmul plus an intercept, so only two things can be wrong: the tensors/intercept, or
the implementation. Score against the *submitted* matrix and the residual sorts perfectly by
absent-probe count and by nothing else:

| cohort | clock | absent probes | max rel |
|---|---|---|---|
| 450K | Hannum | 0 | 6.1e-08 |
| 450K | DNAmCystatinC_GrimAgeV1 | 0 | 8.7e-08 |
| EPICv1 | DNAmTL | 0 | 4.2e-08 |
| EPICv1 | DNAmTIMP1 | 0 | 5.8e-08 |
| EPICv1 | Horvath2 | 0 | 2.5e-07 |
| EPICv1 | GrimAgeV1 | 1 | 7.9e-03 |
| EPICv1 | DNAmPAI1 | 1 | 1.1e-01 |
| 450K | DNAmLeptin | 32 | 5.1e-02 |
| 450K | GrimAgeV1 | 72 | 6.4e-02 |

Agreement to ~1e-8 relative on five independent clocks across both cohorts is not something a wrong
coefficient or a wrong matmul survives. The tensors and the engine are right; the inputs differ.
A relaxed bound would have to span 4.2e-08 to 2.7e-01 -- six orders of magnitude -- so it could
only ever be vacuous. A skip states the truth; a tolerance would launder it.

**What we ruled out, so nobody re-runs it.** Three candidate fill tables all fail: datMini4
`ProbesMedian` (helps some, badly hurts others -- `DNAmLeptin @ 450K` goes 5.1e-2 -> 1.0), the 27k
`overallMeanByCpGacross50data` (barely differs from no fill), and the 21k `goldstandard2` (covers
only the 21k set; `GrimAgeV1` has 47 of 1030 probes there). `betanorm::bmiq_calibration(niter = 5,
nfit = 21368, goldstandard.beta = probeAnnotation21kdatMethUsed$goldstandard2)` **does** reproduce
the BMIQ step -- `Horvath1 @ cohort_450K` improves 2.7e-1 -> 1.9e-3, a 140x gain, and it is the
only clock BMIQ helps -- but EPICv1 stays bad because its 1062 filled probes contaminate the fit.
Upstream's reference script (`aldringsvitenskap/epigeneticclock/StepwiseAnalysis.R`) confirms the
order (impute, then BMIQ, then score `DNAmAge` off the normalized matrix) and names
`goldstandard2` as the fallback fill -- but only above 3000 missing, and only for the 21k
`DNAmAge` path, so it does not cover the surrogates.

**Corrections to the entry below, made the same day.**
1. It claimed the `server_normalization` tier "asserted the oracle saw different betas for 28 pairs
   where it demonstrably did not." Wrong: the oracle *did* see different betas. The field was
   naming the wrong mechanism (normalization) for a real divergence (the fill), and `["none"]`
   is a truthful statement about normalization for those 14 clocks.
2. Mid-investigation we concluded from the run's `Comment` column ("Beta mixture quantile
   normalization was done", present on 71/71 and 80/80 rows) that BMIQ hit every column and that
   upstream's `["none"]` was a metadata bug on 14 clocks. **Also wrong** -- BMIQ was *run* on the
   submission, but only `DNAmAge` consumes the normalized matrix. Measurement settled it: BMIQ
   improves `Horvath1` 140x and makes all 14 others worse by orders of magnitude. Upstream's
   metadata is correct as declared. No upstream bug was filed.

**The skip reads the declared oracle, not a clock list**, so a regenerated fixture retires it
automatically -- the property the `DNAmFitAge` group skip lacked, which is why that one hid 28
*passing* tests for an unknown number of syncs.

---

## 2026-07-25 (later) -- parity: two axes at 1e-10, blocks split on group, correlation banned outright

**Decision.** Reverses the *tier* half of the earlier entry below, and tightens the grade. Parity
now grades **both** `max_abs_diff < 1e-10` and `max_rel_diff < 1e-10`; the blocks are **non-FitAge**
and **DNAmFitAge**, split on `clock_group_id()`; the `DNAmFitAge` group skip is gone; and
"no correlation as a gate" is promoted from a parity-local remark to a package-wide invariant in
`CLAUDE.md`.

**1. Correlation is banned as a rule, not as a local choice.** The earlier entry retired it from
parity only. That left it available to the next person writing a unit test, which is how it got in
the first time. The invariant now names the two failure modes explicitly -- it is blind to
intercept shifts and to rescaling, and it concentrates near 1 so a few catastrophically wrong
samples vanish into a cloud of correct ones -- and extends the same reasoning to `median`/`mean` as
reducers, which is what `rel_diff()` was still doing. A median relative difference passes with half
the cohort arbitrarily wrong. `rel_diff()` is now a `max` and a **gate**, not a printed diagnostic.

**2. One axis was never enough, in either direction.** Absolute-only at `1e-6` was failing clocks
that are numerically correct: measured on the current cohorts, `DNAmCystatinC_GrimAgeV1` sits at
~4e-8 relative on both cohorts and `DNAmLeptin @ cohort_EPICv1` at ~8e-7, and both blew the
absolute bound purely because their units are large (pg/mL-scale biomarkers, not ages). Relative-only
has the mirror hole: a score near zero admits unbounded absolute drift. Requiring both, at the same
harsh `1e-10`, is the only form that is scale-free without being blind. This supersedes the
`detail-plan.md` sec 13 follow-up proposing numpy-allclose (`atol + rtol*|exp|`): that form is a
**sum**, so a generous atol silently forgives a relative miss. Two independent bounds do not.

**3. The `server_normalization` tier split was measuring the wrong thing.** It classified a fixture
as "horvath online" when `server_normalization` was non-empty. But the field is a registry name
list, and **28 of the 30 pairs in that tier declare `["none"]`** -- the server applied no
normalization, so those oracles saw exactly the raw betas we score. Only `Horvath1` (both cohorts)
is genuinely `["BMIQ"]`. The tier label therefore asserted "the oracle saw different betas" for 28
pairs where it demonstrably did not, which is worse than no label: it excused real divergences as
expected. "Non-empty" was in practice a proxy for `oracle == horvath_online` (the field is absent
for every other oracle), which is a different question from the one the tier claimed to answer.

**4. The blocks split on group id instead.** `DNAmFitAge` is one block, everything else is the
other. Membership derives from `clock_group_id()`, so there is still no hand-kept clock list.

**5. Unskipping FitAge found 28 passing tests, not 28 failing ones.** This is the entry's most
useful result and it inverts the previous one. The group skip was added on the belief that the
fixtures came from a truncated Horvath submission; upstream has since regenerated them from author
code (every member now declares `oracle: author_code`, `server_normalization` absent), and the
whole family clears **1e-10 on both axes, 28/28**. The skip had outlived its cause by an unknown
number of syncs and nothing would have noticed -- a skip keyed to an upstream defect has no way to
tell you the defect is fixed. Standing state is now FitAge green, non-FitAge 178/36/2. FitAge keeps
its own block anyway: the retired skip should stay visible until the family has a run of clean
history. `KNOWN_PARITY_GAP_GROUPS` is now empty but **kept as a separate map** -- the namespace
collision that motivated it (`DNAmFitAge` is both a group id and a clock id, as are `EpiTOC2`,
`prcPhenoAge`, `RepliTali`, `SystemsAge`) is unchanged, so folding it into `KNOWN_PARITY_GAPS`
would reintroduce the ambiguity the moment a group needs skipping again.

**6. Tightening to 1e-10 cost 7 new failures, and they were the predicted ones.** `DNAmTL @
cohort_EPICv1` (the old lone tier-2 pass) plus `PCB2M`, `PCCystatinC`, `PCPAI1` on both cohorts.
The PC trio is exactly the class `detail-plan.md` sec 13 flagged as scale-sensitive, so their
appearance at a harsher bound was the instrument working, not a regression -- and it is what
motivated the third block (point 7).

**7. A third block for the external packs, relaxing only the scale-sensitive axis.** `parity_block()`
now classifies each clock as `packs` / `fitage` / `core`, and `PARITY_ABS_TOL` became
`c(core = 1e-10, fitage = 1e-10, packs = 1e-6)`. `PARITY_REL_TOL` stayed `1e-10` for every block
and **must not** become per-block: it is scale-free, so a units argument can never justify moving
it, and per-block relative slack is how a tolerance table turns back into the permission slip
point 3 describes. Splitting the two axes is what makes "relax this block" a bounded statement
about units rather than a blanket loosening.

The pack number is measured, not chosen for convenience. Over all 56 pack rows the worst absolute
miss is 2.05e-8 (`PCB2M @ cohort_450K`, magnitude 3.26e6) -- 6.3e-15 relative, i.e. below the float
representation floor for a value that size before any arithmetic runs. The worst *relative* miss
anywhere in the block is 7.5e-13, 130x inside the untouched relative bound, so the block is graded
strictly on the axis that means something and the relaxation costs no real sensitivity. `1e-6`
clears the observed worst by ~50x and is still ~3e-13 relative at `PCB2M`'s scale. All 56 pack rows
pass; nothing was skipped to get there.

**8. What is left red is now a single, named class.** After the split: core 128 pass / 30 fail / 2
skip, fitage 28/28, packs 56/56, PhysAge 2/2. The 30 failures are exactly the 15 clocks whose
fixtures declare `oracle: horvath_online`, on both cohorts, and nothing else. That is a far more
useful queue than the 36 mixed failures it replaced -- the pack noise is gone, FitAge is gone, and
what remains points at one oracle.

**Not done here.** Root-causing the 36. The measurement of which are scale artifacts vs. real
divergence, and the cohort asymmetry the earlier entry noted, is deliberately left for the
follow-up; this entry only fixes the instrument. One datum worth keeping: our derived absent-CpG
set was checked against upstream's declared `fixtures[].missing` sidecar for all 30 formerly-online
pairs and agreed exactly, so panel resolution is not implicated. That sidecar is still unread by
the test suite.

---

## 2026-07-25 -- the citation verb is `cite_clocks()`, a generic we own

**Decision.** The `citation()` plain function the plans promised is now `cite_clocks()`, an S3
generic with `character` / `mc_result` / `default` methods, returning an `mc_citation` record.
Shape and rationale: detail-plan sec 8.1.

**Both obvious names are taken, and neither is a generic.** `utils::citation` was ruled out on
2026-07-23. `cite` looks free and is not: `utils::cite(keys, bib, ...)` is the Rd `\cite` /
`bibstyle` helper, also a plain function, so `cite <- function(x, ...) UseMethod("cite")` masks it
for every user who attaches us -- the identical trap one name over. Check `find()` before minting a
generic; this is the second time the reflex name was already occupied. `cite_clocks()` is
unclaimed and joins the family the package already has: `list_clocks()`, `calc_clocks()`,
`cite_clocks()`.

**Dispatch earns its keep at the argument, not the output.** The character method routes through
`resolve_clocks()`, so the token vocabulary (ids, groups, tags, `"all"`) and the routed-member
refusal are the ones `calc_clocks()` already teaches -- one namespace, not two. The `mc_result`
method cites `colnames($scores)` and **does not walk `clock_inputs`**: if a composite owes its
components' papers, upstream says so with a `cite_also` row. Deriving that edge from the dependency
graph would be us manufacturing a scientific claim, and it is the same "never search, read the
declaration" line the accessors hold.

**Raw BibTeX, not `bibentry`.** Base R ships `toBibtex` (a writer, and a real generic we may method
onto) and no reader. Real `bibentry` objects would cost either a dependency or a LaTeX-aware parser
of our own -- for a print style. Slicing declared entries out of the vendored file is
dependency-free and byte-faithful, which matters more than it sounds: 12 entries carry UTF-8
author names (accented Polish/German/Spanish surnames), so the export line is
`writeLines(toBibtex(x), "refs.bib", useBytes = TRUE)`. Accepted cost: no `style = "text"`.

**No new export verb.** `print` / `as.data.frame` / `toBibtex` cover console, table and file. A
`write_bib()` or a `file =` argument would be a third spelling of `writeLines`.

**Registration is the maintainer's.** `#' @export` tags are in the source and nothing else was
touched; the NAMESPACE entries this implies are listed at the handoff. Until they exist, dispatch
works under `testthat` (tests run in the package namespace, where `UseMethod` finds the methods
lexically) and **fails from a user session**, which is a green tier that proves nothing about
reachability -- noted in CLAUDE.md so the next person does not read the tier as evidence.

---

## 2026-07-25 -- citations key on the clock -> paper join, not on `meta.pmid`

**Decision.** `sync.R` reads `bibliography/clock_citations.csv` and ships the whole 1:N join as
`mc_citations`. `read_papers_csv()` and the per-clock `entry$bib_key` are deleted; `mc_index$bib_key`
is replaced by `n_citations`. Design and the resulting invariants: detail-plan sec 8.1.

**The bug was invisible and would have stayed that way.** We joined a clock's *scalar* `meta.pmid`
to `papers.csv` to get one `bib_key`. Upstream made clock -> paper 1:N on 2026-07-25: a clock can
cite N papers with `role` `primary` or `cite_also`, and `meta.pmid` only ever names the primary one.
Every one of the 124 rows on master today is `primary` and equals its meta's `pmid`, so the old path
returns the right answer for every clock that currently exists. It would have started silently
under-citing the day a human added one `cite_also` row -- no error on either side of the boundary.
This is why the fix landed with no failing test to motivate it, and why "no test caught it" is not
evidence it was fine.

**The `cite_also` gate stays upstream, and that is the whole of it.** No assertion over the real
table distinguishes a correct reader from a broken one while every row is `primary`, so testing the
1:N contract needs a synthetic fixture -- which upstream already has in `tests/test_gate.py`. We do
not repeat it: `sync.R` is maintainer-side tooling the package neither ships nor depends on, and
reaching its functions from `tests/` means sourcing an Rbuildignored file. A downstream version was
written and removed the same day; `tests/` covers `R/` and the shipped tables only (rule added to
CLAUDE.md "Test altitude"). What we do assert downstream is the *shipped* result: every catalog
clock reaches >= 1 citation, exactly one `primary` per cited clock, aliases resolve through their
donor, and every `bib_key` exists in `clocks.bib`.

**One store, and the alias points rather than copies.** Keys live once, in `mc_citations`. We
rejected a `bib_keys` list-column on `mc_index`: it is the same join in a second place, which is the
duplication upstream keeps deleting. A scalar `bib_key` was worse than absent -- it cannot hold N
keys, so it silently reports the first. Sex-routed aliases are package-minted and have no upstream
row; they carry `donor_clock_id` and cite through the donor, so the alias holds no copied key that
could survive as the last instance of the bug. `pmid` **stays** on the entry and the index: it is
honest primary-paper provenance, just not the citation source. We did not rename it `primary_pmid`
-- `FIELD_REGISTRY` mirrors upstream meta field names, and diverging there costs more than the
ambiguity a comment resolves.

**`clocks.bib` did not move.** Structured join in `sysdata`, human-facing document in `inst/`: users
want a real `.bib` on disk to feed a manuscript. Parsing it into `sysdata` would give us a second
copy of the citation text.

**Also in this sync (7db9f2c), not a package change.** `RetroAge450K` / `RetroAgeEPICv2` ->
`Retroelement_AgeV1` / `Retroelement_AgeV2` (group `RetroAge` -> `Retroelement_Age`): the paper
describes no 450K clock, so the old id asserted an array split the authors never made. Everything is
derived, so `sync.R` needed no change and the panels are bit-identical. We added **no** alias from
the old ids -- upstream's rule is that a consumer pinned to something wrong should break -- and
nothing in `R/` or `tests/` was pinned to them. SystemsAge's primary paper also moved from the 2024
preprint (37503069) to the 2025 publication (40954326) across 13 clocks. The three external packs
rebuilt to **identical** `payload_hash`es, so no release upload is owed.

---

## 2026-07-25 -- parity: drop correlation, split the horvath online tier, skip the FitAge family

**Decision.** Every parity fixture is now graded **exact** (`max_abs_diff < 1e-6`). The
`server_normalization` field no longer selects a *tolerance*; it selects a *tier*. And the
`DNAmFitAge` family is skipped wholesale pending regenerated upstream fixtures.

**1. Correlation was the wrong gate, and loosening it further was the wrong fix.** Server-normalized
fixtures were graded `cor > 0.99`, which is offset- and scale-invariant: it cannot distinguish "we
match the oracle" from "we are uniformly wrong by a constant." An in-flight patch had already
noticed this and bolted a median-relative-difference bound (`< 0.10`) alongside the correlation.
That was still the wrong shape -- a single group-wide relative tolerance has to be set wide enough
for the *worst* clock in the group, and the measured worst (`Horvath1 @ cohort_EPICv1`, median rel
8.2e-2) is seven orders of magnitude looser than what most of the group actually achieves. A
`Hannum` regression from 1e-8 to 1e-2 would have sailed through. Both the correlation bound and the
relative bound are gone; `rel_diff()` survives only as a number printed in the failure label, to
make root-causing faster.

**2. We deliberately let the tier go red.** 29 of the 30 non-FitAge (clock, cohort) pairs now fail
(`DNAmTL @ cohort_EPICv1` is the lone pass). That is the point: the previous policy's job was to
keep the suite green over a divergence nobody had measured, which is how the divergence survived.
Standing failures, one test per pair, are the work queue. Tier 1 (raw-beta oracles) is unaffected
and fully green -- 187 pass, `Zhang2019` skipped -- which is itself evidence the exact tolerance is
achievable and the problem is specific to the online oracle.

**3. What the measurement actually says (and why the tier split is worth having).** The divergence
sorts by *cohort*, not by clock. On `cohort_EPICv1` most of the group is exact to ~1e-8 and only
`Horvath1`, `DNAmPAI1`, `GrimAgeV1` and `GrimAgeV2` move; on `cohort_450K` nearly everything sits
at 1-6% median relative difference. BMIQ cannot explain that -- server normalization was applied to
both cohorts, so it would not spare EPICv1. Absent-probe drop-vs-vendor-fill on the sparser 450K
panel can. Naming the tier "horvath online" rather than "graded loosely" is what makes that
question askable: a tier-2 failure means "the oracle saw different betas", a tier-1 failure means
"our arithmetic is wrong", and they route to different people.

**4. The split derives from the declared field, not a clock list.** `parity_targets(online =)`
filters on `server_normalization` being non-empty and feeds two loops over one shared
`run_parity_target()` body. A hand-kept "clocks that cannot be exact" table was rejected for the
same reason it was rejected the first time: it rots the moment an oracle changes, and it turns a
measurement into a permission slip.

**5. `DNAmFitAge` is skipped as a group, in its own map.** The upstream fixtures were generated by
a Horvath submission that omits the CpG carrying the largest coefficient, which blows the error up;
they need regenerating from the author code upstream before the numbers mean anything. Skipping
by group needed a second map (`KNOWN_PARITY_GAP_GROUPS`) rather than a key in `KNOWN_PARITY_GAPS`,
because group ids and clock ids share a namespace -- `DNAmFitAge` is *both* a group and a clock
(as are `EpiTOC2`, `prcPhenoAge`, `RepliTali`, `SystemsAge`), so a single flat map could not say
which one a key meant.

**Kept.** One file, one duckdb connection per staged cohort. The tier split is two loops inside
`test-fixtures-parity.R`, not a second file -- splitting files would have forced the per-cohort
connection setup into a `helper-*.R` for no behavioral gain.

---

## 2026-07-24 -- one noun for the external payloads ("assets"), a session-wide assets dir, and stale reclaim

**Decision.** The external-pack surface said "cache" and "assets" and "data" for overlapping
things, and one argument meant both a place and a thing. Standardized on **assets** in public
names, split the dir getter from the resolver, added a setter, and taught the clearer to reclaim
superseded packs.

**1. "Assets" is the thing; "cache" is only the storage policy.** In SE terms an asset (GitHub
release asset) or artifact is the payload; a cache is a location holding derived, refetchable
copies. Both were correct, which is exactly why using them interchangeably was confusing -- and
the package had a *second*, unrelated `partial_cache` (the cohort-mean fill, `R/missingness.R`)
threaded through every scorer. Since `partial_cache` is internal-only and never appears in an
exported name, taking "cache" out of the **public** names removes the collision without touching
the scoring code. `mc_data_*` was rejected outright: `R_user_dir(which = "data")` means *user data
you must not delete*, so "data" names the wrong lifecycle.

The public set is now, and is the whole pathing surface:

| `get_mc_assets_dir()` | getter -- the dir in effect |
| `set_mc_assets_dir(path = NULL)` | setter -- `NULL` clears; returns the old value invisibly |
| `list_mc_assets(groups)` | read-only table: size / `downloaded` / `superseded` per group |
| `download_mc_assets(groups, ask)` | fetch kernel: bytes -> disk |
| `load_mc_assets(groups, from, ask)` | reader: fetch-if-missing -> RAM |
| `clear_mc_assets(groups, ask)` | delete every pack we put there, superseded included |

Naming rule these follow: **every public name is `<verb>_mc_<noun>`** -- already the settled
convention (`load_mc_assets` was renamed to the verb form on 2026-07-21), with `mc_data_download()`
the lone outlier, now `download_mc_assets()`. That the maintainer reached for `download_mc_assets`
and `set_mc_path` unprompted, and could not find either, is the evidence the convention was right
and the code broke it in one place. The getter was briefly the bare noun `mc_assets_dir()` on the
theory that getters are nouns; that left one non-verb in a five-function set and an asymmetric
get/set pair, so it is `get_mc_assets_dir()` -- matching `getwd()`/`setwd()` and
`usethis::proj_get()`/`proj_set()`. No bare-noun accessor survives.

**2. The place/thing conflation lived in the `assets` argument.** It accepted a directory path
*or* loaded packs, which is the confusion made executable. Split: `assets` is now always the
things, and the argument naming the source is **`from`** (`NULL` open set, a path = closed set, or
packs in memory) on `load_mc_assets()` / `calc_clocks()` / `sim_DNAm()`.

**3. Dropped the per-call dir override from `download_mc_assets()` and `clear_mc_assets()`.**
There it could only ever be a plain path -- no closed-set meaning -- so once a session setter
exists it is pure duplication of the option layer, and it is the precise footgun behind the
2026-07-23 incident (`clear_mc_cache(assets = my_pack)` resolved a pack to a garbage dir and
deleted all three cached packs). Tests that pointed at a tempdir with `assets =` now use
`withr::local_options(mc.assets_dir = ...)`, i.e. the same mechanism users get. Option and env var
renamed to `mc.assets_dir` / `MC_ASSETS_DIR`; precedence is unchanged (`from` > option > env >
default).

**4. Re-introducing a setter reverses the 2026-07-22 deletion of `mc_set_cache_dir()`, and the
reason is new.** It was deleted as unreferenced, which it was. What changed is the HPC case: the
`R_user_dir()` default is wrong there (home quota, not visible from compute nodes), so users need
to point at scratch. `options()` alone accepts anything and a typo'd path fails inside `mc_fetch()`
*after* the consent prompt; `set_mc_assets_dir()` asserts a single non-empty string, expands,
creates the dir and probes writability, so a bad path fails on line 1 of the job script. It
returns the previous option value invisibly so `on.exit()` / `withr::defer()` can restore it, and
`NULL` **clears** the option rather than pinning the default -- otherwise it would permanently
shadow the env layer and could not round-trip.

**5. The dir must stay `tools::R_user_dir(which = "cache")`.** With "cache" gone from the public
names, the property CRAN cares about -- derived, reclaimable, safe to delete -- is no longer
self-evident from the API. Writing it down here so nobody later "fixes" it to `which = "data"`,
which would be a policy violation. `clear_mc_assets()` remains a real delete for the same reason.

**6. `clear_mc_assets()` removes superseded packs unconditionally -- clear means clear.** Pack
filenames are content-addressed (`<stem>-<sha256>.qs2`), so every `sync()` that moves a
`payload_hash` leaves the old file orphaned: the clearer only ever deleted files it could name from
the current manifest, so the dir grew without bound and the space was unreclaimable -- the exact
thing the CRAN requirement exists to prevent. `mc_stale_files()` now runs on every call and its
results always join the deletion set, labelled `<group> (superseded)`.

**This settled on the third try, and the middle one was wrong.** A first draft defaulted
`stale = TRUE`; that was walked back to `remove_stale = FALSE` plus an informational nudge, on the
theory that a plain `clear` should not delete files the caller never named. The flag is now gone
entirely, because that theory misplaced the safety boundary: **the consent gate is the safety, not
the default.** Nothing is ever deleted without the user seeing the full file list and saying yes
(or passing `ask = FALSE`), so a "safe" default bought nothing and cost the verb its meaning -- a
`clear` that leaves files behind is the actual surprise, and it left the CRAN reclaim story needing
a flag most users would never find. `pak::cache_clean()`, `pip cache purge` and `npm cache clean`
all wipe; that is the settled ecosystem meaning of the word. The residual worry -- someone's own
file caught by the scan -- requires a file inside a package-owned `R_user_dir()` cache named with
our declared stem plus 64 hex plus `.qs2`, which is opting in, not stumbling in.

`mc_delete_summary()` keeps the two kinds visible at the decision point: the prompt, the
non-interactive refusal and the completion message all read "3 downloaded packs and 3 superseded
packs", each part formatted separately so its `{?s}` binds to its own count.

This is a deliberate, narrow carve-out from "accessors never search." The scan is not payload
resolution: it never returns a pack to a scorer, it only enumerates files to delete. It stays
anchored to declarations anyway -- the stem is read off the **declared** `file` field
(`mc_asset_stem()`), never guessed from the group id, and only the 64-hex hash is a wildcard -- so
a foreign stem, an uncontent-addressed file, or a group whose declared name does not match the
expected shape are all skipped rather than guessed at. A stale-*only* mode is deliberately not
exposed: minimal surface pre-alpha, and promoting later is cheap where un-exporting after release
is not.

**7. `list_mc_assets()` exists because every read-only question needed a mutating verb to answer
it.** This reverses the same-day call to keep `mc_cached_files()` internal ("minimal surface,
promote later"), and what changed is stale detection. With three distinct on-disk states
(declared-and-staged, declared-and-missing, superseded) the only way to see any of them was
`clear_mc_assets()` -- a *deleter*, which refuses non-interactively -- and the only way to preview
download sizes was to trigger the consent prompt or actually fetch. Read-only questions get a
read-only answer: a `list_clocks()`-style data.frame of `group_id / n_clocks / n_cpgs / size /
downloaded / superseded / superseded_size`, no prompting, no network, no writes.

**Not** an overload of `download_mc_assets()` (the considered alternative). Its no-arg form already
means "download everything", so re-purposing it silently breaks that; the return type would vary by
arity; and a verb named `download` that sometimes does not download is the same place/thing
polymorphism this entry deleted from `assets`. Sizes are `fs::fs_bytes` columns -- they print as
`8.88M` and still sum -- assigned after `data.frame()` so the class is not stripped. Maintainer-side
plumbing (`release_tag`, `file`, URLs) stays out; it does not inform any user decision. One
follow-on cleanup: `mc_stale_files()` now keys by plain `group_id` and the `" (superseded)"`
decoration is applied at the display layer (`mc_stale_labels()`), so the lister does not have to
re-parse a label to recover the group.

**8. Consent prompts: cli renders, `askYesNo()` asks one line.** The prompts were built with
`paste0()` into a single multi-line string handed to `askYesNo()`, which forwards it to
`readline()` -- whose `prompt` argument is meant to be one short line. Embedded newlines render
malformed on Windows, which is how this surfaced. All prompts now go through `mc_ask_yes_no()`:
`cli_inform()` for the header and dir, `cli::cli_verbatim()` for the aligned size table (verbatim
is the piece that matters -- cli would otherwise reflow the padding away), then a single-line
question. Its `header` arrives pre-formatted via `cli::format_inline()` in the *caller's* frame,
because a plural marker or interpolation evaluated inside the helper would resolve against the
helper's variables, and interpolating the finished string as a value stops it being re-parsed for
braces.

The same pass fixed a latent instance of the documented reflow hazard: the **non-interactive**
refusals interpolated the aligned block into a `cli_abort()` bullet (`" " = "{manifest}"`), where
cli collapsed the whole table onto one line. Conditions cannot use `cli_verbatim`, so those paths
now carry no alignment at all -- `mc_manifest_bullets()` emits one self-contained `group (size)`
bullet per pack. Rule of thumb now in CLAUDE.md: aligned block -> `cli_verbatim`; anywhere reflow
is unavoidable -> one self-contained bullet per row.

**Also.** The five functions are exported for the first time (hand-written `NAMESPACE`, roxygen
still off) -- none of them were reachable without `load_all()`, which is the other half of why
they looked missing. `calc_clocks()` and the rest of the API remain unexported; that is a separate
gap.

---

## 2026-07-24 -- always-on tier trimmed to one golden per scoring path; sim-smoke kept for its configuration, not its clock list

**Decision.** Cut ~114 lines of duplicated always-on tests, and restate why `test-sim-smoke.R`
survives a coverage argument that looks like it should delete it.

**1. One golden per scoring path, not per clock.** 122 of 129 catalog clocks now have parity
fixtures on both cohorts, and the tier passes clean. A **fully staged** run -- both cohorts plus
all three external packs in the runtime cache -- is **492 passed / 0 failed / 0 errors / 0 warnings
/ 2 skipped**, the 2 being the declared Zhang2019 gaps. With the packs absent it is 380 passed /
58 skipped, the extra 56 skips being the external pairs; that partial shape looks healthy and is
easy to mistake for a complete run. The always-on tier had drifted into
re-deriving the same path on several clocks: the vendor-mean fill offset five times, drop policy
twice, linear-mean reduction twice, and the generic `used = present + imputed_full` identity three
times. Each duplicate is a separate breakage on any engine change, for no extra coverage. Kept one
golden per path; the identity now lives only in `test-coverage-report.R`. Also folded
`test-accessors-external.R` into `test-score-external.R` (both defined a `fake_pcclocks_pack`),
dropped three coverage-gate tests that `test-score-dunedin.R` restated with Dunedin, and unpinned
the per-clock fuzzy-suggestion ranking, which breaks on any change of string metric with no
behavior change.

Two in-source comments were justifying dead work and were fixed: `test-score-dunedin.R` claimed
DunedinPACE was "parity skip-listed" (it has two passing fixtures and is not in
`KNOWN_PARITY_GAPS`), and `test-score-miage.R` said "parity owns the numeric golden" directly above
three re-derivations of it.

**What was NOT cut, and why.** The FitAge "vendor-fill from that sex's medians" golden was removed
and then restored: the coverage counts alone prove five CpGs were filled but not *whose* reference
filled them, and selecting the sibling sex's medians is exactly the bug that test exists to catch.
External pack builders stay whole. Parity covers external clocks only when a maintainer has
downloaded the packs into `R_user_dir()` -- machine state, not repo state -- so a fresh clone, CI
without a download step, and CRAN all skip those 56 pairs. Independently of caching, parity never
exercises what those tests assert: the in-memory `assets =` path (it loads from the disk cache),
the closed-set refusal on an absent or wrong pack, PCClocks batching all members in one call, a
subset request returning no extra columns, or SystemsAge whole-group orchestration -- parity
requests one clock at a time. Bundled goldens stay because CRAN
and a default `devtools::test()` never run parity; deleting them on the grounds that "parity owns
the numbers" would ship a package whose own test suite has no numeric gate.

**2. `test-sim-smoke.R` is kept for its configuration, not its clock coverage.** On clock identity
parity dominates it: parity invokes 86 of the 87 clocks smoke does (all but Zhang2019) and adds 29
more when packs are cached, and a `calc_clocks()` blowup surfaces there as an error just as
loudly. So "catches a clock that stopped running" -- the justification the file and `CLAUDE.md`
both carried -- is not a reason to keep it, and someone reading the overlap will try to delete it.
The real reasons: parity is skipped on CRAN and without a staged 670MB cohort, so it contributes
zero to a default run; it scores with `min_clocks_coverage = 0, min_samples_coverage = 0`, so it
structurally cannot exercise the default gate path; and it never calls `sim_DNAm()`. Smoke is 29
lines and 1.0s against a cloned meta repo and ~20 minutes.

Net: always-on tier 1564 -> 1450 lines, 441 -> 414 assertions, parity unchanged. Wall clock was
not the motivation and is not claimed as a result -- the tier was already ~10s, and run-to-run
noise on this machine (7-13s for the same suite) swamps the difference. The win is duplication
surface: one engine change no longer breaks five tests.

---

## 2026-07-24 -- coverage records are written once; accessors read declarations through four helpers

**Decision.** Two cleanups with no behavior change, both verified against the sex-routed path.

**1. `compute_coverage()` builds records last.** It used to build every non-alias record in pass 1
from the *unmasked* member counts, then reach back in pass 3 and overwrite
`per_clock[[id]][["score_imputed_partial"]]` with the masked sum. Correct, but a write-then-patch:
the record was briefly wrong, one field was special-cased outside `coverage_record()`, and any
future field derived from `score_miss` would have silently missed the patch. The order is now
**raw counts -> stitch aliases -> mask members -> build records**, so `coverage_record()` is the
only writer of record fields and each is written once. The two-phase shape is forced and stays
documented at the top of the file: aliases must stitch from counts *before* masking (a row's count
is its member's raw count), records need them *after*.

`coverage_for()` is gone -- it and `coverage_record()` were two functions with one caller between
them, in two files, both computing `clock_impute(...)[["policy"]]`. The policy arithmetic moved
into `coverage_record()`, which `CLAUDE.md` already names as the owner of record fields.
`coverage_record()`'s dead `policy =` default went with it, and its partial sums are now
`na.rm = TRUE` unconditionally, since NA in a miss vector means "this sample was not scored by
this clock".

**2. `construct_mc_result()` reads `normalizes` instead of re-deriving it.** It picked the norm
columns by testing `!is.null(sample_miss$norm[[id]])`. Same answer -- both trace to
`length(nm$needed) > 0L` -- but it was the only reader inferring the declared panel fact rather
than reading it, which is the pattern the "one declared panel fact" rule exists to stop. Now
`isTRUE(per_clock[[id]][["normalizes"]])`, which also drops aliases for free (NULL record).

**3. `accessors.R`: four repeated shapes hoisted, and base `stop()` retired.** `required_field()`
/ `optional_field()` / `require_fields()` / `component_tensor()` + `component_named()` absorb the
declared-field, missing-column and single-component patterns that were restated per clock family
(4 sites for the component-tensor pair alone). All ~15 base `stop(..., call. = FALSE)` calls became
`cli_abort(..., call = NULL)` per the CLI-messages invariant; the shared
`CATALOG_BUG` bullet states the "this is a sync bug" suffix once. Every error path was fired by
hand to confirm no cli pluralization or glue fault hides behind the intended diagnostic.

**Also removed as dead:** `clock_cross_sample_at()` (its only caller, `clock_batch_dependent()`,
went in the entry below; `list_clocks()` reads `batch_dependent` off the index, and the
`cross_sample_at` catalog field stays for when row-chunk streaming lands), and
`linear_predictor()`'s returned `used_cols` / `cached` plus two unused `n <- nrow(DNAm)` locals --
no caller read either field, and `observed_panel()$cached` existed only to feed one of them.

---

## 2026-07-24 -- housekeeping: three tracked dev docs again; DESCRIPTION trimmed to what R/ uses

**Decision.** Two reversals from earlier the same day, plus dead-code removal.

- **`dev/sync-boundary-migration.md` is no longer a tracked doc.** The entry below ("plan/code
  reconciliation") added it as a fourth tracked `dev/` doc and said it was allow-listed in
  `.gitignore`. The allow-list line was never written, so the file has been untracked the whole
  time -- and rather than fix the `.gitignore`, we drop the doc from the tracked set. Reason:
  `data-raw/sync.R` is self-documenting, so a second prose restatement of the upstream contract is
  a copy that can only rot against it. Back to three tracked docs (migration-plan, detail-plan,
  DECISIONS); `.gitignore`'s "three canonical planning docs" comment is now accurate as written.
  The file stays on disk, local-only, as history.
- **DESCRIPTION now lists only what ships.** Dropped `digest` and `rlang` from `Imports` (used
  only in `data-raw/sync.R`, which is `.Rbuildignore`d) and moved `withr` to `Suggests` (used only
  by `R/dev-utils.R`, itself `.Rbuildignore`d, and by tests). `devtools` stays undeclared for the
  same reason. `LazyData: true` with no `data/` dir is a known R CMD check warning, left for the
  maintainer.
- **Deleted `clock_batch_dependent()`** from `R/accessors.R` -- a one-line negation of
  `clock_cross_sample_at()` that nothing called; `list_clocks()` reads the index field directly.

**Also fixed in `CLAUDE.md`** (stale, no decision behind them): `curl` listed as a soft dep though
it appears nowhere in `R/` or `DESCRIPTION`, and `res$provenance$batch_set_id` used as the example
of assertable output though that field was removed earlier today.

**Not changed, deliberately.** The four stray `#' @export` tags in `R/` and the hand-edited
`NAMESPACE` are agent leftovers that `document()` will resolve when roxygen is switched on;
`Config/roxygen2/version` stays as the placeholder for that. The no-roxygen invariant holds --
what it does not license is hand-editing `NAMESPACE`.

---

## 2026-07-24 -- sync.R fails explicitly: closed op set, declared path kinds, stated pack order

**Decision.** Follow-up hardening pass over `data-raw/sync.R` after the boundary migration. The
migration removed the places that *re-derived* a declared fact; this pass removes the places that
**absorbed** an undeclared one. Rule applied throughout: upstream is vendored, so anything the
snapshot does not declare must stop the sync naming the clock and the field, never fall back,
default, or skip.

**What changed, and why each was a silent hole.**

- **`declared_path()` / `vendored_path()` replace `startsWith(path, "weights/")` + silent skip.**
  That one filter was doing real policy work invisibly: 9 of the 11 declared `code_ref`s point into
  `papers/`, which R may never read, and every other declared pointer (60 components, 6 shared, 1
  probe_set, 10 imputation refs, 1 code_dep) is already under `weights/`. So the filter's entire job
  was the `papers/` carve-out -- but expressed as an unconditional `return()`, which means a
  `weights/` pointer that is typo'd, moved, or prefixed would vanish with no message. Now `papers/`
  is skipped **by name**, `weights/` is vendored, and a third prefix, an absolute path, an escaping
  path or an empty string is a stop. Also the single place that normalizes separators, replacing
  five scattered `gsub("\\\\", "/")` copies.
- **A declared path carries its KIND from the field that declared it**, ending the
  `grepl("\\.[Rr]$", rel)` sniff in `build_group_bundles()`. The sniff failed in the dangerous
  direction: a component whose file ended `.R` would become an `r_source` blob, yield no row keys,
  and silently reroute the clock to the score-assembled branch.
- **`KNOWN_OPS` closes the recipe-op vocabulary.** Classifying by membership in `CROSS_SAMPLE_OPS`
  only works if the vocabulary is closed -- otherwise a new upstream op defaults to chunk-safe and
  mis-chunks a clock in silence. This is `score_type()`'s "routing is total, a gap is a hard stop"
  applied on the sync side, and it paid for itself immediately: it caught `epitoc2`, a 13th op that
  a hand audit had missed (the enumerating regex did not allow digits). Classified per-sample.
- **`tensor_row_keys()` stops instead of returning `character()`.** `own_scoring_cpgs()` reads empty
  as "this clock declares no cpg-keyed tensor, walk its recipe inputs", so an absent tensor or an
  unrecognised shape was changing a clock's *routing*, not just its data. It also found the key
  column by searching for one named `"cpg"`; `read_tensor_csv()` has already asserted the declared
  `row_key` against the header, so column 1 is the key column by contract.
- **Explicit stops replace fallbacks** on: duplicate `clock_id` / `group_id` across metas (last one
  silently won), a clock or group meta with no id, a manifest row with no `bundle_hash` (the pack
  staleness key -- `lockfile_hit()`'s `anyNA` was absorbing it as "rebuild forever"), a zero-row
  tensor, an unnamed numeric reaching `align_double()` (dead branch that would mis-align a whole
  coefficient column if reached), a missing `release_tag`, `coef_path()` inventing a path from
  `clock_id` when `group_id` was empty, and `resolve_source()` guessing `origin/master` when
  `origin/HEAD` would not resolve -- a guessed branch name for the source of truth.
- **`attach_sex_routed_aliases()` checks the manifest join.** Upstream validates `routing.sex`
  against the *group*; whether both members survived the manifest cut is ours. Neither vendored is
  a staged-not-done family -- skipped with a message. Exactly one is a half-vendored family that
  cannot be scored -- a stop. Previously a NULL donor produced an alias with NULL fields and
  `clock_inputs` pointing at clocks that do not exist.
- **`covariate_names()` rejects a partially named object.** The old `flatten_names()` kept
  `nms[nzchar(nms)]`, so a half-named covariate object silently dropped its unnamed half.

**`pack_canonical_cpgs()`: the order is declared, so do not sort it.** The retired rule took the
canonical probe order from "the longest named numeric tensor", then swapped in any probe list that
`setequal`ed it. Column *values* could not be wrong (`align_double()` `setequal`-checks every
column), but the order feeds `payload_hash`, so the content address was set by a heuristic over
loaded payloads -- and it was fragile in a specific way: a new member with a longer panel would flip
the winner and silently repack the group while no weight had changed.

The first attempt here replaced it with a radix sort, on the assumption that no declared order
exists (SystemsAge ships `_shared/CpGs.csv.gz`, PCClocks does not, PCBrainAge has no group sidecar
at all) and that **which** order is therefore a packaging choice. **That assumption was wrong, and
checking it is what settled the design.** Upstream writes every cpg-keyed tensor in a group in ONE
shared row order -- verified across all 15 PCClocks tensors, both PCBrainAge tensors, and all 24
SystemsAge tensors, with SystemsAge's `_shared/CpGs.csv.gz` in that same order too. Sorting would
have discarded a declared fact and computed a replacement: precisely the class of hack this pass
existed to remove. (PCClocks' declared order happens to already be sorted; the other two are not.)

So the rule is: **the pack order is the order every declared cpg-keyed tensor agrees on, and a
tensor that disagrees -- in set OR in order -- stops.** That is strictly stronger than the
set-equality check the sort would have needed, has no reference tensor to pick and therefore none to
flip, and reads upstream instead of computing. It also reproduces the existing packs bit-for-bit:
`payload_hash` is unchanged for all three (`a70fc3ed` / `5c04ff6f` / `2c180afc`), so
`EXTERNAL_ENCODING_VERSION` stays at **3** and nothing needs re-uploading.

Retained from that first attempt: radix sort wherever an order feeds `payload_hash` but is genuinely
ours -- `member_coef_files()`'s column order and `stable_external_payload()`'s tensor/clock sorts.
Those were pre-existing locale-collation hazards that would make a content address
machine-dependent. They are value-identical on ASCII ids here, which the unchanged hashes confirm.

**Verification.** Dry-run build in memory against the on-disk snapshot, diffed against the committed
`R/sysdata.rda`: 122 clocks, 129 after aliases, **zero** panel diffs, zero catalog field diffs,
`mc_index` identical, zero bundled-tensor diffs, and all three `payload_hash`es identical to the
committed release tags. Every pack column additionally re-verified against the source tensor read
straight off the snapshot and looked up by cpg. Always-on tiers 441 passed, 0 failed, 0 errors, 0
warnings.

**The publish path got the same treatment.** It only runs under `sync(upload = TRUE)`, but it is the
one outward-facing, hard-to-undo action in the package, so a silent fallback there is the worst kind:

- **`package_release_repo()` stops instead of falling through to `origin`.** A set-but-unparseable
  `MC_RELEASE_REPO` used to be discarded silently and the release published to the package's own
  repo. Measured before the fix: `https://github.com/OWNER/other-repo/` (trailing slash) and
  `other-repo` (no owner) both resolved to `hhp94/methylCIPHERv2`, and `a/b/c` produced an
  incoherent `owner=a, repo=c, slug=a/b/c`. `parse_github_owner_repo()` now strips scheme /
  scp-user / `github.com` host / `.git` / trailing slash and then requires exactly two plain
  segments; everything else is `NULL`, and both callers treat `NULL` as a stop. **The owner is
  never inferred** -- a bare repo name is an error, not "same owner, different repo", because the
  release repo is expected to move.
- **`upload_pat()` reads `MC_UPLOAD_PAT` and nothing else.** The `GITHUB_TOKEN` / `GH_TOKEN`
  fallback meant that whatever broad token happened to be in the environment could publish a
  release. `sync.R` still injects the PAT into the child process under those names -- that is
  gh_upload.py's auth contract -- but the *source* is now single.
- **A violated stdout contract stops.** `gh_upload.py` documents stdout as exactly
  `{"results": [...]}`; R parsed it with `tryCatch(..., error = function(e) NULL)` and then looped
  over `res$results %||% list()`, so garbage on stdout printed nothing and reported success. Now an
  unparseable result, or one that does not account for every asset sent, is an error naming what
  could not be confirmed. Re-running is safe -- publishing is idempotent by `release_tag` + name.
  (`$` -> `[[` here too; `res$results` was the documented partial-match hazard.)
- **`gh_upload.py` routes repo resolution through `die_github()`.** `Github()` / `get_repo()` sat
  outside the `try`, so the two most likely failures -- bad token, inaccessible repo -- escaped as
  raw tracebacks while the careful handler (status + response body + hint) went unused. Added a 404
  hint, since GitHub returns 404 rather than 403 for repos a token cannot see. Malformed manifests
  now return the documented exit-2 contract instead of a `KeyError`.

Verified offline with a stubbed `github` module: 401/404/manifest failures all produce clean
diagnostics and no traceback, the happy path and the idempotent skip both emit the documented stdout
object, and an exit-0-with-garbage-stdout run stops on the R side. Left alone deliberately:
gh_upload.py never re-hashes the file it uploads and decides "already present" by asset **name**.
That is the stated split -- R owns content identity, the script only publishes -- and it holds as
long as nothing hand-edits `data-raw/assets/`.

**Not done here.** `R/sysdata.rda` is not regenerated -- nothing it holds would change (the dry run
proves it), and regenerating is a maintainer step that also stamps `source_git_sha`. Also left
alone: the vendored `r_source` blobs (`weights/MiAge/*.R`,
`weights/prcPhenoAge/score_prcphenoage.R`) ship in `mc_bundles` but nothing in `R/` reads
`r_source` -- MiAge is scored by our own `R/score_MiAge.R`. Dead payload, not a correctness bug.

---

## 2026-07-24 -- plan/code reconciliation after the boundary migration

**Decision.** Audited `CLAUDE.md`, `dev/migration-plan.md` and `dev/detail-plan.md` against the code
after the boundary migration and fixed every disagreement (code is truth). Also added a fourth
tracked `dev/` doc, `dev/sync-boundary-migration.md`, and allow-listed it in `.gitignore`: it is the
only copy of the upstream contract a collaborator can read, since `data-raw/methylCIPHER-meta/` is
gitignored.

**What was stale.** `CLAUDE.md`: the `CUSTOM_GROUPS` registry (deleted) and the `source_git_sha`
lockfile key (now `bundle_hash`). `migration-plan.md`: `depends_on_clocks`, `external_package` as a
live one-off, `CUSTOM_GROUPS`, and the single-cohort parity row. `detail-plan.md`: the same
`depends_on_clocks` / `CUSTOM_GROUPS` / `source_git_sha` items, plus the three-tier panel-resolution
description, the alias's invented `(weights_format = "routed", computation_type = "sex_routed")`
pair, the `sample_scale` batch-op marker (doubly wrong -- `batch_ops` is gone *and* `sample_scale`
was never cross-sample), the `cohort_EPIC` fixture paths, and a follow-up describing a scale-aware
tolerance whose motivating skip class no longer exists.

**Why record it.** Most of these were not drift from neglect -- they were correct until this
session's migration retired the thing they described, and one (`CLAUDE.md`'s sync workflow) was
made stale by the migration itself while its Testing section was updated in the same sitting. The
pattern worth remembering: a `sync.R` change can invalidate prose in three files that never mention
`sync.R`, because the catalog fields it emits are the plans' vocabulary. Grep the retired
identifiers across all tracked docs, not just the one being edited.

**Deliberately not changed.** `CLAUDE.md` still lists `augment` among the mc_result verbs. Only
`rbind.mc_result` exists today; the rest is stated intent for in-flight work
(`tests/testthat/test-mc-result.R`), not a claim about the current build. `dev/id-streaming-plan.md`
(local-only) now notes the gap explicitly.

## 2026-07-24 -- sync.R migrated to the SoT boundary contract (upstream declares, we classify)

**Decision.** Reworked `data-raw/sync.R` to read declared upstream facts instead of re-deriving
them, following `control/DOWNSTREAM_SYNC_MIGRATION.md` (the upstream hand-off, tasks D1-D12).
Retired downstream: the whole-JSON regex path discovery (`WEIGHTS_REF_RE`,
`collect_weights_refs`), column-count shape sniffing in `read_tensor_csv()`, the covariate
archaeology (`linear_sex` / `sex_params` / `sex_stratified` / sex-keyed `imputation.ref`), the
three-tier `covers`-fallback fixpoint in `resolve_group_scoring_probe_sets()`, the MiAge special
case (`CUSTOM_GROUPS`, `miage_site_parameters`, `attach_custom_components`, positional
`(b, c, d)` alignment), `EXTERNAL_FIELDS` / the `external` pin object, and the hardcoded external
encoders (`SYSTEMSAGE_ORGANS`, the PCClocks name-grep, PCBrainAge's two literal paths).

Replacing them: `declared_tensors()` reads the four named pointer fields plus `code_ref` /
`code_deps` and applies the derived `coef_path()` rule; tensor shape is asserted against the
declaring `row_key` / `col_key` (`col_key` is optional -- a rotation matrix's PC columns are
generated, so only `row_key` is checkable there); the scoring panel is a DAG walk over recipe
`inputs`; and external packs are encoded from declared components plus the group's member
`coef_path`s, sorted by `clock_id`.

**Why this is safe rather than a rewrite-and-hope.** Two independent gates, both green:
`assert_declared_n_cpgs()` proves every one of the 122 clocks' derived panel equals its declared
`n_cpgs` with no exemption list (the retired fixpoint gave `DNAmFitAge_Female` the family-wide
627-probe prep panel against a declared 172), and the rebuilt external packs are **bit-identical**
to the old hardcoded encoding -- same matrices, same column order, same `payload_hash` -- so
`EXTERNAL_ENCODING_VERSION` did not need to move.

**Two corrections the D-list did not cover.** Upstream also split `stack` operands into three
disjoint namespaces (`inputs` = other clocks, `internal` = a prior step's `out`, `covariates`).
`physage_surrogates()` and `systemsage_stack_order()` both read `[["inputs"]]`, which is now
`NULL` for PhysAge and SystemsAge -- silently breaking both composites. They now share
`stack_operands()`, which concatenates the three declared lists in the documented key order.

**Enum hygiene.** Sex-routed aliases no longer claim invented `weights_format = "routed"` /
`computation_type = "sex_routed"`; they carry `kind = "sex_routed_alias"` and leave both upstream
enums `NA`, so `mc_index` only ever holds values from `meta_schema.py`. `score_type()` and the
coverage split check `kind` before reading either enum. The dependency edge moved from the retired
`depends_on_clocks` meta field to a derived `clock_inputs` (recipe `inputs`, or an alias's two
members); `clock_depends_on()` reads that.

**Sample-axis.** `cohort_zscore` is the only genuinely cross-sample op. `sample_scale` is a
*within-sample* z-score (Zhang2019) and was wrongly flagged batch-dependent; fixed-target
normalizations (DunedinPACE -> `gold_standard_means`, Horvath-pipeline BMIQ, `noob`) are
per_sample too. The stored fact is now `cross_sample_at` (the first cross-sample step's position,
so row-chunk streaming knows where it may stop), and `batch_dependent` is derived from it in
`build_index()` rather than stored beside it.

## 2026-07-24 -- parity runs on every cohort, and its tolerance is ours, keyed on a declared field

**Decision.** `test-fixtures-parity.R` now runs each (clock, cohort) pair upstream declares a
`fixtures[]` block for -- `cohort_EPICv1` and `cohort_450K`, one duckdb connection each. Grading:
a fixture whose `server_normalization` is non-empty is held to **correlation** (> 0.99); every
other fixture is held to **exact** (max_abs_diff < 1e-6).

**Why keyed on `server_normalization` rather than a clock list.** Upstream retired its `parity`
policy enum -- tolerance is downstream's (sec 12) -- so we needed our own rule, and a hand-kept
list of "clocks that cannot be exact" would rot the moment a clock's oracle changed. The declared
field already draws exactly the right line: it is non-empty iff the golden was produced by the
Horvath calculator on betas the **server** normalized (BMIQ / noob), and we score raw betas and
deliberately vendor no BMIQ gold standard (there is none upstream -- its provenance is unknowable,
so the SoT pins the fixture instead). It partitions the 244 fixture blocks cleanly: 186
`author_code` / `frozen_reference` blocks saw the same raw betas we do and are the exact gate; 58
`horvath_online` blocks did not.

**Why the failures this exposed were not regressions.** Running the tier for the first time against
the new upstream produced 39 failures, and every single one had `oracle: horvath_online` -- zero
from the other two oracles. Cross-checking the previous committed `sysdata.rda`, 18 of the 22
distinct clocks had been graded `correlation` or `skipped` by upstream's own retired policy, i.e.
never held to exact tolerance; the remaining 4 failed only on `cohort_450K`, a cohort that did not
exist before. So the migration changed no number that was previously being checked.

**Skip list shrank 8 -> 1.** Seven inherited `KNOWN_PARITY_GAPS` entries were GrimAge/FitAge
surrogates skipped for sub-1e-5 exact-tolerance drift; under the correct grading they pass
outright and are now real tests. Only `Zhang2019` remains (its `sample_scale` moments are taken
over the needed-CpG subset rather than the full panel, so exact parity is unreachable by
construction).

---

## 2026-07-24 -- `batch_set_id` provenance field removed

**Decision.** Dropped `$provenance$batch_set_id` and its computation in `calc_clocks()`
(`digest(sort(sample_id))` when any clock was `batch_dependent`). `digest` is no longer used
anywhere in `R/`. The `batch_dependent` catalog flag and `clock_batch_dependent()` accessor stay --
they still drive `list_clocks()` discovery and the future kind-1/kind-2 split
(`dev/id-streaming-plan.md` Phase 2).

**Why.** `batch_set_id` was a fingerprint of the record's own `sample_id`, and the full id vector
already lives at `$provenance$sample_id`. Comparing two `batch_set_id`s is identical to
`setequal()` on the ids, so it added nothing: `cbind` compatibility reads the id set directly
(equal sets), and a future streaming `rbind` wants disjoint sets -- both computable from the ids.
Worse, it hashed the *current* ids, so after a `[` row-subset it could not detect the one thing it
existed for ("same current ids, but the batch-dependent column was computed over a different
cohort"): it would silently agree. That cohort-mismatch concern is dissolved anyway by Phase 3
(move every cross-sample op to `augment`, leaving the scoring loop batch-invariant, so no cohort is
baked into a score to fingerprint). Surfaced while reasoning through the streaming plan's rbind
gates; see the 2026-07-24 mandatory-rownames entry for the sibling cleanup.

**Blast radius.** `calc_clocks.R` (computation + `construct_mc_result` param + provenance field);
`test-score-physage.R` and `test-score-fitage.R` (dropped the `batch_set_id` assertions);
`detail-plan.md` sec 1.3/6/7/7.1 (`cbind` now two gates, no batch-cohort gate). This reverses the
sec 6 / sec 7.1-gate-3 "store sample-set id, cbind rejects incompatible sets" design.

## 2026-07-24 -- DNAm rownames are mandatory; positional/synthetic ids deleted

**Decision.** `calc_clocks()` now requires `rownames(DNAm)` and never invents sample ids. This
reverses the 2026-07-17 design (rowname-less DNAm -> stamp `sample1..N`, mark non-bindable, refuse
at cbind). Gone: the `allow_positional_ids` param, the positional stamp block, the `positional`
plumbing threaded through `check_pheno()` / `warn_missing_covariates()` / `resolve_pheno()` (the
row-order-join branch is deleted), and the write-only `$provenance$positional_ids` field (nothing
read it). `check_DNAm()` already asserted non-NULL unique rownames -- the stamp existed only to
satisfy that assertion before it fired, so this is a net deletion; `check_DNAm()` now leads with a
cli message handing the user the one-liner to name anonymous rows themselves.

**Why.** The positional apparatus existed to guess ids and then carry a `synthetic`/`bindable` flag
everywhere to remember the guess was meaningless. A real methylation matrix essentially always has
sample identifiers; a nameless one is a degenerate input, and the honest response is an error, not a
silent guess. Mandating rownames does not remove the `sample1..N` capability -- it relocates it to a
caller-written one-liner (`rownames(DNAm) <- paste0("sample", seq_len(nrow(DNAm)))`), moving the
"these rows are positional" admission to where the knowledge is. It also kills a genuine footgun in
the planned id-resolution contract (adopting `pheno[[id]]` onto rowname-less DNAm by row order,
which nrow-checks length but not order). Pre-alpha, so the public API break is acceptable.

**Streaming (future, `dev/id-streaming-plan.md`).** This subsumes that plan's Phase 1. The Phase 4
rbind refusal loses its `synthetic_ids` gate, but the id-collision gate it already lists catches the
same failure mode (lazy per-chunk `sample1..N` collides on `sample1`). Phases 2-4 stand unchanged.

## 2026-07-24 -- housekeeping pass: rbind refuses, NAMESPACE trimmed, pick_one re-hoists "exactly one"

**Decision.** A trim/hoist sweep, scores untouched (442 always-on tests green before and after):

- **`rbind.mc_result` now refuses** via `cli_abort` instead of an empty `TODO` stub that silently
  returned `NULL`. This is the documented invariant ("rbind refuses"); it was registered but unimplemented.
- **NAMESPACE:** dropped `S3method(as.matrix,mc_result)` and `S3method(print,mc_result)` -- neither
  method is defined yet, and R warns "declared in NAMESPACE but not found" at load. Re-add each line
  when its method lands. Also corrected the false "Generated by roxygen2" header (NAMESPACE is
  hand-managed, see CLAUDE.md).
- **Dead params trimmed:** `random_betas(seed=)` (never passed; contradicts the unseeded-inputs test
  policy) and `build_partial_cache(cores=)` (never passed; always 1L, now inlined). `%||% NULL`
  no-op in `accessors.R` removed.
- **`pick_one()` re-hoisted.** The ~11 "`Filter(pred, ...)` -> `length != 1L` -> `stop()`" blocks in
  `accessors.R` collapse into one `pick_one(items, pred, what, id)` that cli-stops when a declared
  pointer is absent or ambiguous. This looks like the `only_one` helper removed from `score_GrimAge`
  earlier today, but that removal was a base-`stop` -> cli migration, not a rejection of consolidation:
  `pick_one` keeps the consolidation *and* is cli-native. `miage_params` gains the exactly-one check it
  lacked. Sites with non-strict counts (`recipe_step_out`, `physage_poly_coef`, allow 0) stay as-is.
  `score_GrimAge`'s inline site was left untouched (active uncommitted work).

**Why now.** User-requested housekeeping ("scope creeped a bit"). Kept the fuzzy `did_you_mean`
suggestion machinery in `resolve_inputs.R` -- helpful, tested, low-maintenance.

---

## 2026-07-24 -- coverage_table(x, by) splits into clocks_coverage() + samples_coverage()

**Decision.** The single `coverage_table(x, by = c("clock", "sample"))` planned in `detail-plan.md`
sec 4 and named in CLAUDE.md is built instead as **two** plain functions: `clocks_coverage(x)` (the
per-clock wide table) and `samples_coverage(x)` (a per-(sample, clock, panel) long table). Both are
pure formatters over `$coverage` + `$provenance` in `R/coverage_report.R`; neither re-touches beta.
This is the last piece of the coverage hoist (the three data-model commits landed earlier today).

**Why two, not one `by`.** The two granularities do not share a shape. The clock table is wide and
sample-invariant (one row per clock, panel set sizes). The sample table must be **long**, keyed
`(id, clock_id, panel)`, because only long carries a per-row denominator -- a sex-routed alias needs
`n_needed` = 172 for a woman and 190 for a man on the same column, taken from the member that scored
that row. A `by` switch returning two unrelated frame shapes from one name buys nothing over two
named functions, and two names document themselves.

**Coverage is literally `row_coverage()`.** `samples_coverage()$coverage` on a clock's row-gate panel
is the value `row_coverage()` returns (`R/resolve_inputs.R`), the same source `check_row_coverage()`
warns from, so the table and the `min_samples_coverage` warning cannot disagree. For an alias (NULL
record) the ratio is stitched from each member's `row_coverage()` over its own masked rows; for the
non-gate score panel of a normalizing clock (DunedinPACE) it is the same `(present - miss) / needed`
formula on the score panel, which has no warning to disagree with.

**`role`, two values.** `clocks_coverage()` tags each row `returned` (any score column, including a
sex-routed alias whose panel fields are NA and dependency columns like `GrimAgeV1`) or
`routing_target` (a member kept for coverage, never a column). A third `dependency` value was
considered and dropped -- `provenance$dependencies` already records that fact, so duplicating it in
the table adds no information. Per-clock imputation counts live here (per-clock constants, next to
`policy`), not in `samples_coverage`.

**Scores untouched.** Formatters read a finished record only; parity is unaffected by construction.

---

## 2026-07-24 -- per-panel per-sample counts; the norm/score panel is one declared fact

**Decision.** Commit 3 of the coverage hoist. Per-sample miss is now counted **per panel role**, not
over one implicitly-chosen panel: every clock has a score panel; a normalizing clock (only
DunedinPACE today) also has a norm panel. `$coverage$sample_miss` changes from an `n x k` matrix to
`list(score = <n x k>, norm = <n x k' over just the normalizing columns>)`, and the per-clock record
splits `score_imputed_partial` (score panel) from a new `norm_imputed_partial` (norm panel). The
score-panel per-sample count for DunedinPACE, which no code computed before, now exists. This is a
**public output change** (chosen deliberately over the output-stable option), so the two
`$coverage$sample_miss[, clock]` test references moved to `$score`.

**Why now / what it enables.** A future normalization group (BMIQ) must report norm-panel
missingness separately from score-panel missingness, and the old single implicit-panel vector could
not express it. The data model now can. It also removes real duplicated work: counts ride
`resolve_cpgs()`'s existing distinct-panel index (`panel_index`) and are computed once per distinct
panel, then fanned out -- FitAge's 7 aliases + 14 members and GrimAge's shared component panels no
longer recount identical panels.

**One declared panel fact (sec 3.3 drift closed).** "Does this clock count over a norm panel" was
derived three independent ways -- `score_Dunedin` (`clock_norm_scheme == "quantile"`),
`clock_coverage` (`length(norm_needed) > 0`), `row_coverage` (`norm_needed > 0`). `resolve_cpgs()`
now computes `normalizes` once into each per_clock entry; `coverage_record()` copies it onto the
record; `score_Dunedin` (panel + fill ref) and `row_coverage` (row-gate denominator) both read it.
No reader re-derives it. `count_sample_miss()` was also normalized to always return integer, so every
miss vector and partial-fill sum is integer end to end.

**Scores untouched.** Parity stays green (0 fail); the widening moves only which panel a count is
attributed to, never a coefficient or a score.

---

## 2026-07-24 -- coverage is hoisted upstream; a branch returns only its score

**Decision.** Reverses the "every branch returns the same record `list(score, coverage,
sample_miss)`" invariant. Coverage/QC depends on no score, so it is now computed once, upstream of
the scoring loop, by `compute_coverage()` (new `R/coverage.R`), keyed by clock id, and merged in
`construct_mc_result()`. Every scoring branch returns just its `<n x 1>` score matrix; `results` is
a list of score matrices, and composites read a dependency's score as `as.numeric(results[[dep]])`.
This is commit 2 of the coverage hoist (gitignored `dev/coverage_hoist_plan.md`), and it depends on
commit 1 having made the counts a pure function of policy.

**Why it is now possible / what it buys.** After the 2026-07-24 `used` unification, every branch's
coverage counts were the same function of `(policy, present, absent)`, computed inline in ten
places. Hoisting collapses those ten copies to one: `clock_coverage()` counts `sample_miss` over the
norm panel when the clock normalizes (only DunedinPACE) else the score panel, and derives
used/imputed_full/dropped from `clock_impute(id)[["policy"]]` (`vendor_mean` fills, else drops). It
does **not** move the throwing gate earlier -- `check_coverage()` already ran before scoring; the
win is one source of truth and a data model that can later express per-panel counts (commit 3).

**What stayed where.** Scores are untouched, so parity is unaffected by construction, and the 33
output-level coverage/`sample_miss` assertions pass byte-identical. The sex-routed split is
preserved exactly: an alias keeps a `NULL` per_clock record and a per-sample miss **stitched** from
its members (each row's count is the count of the one member that scored it), and routed members
are **masked** to their own samples (wrong-sex rows blanked, `score_imputed_partial` recomputed) --
this logic moved from `mask_routed_members()`/`score_sex_routed()` into `compute_coverage()`; the
now-inert score-masking of members (never a column, never read post-loop) was dropped.

**sec 3.2 resolved package-side.** `pack_linear_coverage()`/`systemsage_composite_coverage()` (which
hardcoded `policy = "vendor_mean"`) are gone; the one upstream computation reads the catalog policy,
confirmed non-NULL for every record-building clock (see 2026-07-24 `used` entry). Dunedin's fill-ref
selection stays in its scorer untouched -- its `vendor_mean` label names the QN-gold mechanism, not
`clock_impute_ref()`, and that fold-in waits for commit 3.

---

## 2026-07-24 -- score_used counts vendor-filled CpGs everywhere (used = present + imputed_full)

**Decision.** `score_used` in every coverage record now means "terms that entered the linear
predictor" = `score_present + score_imputed_full`. A vendor-mean-filled absent CpG is counted as
used, because it does enter the sum; an omitted/dropped CpG is not. This is commit 1 of the coverage
hoist (see gitignored `dev/coverage_hoist_plan.md` sec 3.1 / 8); centralizing coverage in commit 2
needs one definition, so the divergent ones are unified first, alone, as a visible behavior change.

**Why / what actually changed.** The engine (`linear_score`) already defined `used = present +
vendor_filled`. Auditing all ten branches against that, only **two** disagreed -- the plan's sec 3.1
overstated it as "the composites": `score_DNAmFitAge` (the KDM composite,
`DNAmFitAge_{Female,Male}`) and `systemsage_composite_coverage` (`Age_prediction`, `SystemsAge`)
both reported `used = present` while also reporting `imputed_full = absent`, so `used != present +
imputed_full` exactly when a fill happened. `score_PhysAge` already reported `used = score_needed`
(correct); `score_GrimAge`/`EpiTOC2`/`MiAge`/`Zhang2019` omit (imputed_full = 0, so `used = present`
is already `present + 0`); `pack_linear_coverage` and `Dunedin` were already `present + absent`.
Fixed the two offenders; `score_used` is now asserted (it was asserted nowhere before) via an
always-on FitAge composite golden and the SystemsAge external test. Scores are untouched, so the
parity tier is unaffected by construction.

---

## 2026-07-23 -- routed members leave the output too; the alias is the family's only column

**Decision.** Reverses the "still columns in the result" half of 2026-07-22. The 14 sex-resolved
DNAmFitAge members are now internal end to end: scored, consumed by their alias, kept as coverage
rows, and **dropped from `output_ids`**. A sex-routed family comes back as one column per alias,
finite for every sample in `pheno`. `drop_routed_members()` sits next to `sex_routed_members()` so
the callable pool and the output filter read the same table.

**Why.** Asking for `DNAmFitAge` on a mixed cohort returned 30 columns, 14 of them half-`NA`:
`DNAmFitAge_Female` populated for the women, `DNAmFitAge_Male` for the men, and the alias --
the actual answer -- buried among them. Every one of those NAs is structural, not a data problem,
but nothing in the object says so, so the shape reads as a failed run. It also forces every
downstream consumer (`as.data.frame()`, a `cbind` into a phenotype table, a regression) to know
which columns to ignore. The masking from 2026-07-22 fixed the *values* footgun -- no plausible
wrong-sex numbers -- without fixing the *shape* one. Members were only ever columns because they
were the place coverage lived, which is a reason for them to be coverage rows, not columns.

**What crosses the routing split is per-sample, not per-panel.** This is the line, and it is not
"scores are stitched, coverage is not" -- we tried that phrasing and it cuts in the wrong place.

A routed sample was scored by **exactly one** member. So anything indexed *by sample* routes with
the score and remains true of its row: `$scores`, and the tier-2 `sample_miss`, both stitched in
`score_sex_routed()`. Anything indexed *by panel* does not, because the two panels differ (172 F /
190 M) and such a count is only meaningful against the panel it was counted over -- "4 imputed" is
4 of 172 for a woman and 4 of 190 for a man. Merge those numerators into one figure while the
denominators stay on two separate member rows and you get a number nobody can divide. So
`per_clock[[alias]]` stays `NULL`, the panel counts stay on the member rows -- which now outlive
their columns, `construct_mc_result()` taking the full `results` for `per_clock` and only
`results[output_ids]` for `$scores` / `$sample_miss` -- and the `NULL` is not a gap to fill in. It
is the correct answer to "what is this family's panel coverage" when the honest reply is "it
depends which model scored you, and here are both rows."

The residual sharpness, worth stating so nobody trips on it: a stitched `sample_miss` column has a
**per-row denominator**. Comparing 4 against 4 down that column is not comparing like with like if
the two rows routed to different members. The column is a per-sample QC flag ("this score leaned
on N imputed CpGs"), not a coverage rate; the rate needs the member's `per_clock` row.

**`check_row_coverage()` therefore runs over every clock computed, not `results[output_ids]`.**
The row floor *is* a per-panel division, so it can only be evaluated on the member records -- and
those have no columns. Gating the returned columns alone would have silently switched the floor
off for the whole family. A warning can name `DNAmGrip_wAge_Female`, which is not a column --
deliberate: it names the model that was thin and the `per_clock` row holding its numbers, and
`list_clocks()` still shows it with `callable = FALSE`. Masking (`mask_routed_members()`) is kept
for the same reason -- it is no longer about output values, it is what keeps a member's coverage
describing only its own sex.

**`batch_set_id` moved to `clock_sequence`** for the same class of reason: batch dependence is a
property of what was computed, not of what is returned. No behavior change today (the two sets
coincided until now, and no FitAge clock is batch-dependent), but a batch-dependent routed member
would otherwise have lost the stamp on the way out.

**Not extended to ordinary dependencies.** `GrimAgeV1` and its 8 surrogates still come back as
columns from a `DNAmFitAge` request. They are callable clocks with their own meaning and their own
coverage, and a caller who wanted only the alias can subset. The rule is narrow: *routing targets*
are not results.

---

## 2026-07-23 -- the record carries the aligned pheno; rbind is refused, not built

**Decision.** Three parts of the `mc_result` S3 surface, settled together as Phase 1 of the
verb build-out.

**1. The record carries the aligned pheno at top-level `$pheno`**, narrowed by
`resolve_pheno()` to `c(pheno_id, <covariate union>)` -- for the current catalog that is at most
`ID`, `Age`, `Female`, since those are the only two covariates any clock requires. Columns the
user supplied but no clock needed are dropped, so an arbitrary clinical table cannot ride along.

We first designed this as a per-covariate `digest()` keyed by `sample_id`, to keep subject data
out of the object. Values win on three counts. A digest **cannot be recomputed after a row
subset**, so any future `[` would silently void the `cbind` shared-covariate gate (sec 7.1 gate 2)
for subset records, which is exactly where a mistake is most likely. A digest can only report
*that* `Age` disagrees, never which samples. And a saved record could not answer "what Age
produced this score", which matters when a clock fed age in months looks merely implausible
rather than wrong. Values are the strictly stronger representation -- a digest is derivable from
them, never the reverse -- so a comparison can always reduce both sides to the weaker form, and a
future scrubber that swaps values for digests stays compatible with the gates.

**The privacy objection is real but weaker than it looks here.** `$provenance$sample_id` already
keys the record to individuals, and the scores are themselves health data, so `Age` / `Female` is
a marginal rather than categorical increase. More decisively, the export path out of this object
(`as.data.frame()`) is **scores plus id only** -- pheno never rides it. Leaking pheno requires
deliberately `saveRDS()`-ing the whole record, and `print.mc_result` shows a `<$pheno>` block so
nobody can be unaware it is in there. Keep `as.data.frame()` scores-only; that is what makes the
rest of this sound.

**Scores-only is unchanged.** The invariant governs `$scores` and the `as.data.frame()` shape, not
whether the record remembers its own inputs.

**2. `rbind.mc_result` refuses.** It was a registered stub with an empty body, so `rbind(a, b)`
returned `NULL` silently. Deleting the method is not the fix: `rbind` is an internal generic, so
with no method the default fires and yields a 2 x N list-matrix -- still garbage, merely
differently shaped. It is now a `cli_abort` naming the reason: binding samples orphans
`batch_set_id` and leaves any cohort-dependent score (`cohort_zscore`, and residual acceleration
later) computed against the wrong cohort. The supported path is one `calc_clocks()` per cohort,
then bind the `as.data.frame()` outputs -- the same division of labour broom draws, where three
`lm`s mean three `augment`s and the caller owns the bind.

**3. `summary()` is dropped in favour of `coverage_table()`** (not yet built). `summary()` on an
object wrapping an n x k numeric matrix sets up the expectation of score distributions, and
handing back a CpG-coverage QC table instead spends the conventional name on the unexpected
thing, leaving nothing for the expected one. `coverage_table(x, by = c("clock", "sample"))` also
fits the record better: `$coverage` genuinely holds **two** shapes -- `per_clock` and the
samples-x-clocks `sample_miss` -- and one function with an orientation argument covers both.
This reverses the sec 1.3 method table, which listed `summary`.

**Also corrected while here.** `check_DNAm()` ran *before* the positional-id stamp, and it asserts
non-NULL rownames -- so rowname-less DNAm always errored there and `allow_positional_ids = TRUE`,
plus everything downstream of `$provenance$positional_ids`, was unreachable. The stamp now runs
first, guarded by `is.matrix()` so a non-matrix still gets `check_DNAm()`'s error rather than a
confusing failure inside the rownames assignment. Detail-plan sec 5.1 already specified this
order; the code had it backwards.

`citation` / `codebook` are **not** methods: `utils::citation` is a plain function, not an S3
generic, so `citation.mc_result` would never dispatch and defining our own generic would mask it.
They become plain functions if built. Deferred, along with `augment` and `[`.

---

## 2026-07-23 -- cli pluralization is a hard rule, not a style preference

**Decision.** CLAUDE.md now carries a "CLI messages" section requiring an explicit `cli::qty()`
on any `{?}` marker whose quantity is not the interpolation immediately before it.

**Why this is an entry.** Four shipped sites threw from inside cli instead of reporting the
condition they were written to report -- `score_default.R` and `score_pack.R` ("Add {?it/them} to
`pheno`"), `missingness.R` ("Remove or fix {?it/them}"), and `mc_data.R` ("Pack{?s} ... not found
in {.path {dir}}", which had two forward candidates). A user hitting a missing-pheno-column or a
missing-pack error saw `Cannot pluralize without a quantity` and no diagnostic at all. None of
this was caught because the tier that covers these paths uses bare `expect_error(expr)` per the
test-altitude rule, and a crash inside the error handler is still an error.

**The silent variant is the one to watch.** `"{.val {id}} needs pheno column{?s} {.field {need}}"`
binds the marker to `id`, not `need`, so it read "needs pheno column Age and Female" forever.
Same shape in `resolve_inputs.R`'s "points at missing token{?s}". These do not fail any test and
never will -- only reading the rendered string catches them.

**Not fixed by a test.** Asserting message wording would violate the test-altitude rule and would
not generalize to the next message anyway. The rule lives in CLAUDE.md instead, where it is read
before the message is written.

**Also recorded here: cli reflows whitespace.** `mc_manifest()` builds a column-aligned size table
for the interactive `askYesNo()` prompt; interpolating it into a `cli_abort()` bullet collapsed it
to a single run-on line. The non-interactive consent aborts now use `mc_manifest_bullets()` (one
bullet per pack plus a total) and `mc_manifest()` stays for the plain-text prompt where the
alignment survives. Two renderings of one manifest is deliberate, not duplication.

---

## 2026-07-23 -- the norm-panel-warns-only decision is finally in the code (it never was)

**Decision.** `check_coverage()` now does what the "scoring panel stops, normalization panel warns
only" entry below already decided: `classify()` grades the **scoring** panel alone (stop under
`min_clocks_coverage` or at 0 observed, warn within 10% of the floor), and a thin **normalization**
background is a second, warn-only pass that cannot stop the call. `detail-plan.md`'s column-floor
section is corrected with it.

**Why this is an entry and not a silent fix.** The decision was written and the code was not
touched -- `classify()` kept taking the worst of `{scoring, normalization}`, so a clock whose norm
panel fell below the floor stopped exactly as before. It survived because the only clock with a
norm panel is DunedinPACE, and the one test that covers the case
(`test-coverage-gate.R`, "a sparse normalization panel warns but still scores") is behind
`skip_if_not_installed("betanorm")`, so most runs never executed it. Recording the reconciliation
rather than editing the old entry keeps the gap visible: a decision entry is not evidence the code
follows it, and a skip-gated test is not evidence either.

**Also fixed here.** `check_row_coverage()`'s message printed `score_needed` as its denominator
while `row_coverage()` had already computed the fraction over the `norm_*` pair for a normalizing
clock -- "worst 50.0% of 173 CpGs" when the 50% was 10000/20000 gold CpGs. `row_coverage()` now
returns the denominator it used alongside the fraction, so the number and the panel it names
always agree.

**No behavior change for any other clock.** Every non-Dunedin clock has an empty
`clock_norm_cpgs()`, so it was already graded on its scoring panel alone.

---

## 2026-07-23 -- `clocks` discovery: cli, nearest-match suggestions, `list_clocks()`

`clocks` is one argument over 129 catalog entries in two namespaces, which made a typo a dead end.
Three decisions, in the order they were argued.

- **No token normalization -- resolution stays exact.** Case folding was the tempting fix
  (`"systemsage"` -> `SystemsAge`) and is unsafe here: `prcPhenoAge` (a group of 2) and
  `PRCPhenoAge` (a clock) differ **only** in case, so a case-folded resolver has to silently pick
  between scoring 1 clock and scoring 2. That is live in today's catalog, not hypothetical. It
  would also open a second namespace to keep collision-free on every sync, forever, and buys
  nothing in output (column names come from the catalog either way). Case folding therefore lives
  in `did_you_mean()` and nowhere else: `"SYSTEMSAGE"` earns a suggestion, never a resolution.
- **A near miss is an error, never an auto-correct.** No "assuming you meant X" + warning: warnings
  get swallowed in scripts and knitted output, and the failure mode is scoring the wrong clock on
  someone's cohort. `resolve_clocks()` still batch-collects every bad token and errors once, so
  `c(good, bad1, bad2)` reports both.
- **Nearest-k with no distance cutoff, and the two namespaces kept apart.** A threshold was tried
  first and needed tuning per token length (`"Horvat"` found nothing at `floor(nchar/4)`); ranking
  and taking the top 5 drops the knob entirely. `utils::adist(partial = TRUE)` -- the token scored
  against the closest *substring* of each candidate -- is what makes a half-remembered fragment
  work: under full distance `"Zhang"`, `"PACE"` and `"cortical"` all returned pure noise, because
  plain Levenshtein charges for the length gap. Ties break toward the shorter id so a plain name
  outranks its decorated relatives (`phenoage` -> `PhenoAge`, not `HRSInChPhenoAge`). Clocks and
  groups get their own pool and their own line, because they resolve differently -- and the split
  is what makes `"grimage"` legible: `GrimAge` is *only* a group id, `GrimAgeV1` is the clock.
  `"all"` is deliberately not in either pool; nobody mistypes it into a clock name and it only
  crowds out real candidates.

**Group-vs-clock ambiguity is documented, not fixed.** 27 of 40 group ids are also clock ids; for
23 the group is the singleton `{itself}`. Four are not -- `SystemsAge` (13), `DNAmFitAge` (7),
`EpiTOC2` (2), `RepliTali` (2) -- so those tokens score the group, not the same-named clock. No
`clock:` / `group:` prefix grammar was added: the group is always a superset containing the
same-named clock, its members share a panel/pack so the extra columns are nearly free, and
`mc_result` is subsettable. `list_clocks()` surfaces it as a `group_size` column instead.

`list_clocks()` derives from the same pool `resolve_clocks()` accepts, so what is listed is what is
callable; the 14 sex-routed members appear with `callable = FALSE` and `request_as` naming their
alias.

**Keywords (`MC_TAGS`) are a runtime registry, not a sync artifact.** `"all"` was already a keyword
that expands to a set, so `gestational` / `mitotic` / `mortality` generalize an existing special
case rather than adding a mechanism. They live in `R/tags.R` and **not** in `sysdata.rda` because a
tag is not part of the scoring contract -- no coefficient, panel, or routing decision reads one, so
routing it through `sync()` would make a keyword typo a maintainer-only round-trip. The hard-stop
property that a build step would have bought is bought instead by `resolve_clocks()` erroring on a
tag member the catalog no longer has, plus a test that every tag is non-empty and every member
resolves. Membership is **hardcoded package-side for now and marked `TODO`**: which clocks are
"mitotic" is a scientific claim and belongs upstream in `methylCIPHER-meta` beside `papers.csv` /
`clocks.bib`. Starting package-side is reversible -- the resolver reads a declaration either way,
so moving the source is a `sync()` change, not a resolver change. Because the tag set is small and
closed, the unknown-clock error **enumerates** the keywords instead of ranking them: one footer
line for the whole set, rather than a third suggestion line per bad token.

**cli + rlang added to Imports.** `cli::cli_abort()` / `cli_warn()` / `cli_inform()` replace
`stop()` / `warning()` / `message()`, with `call = NULL` throughout to reproduce the old
`call. = FALSE` (the alternative leaks internal helper names like `check_coverage()` into
user-facing errors). `rlang` is required because `cli` only **Suggests** it and `cli_abort()`
dispatches to `rlang::abort()`. Converted so far: `resolve_inputs.R`, `calc_clocks.R`, `utils.R`.
`accessors.R`, `mc_data.R`, `missingness.R`, `score_*.R` and `sim_DNAm.R` still use base
conditions -- a mechanical sweep, deliberately left to its own commit.

---

## 2026-07-23 -- plan/code reconciliation after the coverage-floor work

Swept `dev/detail-plan.md` against the code (every backticked `fn()` in the plans checked for a
definition in `R/`). Code is truth; the plan was corrected. Recorded here because two of these are
easy to re-introduce.

- **A composite's panel is the declared one, not the closure over its dependencies.** The plan said
  `DNAmFitAge_{Sex}` qualifies to record coverage because "the union panel describes each sample
  truthfully", which reads as *union including `GrimAgeV1`*. It is not: the declared panel is 172
  (F) / 190 (M) CpGs, exactly the union of the three fitness members, sharing **one** CpG with
  `GrimAgeV1`'s 1030. The composite still qualifies (all components hit all samples), but its
  coverage numbers describe the fitness panel only, and the `GrimAgeV1` share is accounted on
  `GrimAgeV1`'s own dependency column. Both floors read the same declared panel, so they agree with
  each other -- which is why this stayed invisible.
- **sec 2.4's "intended direction" was already the code.** `clock_needs_full_panel()` and
  `check_DNAm_extra()` no longer exist; `resolve_DNAm_extra()` is the single hard-coded
  `"Zhang2019"` case the plan proposed. Forward-looking text left in after the work landed.
- `score_systemsage()` -> `score_systemsage_group()`; the `sample_coverage(x)` method was never
  built and the matrix is reached as `$coverage$sample_miss`.
- sec 4.2 now points at `check_row_coverage()` as the consumer of the tier-2 matrix -- the reason
  no branch takes `min_samples_coverage`.

**Not reconciled, flagged instead.** sec 1.1 "Exported surface" lists `list_clocks()`,
`get_clock()`, `get_clock_probes()`, `summary()`, `augment()`; none exist and `NAMESPACE` exports
none of them (only `rbind.mc_result` + `print.mc_sim`). `get_clock` / `get_clock_probes` have built
equivalents under the `clock_*` names (`clock_entry()`, `clock_scoring_cpgs()` +
`clock_impute_ref()`); `list_clocks()` has no equivalent at all. Left as a design target rather than
silently renamed, since deciding whether that surface is still wanted is a human call. `CLAUDE.md`
names `get_clock` in its accessors invariant and lists the same unbuilt `mc_result` verbs.

---

## 2026-07-23 -- coverage gates split by axis: columns stop, rows warn and never NA

**Decision.** Three changes, one idea.

**1. `min_coverage` becomes `min_clocks_coverage` + `min_samples_coverage`** (both 0.75, so default
behavior is unchanged; no deprecation shim, pre-alpha). The 2026-07-13 entry below records the
reuse -- "the gates reuse the existing `min_coverage` argument ... so no new per-clock argument was
added". That was right when the front end had no gate of its own, but the graded `check_coverage()`
front end (2026-07-21) made one number drive two genuinely different tests, and the reuse became a
trap: a caller lowering the floor to admit an under-covered *panel* also silently relaxed the
*sample* check. The old `min_coverage = 0` in the parity tier was doing exactly that -- one flag
disabling two gates, only one of which parity meant to disable.

**2. The two axes get different severities.** Columns **stop**, rows **warn**. A missing column is
missing for every sample alike -- it is a property of the input matrix, no sample escapes it, and
scoring on it produces a uniformly wrong number, so the call stops. A sparse row is one bad sample
among good ones. The rest of the cohort is fine and the caller may well still want the number, so
the score stands and the warning names the sample. This is why there is no "row stop" band and no
graded warn band on the row side: the row check has exactly one thing to say.

**3. Dunedin stops NA-ing under-covered rows.** `score_Dunedin()`'s `low_sample` gate -- which
blanked a sample's score and warned -- is deleted along with its floor argument. It was the only
branch that did this, so the same cohort got different treatment depending on which clock was
requested: `DunedinPoAm38` returned `NA` for a sparse sample while `Hannum` returned a
cohort-mean-imputed number from the same row. Consistency wins, and it is the better default --
the imputation machinery exists precisely so a partially observed sample still scores.

**Why the check lives in `calc_clocks()`, not in `count_sample_miss()`.** The obvious home for a
per-sample warning is the hoisted helper that already computes per-sample miss counts. It does not
work: `count_sample_miss()` sees only the cached columns, not `score_needed` / `score_absent`, so
it cannot form a *fraction* -- it would need four more arguments threaded through ten call sites,
several of which are per-*component* (GrimAge, DNAmFitAge) rather than per-clock, and scoring
`"all"` would emit ~85 separate warnings. Instead `check_row_coverage()` runs once after scoring
over the record every branch already returns (`$coverage` + `$sample_miss`), and derives coverage
as `(score_present - sample_miss) / score_needed`. One warning, listing every affected clock,
capped like `check_coverage()`'s. A new branch is covered the moment it returns the standard
record -- there is nothing to remember to call. The denominator switches to the `norm_*` pair when
a clock has a normalization panel, because that is the panel `count_sample_miss()` counted over
for DunedinPACE; mixing them would report a coverage above 1 or below 0.

**Validation is up front.** `check_coverage()` asserts `min_clocks_coverage` and already runs before
scoring, but `check_row_coverage()` runs after, so `calc_clocks()` asserts `min_samples_coverage`
early. A bad argument should not cost a full scoring pass first.

---

## 2026-07-23 -- the cache consent gate fails closed: strict `ask`, path-only `assets`, fs for file work

**Decision.** Two fail-open holes in `R/mc_data.R` are closed, and the file now does its path/size/
move/delete work through `fs`.

**1. `ask` is a strict flag.** `mc_consent()` and `mc_consent_delete()` tested `!isTRUE(ask)` and
treated *everything that was not exactly TRUE* as consent -- `NA`, `NULL`, `"yes"`, `c(TRUE, TRUE)`
all skipped the prompt. Verified: `clear_mc_cache(ask = NA)` deleted a staged pack silently. Only
`FALSE` is the consent signal, so the three entry points (`load_mc_assets()`, `mc_data_download()`,
`clear_mc_cache()`) now `checkmate::assert_flag(ask)` and the helpers branch on the validated
value. An `NA` propagated from a config is a bug, not permission to delete a user's cache.

**2. `assets` must be a path in the cache verbs.** `mc_cache_dir()` accepted only a single non-empty
string and *silently fell back to the default cache dir* for anything else. But `assets` is
documented (sec 9.4) as accepting a loaded pack or list of packs, and `mc_data_download()`,
`mc_cached_files()` and `clear_mc_cache()` passed the user's raw value straight in -- only
`load_mc_assets()` canonicalized first. So `clear_mc_cache(assets = my_pack)` resolved to
`R_user_dir()` and deleted the real cache. Found by an audit probe that did exactly that (38 MB,
all three packs; re-downloaded, content-addressed so no loss). `mc_cache_dir()` now stops on a
non-NULL non-path `assets` and names the mistake. A loaded pack has no directory, so guessing one
is never right.

**Why a stop and not a coercion.** The two arguments are the only things standing between a scripted
call and a deleted shared cache, and both failures were *silent widenings* of what counts as
permission. A wrong `ask` or a wrong `assets` should cost an error message, never data.

**3. `fs` moves from `data-raw/` into `R/`** (it was already in Imports while only sync.R used it).
It buys four things that are all defect fixes rather than taste: `fs_bytes()` is vectorized and
zero-length safe, where the old `mc_bytes()` (`format.object_size`) errored with "argument is of
length zero" on an empty or absent `size_bytes`; `file_move()`/`file_delete()` **error** where
`file.rename()`/`unlink()` return FALSE unchecked, which is the pack-lands-nowhere case;
`path()`/`path_expand()` tidy separators, so prompts stop printing
`C:\Users\x\AppData\Local/R/cache/...`; and `file_exists()`/`file_size()` are vectorized.

**Also in this pass.**

- **`download.file()`'s status code is checked.** It reports some failures only through its
  (invisible) return value rather than by signalling, so the `tryCatch(error =)` wrapper alone could
  let a failed transfer reach the qs2 read and surface as a confusing checksum error with no URL.
- **The fetched pack is no longer decompressed twice.** `mc_fetch()` returns
  `list(path =, payload =)` -- the payload from the validating read it already performs -- and
  `load_mc_assets()` reads from disk only for the groups it did not just fetch. The validating read
  stays: it is what keeps a truncated transfer from landing in the cache.
- **Consent prompts are tabular** (aligned group + size, total only when there is more than one
  pack) and the delete prompt labels rows by **group id**, not the 78-character content-addressed
  file name. `mc_cached_files()` names its result by group id to carry those labels.
- **`mc_set_cache_dir()` never existed**, though sec 1.1 and 9.4 both named it. Plans corrected to
  `options(mc.cache_dir = )` per code-is-truth; no function added, since the option and
  `MC_CACHE_DIR` already cover it.

**Tests.** `test-mc_data.R` gains a full round trip (download via `file://` -> load in memory ->
`clear_mc_cache()` -> assert gone, and a following load refuses because it would have to
re-download) plus a regression test per hole: a pack passed as `assets` errors and leaves the staged
file untouched, and every non-flag `ask` errors without fetching or deleting. Both regression tests
pin `mc.cache_dir` at a tempdir so a future regression cannot reach a real cache.

---

## 2026-07-23 -- score_type() stops on an unroutable clock; the "unsupported" tag is gone

**Decision.** `score_type()` no longer returns `"unsupported"` as a thirteenth tag. A catalog entry
that no branch claims now hits `unroutable()` and stops, naming the clock's group /
`weights_format` / `computation_type`. The tag is deleted everywhere, including the three test
sites that filtered on it.

**Why.** `"unsupported"` was a sentinel that travelled a long way before anything went wrong -- the
value flowed out of the classifier, through `pack_groups_needed()`, into the dispatch `switch`, and
only errored there. Worse, two test tiers *filtered on it*: `test-sim-smoke.R` and the parity
target list both dropped unsupported clocks, so a sync that introduced an unroutable entry would
have silently **shrunk** the smoke tier rather than failing it. A scoring gap is a gap to look at,
not one to score around.

**The filters are now assertions.** Removing them costs nothing today (all 129 catalog entries
route -- verified, so `score_type()` is total over the whole catalog, not merely over `mc_index`'s
callable pool) and means the next unroutable entry fails the always-on tier at pool construction.
`test-score-reduction.R`'s "known tag" test gets stronger for free: mapping the catalog through
`score_type()` now proves routability as well as tag membership.

**Considered and rejected: keeping the classifier total because tests filter on it.** That is the
scaffolding dictating the design. The tests existed to tolerate a gap, and tolerating it is the
behavior being removed.

**The dispatch `switch` default still stops**, but now means something different -- a tag added to
`score_type()` with no branch below it (developer error), not a catalog gap. `switch()` with no
default returns `NULL` silently, so the guard stays regardless of being unreachable by
construction.

---

## 2026-07-23 -- clear_mc_cache() actually deletes, under the download consent gate

**Decision.** `clear_mc_cache()` is no longer a report-only stub. It gains `ask = TRUE` and
mirrors `mc_consent()` exactly: interactive prompt, hard refusal non-interactively, `ask = FALSE`
as the explicit-consent signal. This reverses the 2026-07-21 "kept as a must-have stub, delete flow
deferred" call.

**Why.** CRAN requires that a package which writes to `R_user_dir(.., "cache")` give the user a
supported way to reclaim it; "remove these files manually" is not one, and the stub was the only
thing standing between us and that requirement. Reusing the download consent shape means there is
one rule to learn for both directions of the cache, and the non-interactive refusal keeps a
scripted `clear_mc_cache()` from silently nuking a shared cache.

**Declining is a no-op, not an error** -- unlike a declined download, nothing downstream depends on
the deletion, so there is no failure to propagate. Deletion is confined to registry-named files in
the resolved cache dir (via `mc_cached_files()`), and a file that survives `unlink()` is an error
rather than a silent partial clear.

**Also:** `load_mc_assets()` now resolves `groups = "all"`, which previously died in `mc_asset("all")`
while every sibling verb accepted it. An empty `groups` still returns an empty registry -- it comes
from `pack_groups_needed()` and means "no external pack needed", so it must not collapse into "all"
the way `mc_resolve_groups()`'s does.

---

## 2026-07-23 -- group-specific scorer tags/functions named by group_id

**Decision.** Every group-specific `score_type()` tag and its scorer function is now named for the
clock's `group_id` rather than an ad-hoc lowercase word: `zhang`/`score_zhang` -> `Zhang2019`/
`score_Zhang2019`, `dunedin` -> `Dunedin`, `grimage` -> `GrimAge`, `physage` -> `PhysAge`,
`epitoc2` -> `EpiTOC2`, `miage` -> `MiAge`, `fitage_composite` -> `DNAmFitAge`. The shared/mechanism
tags keep descriptive names -- `linear` (shared engine, ~83 clocks), `sex_routed`, `pack_linear`
(spans PCClocks + PCBrainAge), `pack_systemsage`, `unsupported`. The dev tracker
(`dev/clock_tracker.csv`) keys scorers on `group_id`, so the code now matches that scheme.

**Why + dead code removed.** Dunedin and Zhang2019 are `cpg_coefficient`/`linear`, colliding with the
generic linear rule, so each had its own early-return `if` grafted above the group `switch()` -- a
second, parallel dispatch path. Folding both into the one group switch (which now resolves *before*
the generic linear fallback) removes that duplicate structure: there is a single group switch, then
the `cpg_coefficient` linear fallback, then `unsupported`. Naming the tag == `group_id` also means
the `calc_clocks()` dispatch reads `Zhang2019 = score_Zhang2019(...)` etc., so tag, function, and
catalog group are one word. Routing is unchanged (verified: same tag for all 129 clocks, full suite
green); this is a naming + structure cleanup only. Accessor/recipe vocabulary (`grimage_cox_coef`,
`fitage_kdm_params`, `epitoc2_params`, the `sample_scale` recipe op) is untouched -- those name
recipe steps, not dispatch branches.

## 2026-07-23 -- coverage gate splits: scoring panel stops, normalization panel warns only

**Decision.** `check_coverage()` no longer takes the worst of a clock's {scoring, normalization}
panels. The **scoring panel is the sole stop/warn gate**, applied identically to every clock (stop
under `min_coverage` or at 0 observed CpGs, warn within 10% of the floor). A **normalization
background is a separate, warn-only check** -- it never stops. The Dunedin per-sample NA gate
(`score_Dunedin()`) likewise now measures every sample against the **scoring** panel, not the
norm panel.

**Why.** Only **DunedinPACE** actually carries a normalization CpG panel (the 20,000-CpG QN gold
background). The other norm schemes in the catalog -- BMIQ (the Horvath-online family), noob -- are
scheme *annotations* with no background panel (`clock_norm_cpgs()` is empty), so they were already
gated on scoring alone. So the worst-of-both rule and the per-sample norm gate affected exactly one
clock, PACE, which meant "under-covered" was defined on a different denominator (20k) than every
other clock (its scoring panel). That was the whole inconsistency.

**Why norm is warn-only, not a second stop.** A sparse QN background degrades gracefully: absent
gold CpGs fill to their gold mean before quantile-normalizing, so PACE still returns a score. That
is a different failure mode from missing the CpGs the linear sum actually reads, which is what the
scoring gate protects. Folding both into one stop under one threshold conflated "cannot score" with
"normalized against a thin background". The coverage record already reports `norm_present /
norm_needed` separately, so nothing is lost by making the check mirror that split.

**Refines** the 2026-07-22 graded-gate entry, whose gate "gains the norm panel (and an
unconditional zero-coverage clause)". The zero-coverage clause and the graded stop/warn bands are
unchanged; only the norm panel moves from the stop gate to its own warn. The measured stop/warn
tables in that entry are unaffected except for PACE (its norm panel can no longer stop the call).

---

## 2026-07-23 -- Zhang2019 gets a sample_scale branch; it was silently mis-scored as plain linear

**Decision.** `Zhang2019` now routes to its own `score_zhang()` (`R/score_zhang.R`) via a group check
in `score_type()` placed *before* the cpg_coefficient/linear predicate (same pattern as Dunedin).
The scorer computes a per-sample z-score (`sample_scale`) and then the EN linear sum.

**Why it needed a branch -- and why this is a correctness fix, not a cleanup.** Zhang's catalog tags
are `cpg_coefficient` + `linear`, so the general predicate claimed it and it ran through
`linear_score()` as `b0 + sum(coef*beta)` on raw betas -- with **no z-score at all**. Its recipe
carries a `sample_scale` op the linear engine does not execute, and the miss was masked because the
clock sat in `KNOWN_PARITY_GAPS` (skipped, never numerically checked). So the shipped scores were
wrong, not merely imprecise.

**The scoring contract (author pred.R v1.0.0).** `sample_scale` is a within-sample z-score over
**every probe in the input matrix**, then the linear sum over the EN probes present after scaling
(policy `omit`). The moments MUST span the full array, not the 514-CpG subset -- subsetting first
changes each sample's mean/sd (the methylCIPHER v1 `calcZhang2019` bug). `calc_clocks()` already
hands each scorer the full user `DNAm`, so real usage gets full-panel moments; passing only the
subset yields approximate scores (the pre-existing `resolve_DNAm_extra()` message already says so).
The scorer uses the algebraic form `b0 + (sum_j coef_j*beta_ij - m_i*sum_j coef_j) / s_i`, so only
the present-CpG linear term is materialized; `m_i`/`s_i` come from `matrixStats::rowMeans2`/`rowSds`
(na.rm, matching base `scale`). **New dependency:** `matrixStats` (Imports) for the row moments.

**Parity.** The fixture oracle ran on the full EPIC panel, so the parity tier now loads the whole
array for Zhang (`FULL_PANEL_CLOCKS` + `cohort_betas_full()`), and Zhang leaves `KNOWN_PARITY_GAPS`.
Exact policy passes at `max_abs_diff` ~9.6e-13 vs the author enpred.

---

## 2026-07-23 -- CellDRIFT is scored via the shared linear engine, not a named branch

**Decision.** `CellDRIFT` is now supported. `score_type()` gains a `CellDRIFT` group branch that
maps `computation_type = "reference_code_required"` to the `"linear"` tag, so it runs through
`linear_score()` -- no new scorer file. With CellDRIFT and MiAge both routed, every catalog clock
now maps to a scoring branch.

**Why linear and not its own branch like MiAge / EpiTOC2.** `reference_code_required` marks how the
weights were *derived* (author code, not a self-contained recipe), not that the math is non-linear.
CellDRIFT's author code is a dual-module PCA + elastic net, but the maintainer collapsed it
**upstream** to one per-CpG coefficient vector plus an intercept (`weights_format =
cpg_coefficient`, `output_transform = identity`, `normalization = none`, `intercept =
11.2177...`). Post-collapse the score is exactly `intercept + sum(coef * beta)` over 2322 CpGs with
`vendor_mean` fill for absent probes -- identical in shape to any linear clock, so the correct
branch is the existing engine. This is the mirror image of the two clocks sharing the same tag:
MiAge (`weights_format = custom`) is a genuine per-sample L-BFGS-B optimization that cannot flatten,
and EpiTOC2 folds a ground-state offset into a mean -- both earned bespoke scorers; CellDRIFT does
not. The tag routes on the `(weights_format, computation_type)` pair per clock, and a
`reference_code_required` clock whose weights are `cpg_coefficient` is linear.

**Verification.** Parity against the author fixture on `fixtures/cohort_EPIC` passes exactly
(policy `exact`, `max_abs_diff` ~5e-14, matching the meta's recorded value); 8/2322 clock CpGs
absent from the cohort are `vendor_mean`-filled from the module-mean ref, as documented. Smoke
(bundled + `score_type != "unsupported"`) and the exact-policy parity fixture both enrol it
automatically -- no test edit.

---

## 2026-07-22 -- the `[[` rule is enforced by data provenance, not by file

**Decision.** The `[[`-not-`$` invariant now covers every read of an
**upstream-owned** structure, repo-wide: the seven remaining `clock_impute(id)$policy`
sites, `clock_group_bundle(id)$group_id`, the GrimAgeV2 `comp$file` / `$intercept` /
`$covariates` component reads, all pack payload reads (`pack$cpgs`, `$organs`,
`$systems`, `$age`, `$impute`, `$coefficient_matrix`), the asset/provenance rows in
`R/mc_data.R`, the `mc_index` data-frame columns, and ~230 sites in `data-raw/sync.R`,
which parses the meta JSON directly. Package-internal lists keep `$`.

**Why the boundary is provenance and not "all of `R/`".** The bug this prevents is an
upstream rename silently resolving to a longer sibling field. That can only happen where
the field names are set by a repo we do not control. `cpgs$score_present`,
`lp$cpg_contrib`, `design$sample_miss` and friends are containers this package builds and
consumes a few lines later; renaming one is a local refactor a test catches immediately,
and converting them buys nothing while adding noise to every scorer. Environments are
excluded outright -- `$` does not partial-match on an environment.

**Why `sync.R` mattered most.** It is the only code that reads raw upstream JSON, so it
is where a contract rename lands first, and it is *not* covered by the test suite. Its
sweep was verified three ways: the in-memory build diff (all 129 scoring panels, all 37
group tensor sets, `mc_index` and `mc_groups` bit-identical to the committed
`R/sysdata.rda`), an old-vs-new build diff (whole-object identical), and a static proof
that the file's entire diff reduces to the `$` <-> `[["..."]]` rewrite -- which also
covers the pack-encoding, lockfile and upload paths the dry-run never calls.

**Not currently a live bug anywhere.** Every field set was checked for prefix collisions
and none exists today (`imputation` is `{policy, ref}`; components are
`{name, file, row_key, col_key, intercept, covariates}`; packs are
`{cpgs, coefficient_matrix, organs, systems, age, impute, tensors}`). This is a
prospective guard: `covariates` on a component is one upstream `covariates_required`
away from replaying the exact FitAge failure.

---

## 2026-07-22 -- sex-routed aliases are minted package-side; coverage never lands on one

**Decision.** The un-suffixed stems (`DNAmFitAge`, `DNAmGrip_wAge`, ... -- 7 of them) come back as
**alias clocks**, minted by `sync.R` from `_group.meta.json` `routing.sex`. An alias owns no
weights and no panel; `score_sex_routed()` gives each sample the score of the member matching its
`Female`. The 14 sex-resolved members leave the **callable pool**: asking for one is a hard error
that names its alias. They are still scored, still columns in the result, and still where coverage
lives.

**Why mint them here and not upstream.** Upstream deliberately *deleted* the pooled `clock_id`s;
asking for them back re-opens a split they made for good reasons (near-disjoint coefficient sets,
separate intercepts, separate median tables, separate KD constants). Routing is a consumer concern
and they said so in the routing note. This is the second member of the pattern `CUSTOM_GROUPS`
established for MiAge: a small closed registry in sync that adapts an upstream contract to what a
caller needs, with nothing downstream special-cased.

**Why the members are not callable.** Scoring `DNAmGrip_wAge_Female` on a male sample returns a
plausible number, not an error -- the female model evaluates fine on anyone's betas. A pool that
allows it is a silent-wrong-answer surface. The refusal, the callable pool and the suggested
alias all derive from one function over the routing tables (`sex_routed_members()`), so they
cannot drift apart; the planned "available options" helper and its Levenshtein suggestions should
read the same source, which also lets a rejected member answer with its alias rather than a
distance guess.

**Why members are blanked outside their sex rather than scored for everyone.** Scoring both
members on all rows and selecting at the end is simpler, but it moves the footgun from the request
to the output: the dependency columns would then hold plausible wrong-sex numbers. `mask_routed_members()`
runs after all scoring -- so the alias and any downstream composite still read full columns -- and
then `NA`s the rows a member does not apply to, recomputing its `score_imputed_partial` over what
is left. The wasted compute is a few hundred CpGs.

**Coverage.** The rule is narrower than "composites do not report coverage", which is where this
started and which would have wrongly stripped `GrimAgeV1`. A score-assembled clock keeps coverage
**iff every component contributes to every sample**: true for `GrimAgeV1` and for
`DNAmFitAge_{Sex}`, false for an alias, whose two members cover disjoint halves of the cohort.
So the alias reports `per_clock[[id]] = NULL` and an all-`NA` tier-2 column. Separately, the
composites' panel was wrong for a second reason -- the resolver reached the group's shared bare
CpG list before recursing into dependencies, handing them the 627-probe family-wide `data_prep2`
prep panel. Dependencies now outrank the shared list, giving 172 / 190, which match upstream's
independently derived counts exactly. In-group inputs only: `GrimAgeV1` is excluded because it is
its own column with its own coverage row, and its 1030 probes would both double-count and swamp
the family's own figure.

---

## 2026-07-22 -- probe_sets are kept for every weights_format; the Dunedin tensor grep is gone

**Decision.** `prune_clock_meta()` no longer drops declared `probe_sets` for clocks that are not
`external_package`, and `materialize_probe_set()` carries the `file` pointer through alongside the
resolved CpGs. `dunedin_gold_means()` then reads the declared
`quantile_normalization_background` probe_set and errors when it is absent or ambiguous.

**Why.** The accessor had been locating its QN target tensor by `grep`-ing bundle tensor names for
`(^|/)<id>/gold_standard_means`. That was compensating for a regression, not a missing
declaration: DunedinPACE declares the pointer properly, but it migrated off `external_package` to
`cpg_coefficient` upstream (`8fed519`) and the prune rule silently discarded the field on the way
in. A regex over tensor names is the wrong instrument regardless -- it is anchored on an id
substring, so it would happily match a sibling clock's tensor and return the wrong QN target
without failing. An accessor that cannot find its declaration has done its job by stopping.
`clock_norm_cpgs()` also loses its `dunedin_gold_means()` fallback, which existed only to route
around the same gap.

**Blast radius, measured before syncing.** The build was dry-run against the staged clone first:
exactly one `quantile_normalization_background` probe_set appears across all 122 clocks, and no
other clock's scoring panel changed.

---

## 2026-07-22 -- the sex-split branch is deleted, not ported; accessors read with `[[` only

**Decision.** Upstream split the seven DNAmFitAge members into fourteen sex-resolved `clock_id`s
(meta `eb996f9`), retiring the `linear_sex` recipe op, the `{female, male}` imputation ref and the
`sex_stratified` / `sex_params` blocks. Downstream we **delete** the corresponding branch rather
than port it: `score_fitage_member()` and the five `fitage_*` accessors that fed it are gone, the
twelve biomarker members route through the shared `linear_score()`, and `score_fitage_composite()`
loses its sex loop entirely (one KDM table per sex, rows keyed by input `clock_id`). The closed
branch set in the detail plan loses its "Sex-split" row.

**Why deleting beats porting.** The branch existed only to pick one of two coefficient vectors
from `pheno$Female`. With sex in the key there is nothing to pick: each member is one tensor, one
intercept, one string `ref` -- exactly the shape `linear_score()` already handles. Porting would
have kept a named branch whose whole body was a `switch` on a field the catalog no longer has.
Net effect is 129 lines of scorer and 85 lines of accessor removed for no behavior change, and
`Female` now enters the family through exactly one door (GrimAgeV1's `stack` covariate).

**The trap this exposed, and why the fix is repo-wide.** The new member entries carry their `Age`
coefficient on the recipe step, not at top level. `clock_covariates_coefs()` read
`entry$covariates`, and `$` on a list **partial-matches** -- it silently resolved to
`covariates_required`, the unnamed string `"Age"`, which `covariate_coefs_from()` maps to
`numeric(0)`. Members would have scored without their Age term (-0.222 on `DNAmGrip_wAge_Female`)
and never errored: a wrong number, not a crash, and no fixture would have caught it on a cohort
where the term is small. So every catalog read in `R/accessors.R` now uses `[[`, not just the one
that bit. `options(warnPartialMatchDollar)` was rejected -- it is a session-global the package
cannot set for its users, it does not fire under `R CMD check`, and it warns rather than fails.
`[[` makes the class of bug unrepresentable instead of detectable.

**Consequence for parity.** The family goes from wholly skipped to ten of fourteen exact. The four
Grip clocks stay skip-listed for two distinct, now-separated reasons: `Grip_noAge_Female` (4e-6)
and `Grip_wAge_Male` (3e-6) are float smear against a 1e-6 tolerance, the same open
exact-tolerance policy question as DNAmADM / DNAmPACKYRS / GrimAgeV2; `Grip_noAge_Male` (~9) and
`Grip_wAge_Female` (~0.03) are the only two members with any cohort-absent CpG (6 and 3), where
the `horvath_online` oracle zero-filled and our contract fills from that sex's medians. Upstream
measured and recorded the same ~9 independently. Both are oracle-vs-contract gaps, not scorer
defects, and they are listed separately so the first class can be closed by a tolerance decision
without touching the second.

---

## 2026-07-22 -- MiAge: a reimplemented optimizer branch, and one closed hook for custom payloads

**Decision.** `MiAge` scores through `score_miage()`, a from-scratch reimplementation of the author
`mitotic.age()` search: per sample, `n = argmin` over `[10, 10000]` of
`sum_i (c_i + b_i^(n-1) * d_i - beta_i)^2`, minimized by bounded L-BFGS-B (`factr = 1`) from the
author's five starts, best fit wins, earliest start breaks ties. It is a named branch, not an
interpreter, and it leaves `CellDRIFT` as the only remaining `reference_code_required` clock
without one (it has coefficients but its own non-linear reference code, so it is its own job).

**Why a reimplementation and not the vendored author code.** `weights/MiAge/score_miage.R` is
sourced-and-`sys.source`-ing frozen author code against an `.Rdata` blob. Shipping that would put
`load()` + `sys.source()` of third-party code on the scoring path, defeat the double-precision
tensor contract, and make the branch unprofileable and untestable at the altitude every other
scorer is tested at. The maths is 8 lines; the wrapper around it was the only hard part.

**Why it is safe to call it exact.** Parity through `calc_clocks()` on the EPIC cohort is
**bit-identical** to the author fixture (`max_abs_diff = 0`, not merely under tolerance). That is
not luck: the residual sum is order-sensitive in floating point, so the branch keeps the *panel*
order (the `Additional_File1` order the author's `match()` produces) through the DNAm subset rather
than the engine's usual `c(cached, raw)` order. Verified separately that present-only summation and
the author's NA-padded 268-row summation give identical bits, which is what makes dropping absent
CpGs (policy `omit`) equivalent to the author's `na.rm = TRUE` rather than merely close.

**Why partial NA still goes to the cohort cache.** The author's `na.rm = TRUE` omits an NA *cell*
from one sample's objective; the package imputes it from the cohort mean like every other clock.
These differ in principle, and the EPIC cohort cannot tell them apart (240 present CpGs, zero NA
cells). Imputation living in one place beats per-branch fidelity to an author's incidental NA
handling -- the invariant is package-wide, and a divergence here would be invisible to the fixture
either way.

**Why `sync.R` grew `CUSTOM_GROUPS` instead of asking upstream for a CSV tensor.** The meta repo
deliberately classifies MiAge as `weights_format = "custom"`: the frozen code *is* the declared
definition, so nothing in the meta JSON references the parameter blob and the generic `weights/`
ref crawl bundled only the R source text -- a catalog entry with no panel and no coefficients. The
alternative was a meta change emitting `site_parameters.csv.gz` plus a `components` entry (the
EpiTOC2 shape). That is a cleaner contract and remains open, but it re-litigates upstream's format
call for one clock and blocks the package on a second repo's extraction/ledger cycle. The registry
is deliberately closed and tiny -- a loader plus a component declaration per custom group, applied
before the bundle crawl -- so everything downstream (probe sets, accessors, coverage gate, sim,
parity) sees an ordinary cpg-keyed tensor and knows nothing about MiAge. If upstream later ships
the tensor, deleting the registry entry is the whole migration.

**Consequence for the coverage gate.** MiAge now has a 268-CpG panel, so it is gated like every
other clock and no longer relies on the "clock with no panel, fails later on the unimplemented
scorer" escape hatch named in the graded-coverage-floor entry below. No clock in the catalog is
exempt for having an empty panel today.

---

## 2026-07-22 -- package renamed to methylCIPHERv2; v1 repo severed without a redirect

**Decision.** The package is `methylCIPHERv2`, published at `hhp94/methylCIPHERv2`. The old
`hhp94/methylCIPHER` is retired to private (deletion later), **not renamed**.

**Why not a GitHub rename.** A rename leaves a permanent redirect on the git URL, the API, and
release-asset downloads. The requirement here is that the old coordinate stop resolving, so a
redirect defeats the purpose. Private gives the same public 404 as deletion while keeping the
existing releases and history recoverable during the transition.

**Why a new name at all.** Merging the rewrite back into v1 under the v1 name would break every
existing caller's API. A distinct name lets v1 keep its contract and v2 keep its own.

**Why `mc_` and not `methylCIPHERv2_` for internals.** The S3 classes are now `mc_result` /
`mc_sim`, options `mc.cache_dir` / `mc.release_repo`, env vars `MC_*`. `mc_` was already the
pervasive internal prefix (`mc_catalog`, `mc_cache_dir`, `mc_asset`, ...), and it does not have to
be renamed again at v3. The R constant became `MC_DEFAULT_RELEASE_REPO` so it does not read as a
duplicate of the new `MC_RELEASE_REPO` env var in `sync.R`.

**Exception: `R_user_dir` keeps the package name.** `mc_default_cache_dir()` uses
`tools::R_user_dir("methylCIPHERv2", "cache")` -- CRAN expects the real package name there, not the
internal shorthand. This orphans any existing cache; one re-download.

**No `sysdata.rda` rebuild.** `mc_provenance$external_assets` holds only `release_tag` + `file`;
the repo is assembled at runtime in `mc_asset_url()`. Only the constant moved, so the compiled
catalog is untouched by the rename.

**History kept.** The 184-commit history (inherited tree plus the rewrite) is pushed as-is. It is
the record of which commits are the rewrite.

**Out of scope, still open.** Licensing. `LICENSE` still names Yale University as copyright holder
and `DESCRIPTION` still carries `License: BSD_3_clause` and the original `Authors@R`, none of which
were touched here. Separately, the catalog's per-clock `license` field -- read by no R code -- spans
14 restrictive clocks (12 `non-commercial`, `Mayne` CC BY-NC-ND, `DunedinPACE` research-use-only),
8 copyleft (GPL-2/2+/3), and 54 with no clear grant. That is a CRAN blocker to resolve before alpha.

**Entries below this line predate the rename** and refer to the old names (`methylCIPHER.cache_dir`,
`METHYLCIPHER_PARITY`, class `methylCIPHER`, `R/methylCIPHER_data.R`). This log is append-only, so
they were left as written.

---

## 2026-07-22 -- the coverage floor becomes a graded stop/warn gate at 0.75

**Decision.** `warn_low_coverage()` becomes `check_coverage()`, a two-band gate on a clock's worst
panel: **stop** under `min_coverage` (default lowered 0.8 -> **0.75**) or at 0 observed CpGs,
**warn** when it clears the floor by under 10% (`< 1.1 x min_coverage`, capped at 1), silent above.
Same call site in `calc_clocks()` -- after `resolve_cpgs()`, before `build_partial_cache()` and the
scoring loop, so nothing expensive runs first. It gains the norm panel and an unconditional
zero-coverage clause.

**Why graded.** A single boundary made one number carry two jobs: "this cannot be scored" and "this
is thinner than you probably want". Measured on EPICv2 (below), a hard floor at 0.8 stopped 12
clocks, of which 8 were within 5 points of the line -- usable scores killed by a threshold, which
pushes callers to lower `min_coverage` globally and lose the signal entirely. Two bands let the
default stay honest without being the only lever.

**Why the band is relative, not a second absolute number.** It tracks a caller who moves the floor
(0.6 -> warn to 0.66) instead of needing its own knob, and it keeps the API at one argument. Full
coverage never warns whatever the threshold, so a strict `min_coverage = 1` stays silent on
complete panels rather than warning about everything.

**Why 0.75 and not 0.6.** 0.6 would score `DunedinPoAm38` (0.67 on EPICv2) and the wrong-array
`RetroAge450K` (0.64) by default, which reads as the package endorsing those numbers. At 0.75 a
caller who wants them has to pass a lower `min_coverage` explicitly and own the call.

**Why.** Warning-and-proceeding was emitting scores that look publishable and are not. Measured on
a foreign panel at zero coverage: `Hannum` returned `0`, `EpiTOC` returned `NaN`, and `Horvath1`
returned **42.4, 39.1, 41.3, 42.3 years** -- `anti_trafo(intercept)`, a fabricated age per sample,
from no observed CpGs. Three different degenerate behaviors for one degenerate input, and the
default `min_coverage = 0.8` did not stop any of them. A number a user can paste into a paper is
not an acceptable output for "we observed nothing".

**Why the zero clause is unconditional.** For any `min_coverage > 0` the ratio test already
subsumes it (0 < threshold). The gap was `min_coverage = 0`, where the old function early-returned
without checking anything -- and that is exactly the setting someone reaches for when coverage is
already bad. It costs one `||` in the existing predicate, not a second mechanism. Clocks with an
empty panel (`MiAge`) are exempt so they keep failing on the clearer unimplemented-scorer error.

**Why not auto-drop the failing clocks (`setdiff`) and warn.** Mass low coverage is near-always one
root cause -- transposed matrix, 450k data against EPIC panels -- not N independent problems.
Auto-dropping turns one root cause into N warnings, or worse, an empty result that looks like a
successful run. Same reason the error caps its listing at 10 clocks.

**Why one global scalar and not per-clock thresholds.** The gate is monotonic: lowering the
threshold can never make a passing clock fail. So a caller who wants 0.6 for one clock can set 0.6
for the call and lose nothing -- the other clocks already cleared a stricter bar, and their actual
coverage is still reported per clock in `$coverage`.

**Consequences.** `score_dunedin()`'s whole-clock gate (`clock_gated`) is now unreachable and is
deleted, along with its warning -- `model_cov` was byte-for-byte the ratio the front end computes.
Its **per-sample** gate (`low_sample`) stays: that reads per-row NA counts, which the
sample-invariant front end cannot see. The parity tier now passes `min_coverage = 0` explicitly --
the EPIC cohort under-covers some panels (CausAge 420/585) and the author oracle ran on that same
subset, so parity gates numbers, not coverage policy.

**Measured.** Against Illumina manifests (`slideimp.extra::ilmn_manifest(chip, deduped = TRUE)`),
over the 85 bundled scoreable clocks, at `min_coverage = 0.75` (warn band 0.75-0.825):

| Array | probes | stop | warn | ok |
|---|---|---|---|---|
| 450K | 485,577 | 16 | 0 | 69 |
| EPICv1 | 865,918 | 1 | 1 | 83 |
| EPICv2 | 930,658 | 4 | 14 | 67 |
| MSA | 269,094 | 41 | 0 | 44 |

EPICv2 stops shrink from 12 (at 0.8) to **4**: `RetroAge450K` 0.64, `DunedinPoAm38` 0.67,
`DNAmFEV1_wAge` 0.69, `DNAmGrip_wAge` 0.69. The GrimAge/DNAmFitAge cluster that bunched at
0.75-0.80 now warns and scores. `RetroAge450K` stopping on EPICv2 is the gate working correctly --
it is a 450K-specific weight set with an EPICv2 sibling, so scoring it on EPICv2 data was always
wrong; the same logic explains 450K's 16 (EPIC-era clocks) and MSA's 41 (a 269k array).
EPICv1's single stop is `CausAge` 0.72, which matches the 420/585 = 0.718 seen on the EPIC cohort
fixture.

**Non-issue, checked.** The `clock_id == group_id` collision (SystemsAge, DNAmFitAge, RepliTali,
EpiTOC2 -- `resolve_clocks()` gives the group precedence, so requesting `"EpiTOC2"` also scores
HypoClock) does not strand any flagship in practice. HypoClock's 678 solo-WCGWs cover **100% of
450K, 96% EPICv1, 85% EPICv2, 96% MSA** -- it clears the floor, and the warn band, on every real
array. The collision remains a latent API gap (no token names the member clock alone), not an
active bug. The whole epiTOC family is "ok" on all four arrays; no clock there is hardcoded out,
and `dev/hardcoded.md` was not created -- it would have had no entries. If a disabled-clock list is
ever needed, the pattern to copy is `KNOWN_PARITY_GAPS` in `test-fixtures-parity.R`: a named
`clock_id -> reason` constant next to the code that reads it, which cannot drift the way a separate
doc would.

---

## 2026-07-22 -- EpiTOC2 scores the full model only, as a closed branch

**Decision.** `score_epitoc2()` (`R/score_epitoc2.R`) implements one scorer for the `EpiTOC2`
clock: `2 * mean_j((beta_ij - beta0_j) / (delta_j * (1 - beta0_j)))` over represented CpGs. It
reuses `linear_predictor()` for the CpG term (so the partial-NA cache path is shared) and folds
the ground state in as a constant offset inside the mean.

**Why not the other author outputs.** `epiTOC2.R` returns eight things. Four are already covered
or out of scope here: `pcgtAge` is the separate catalog clock `EpiTOC`, `hypoSC` is `HypoClock`
(both plain `linear_mean`, both at exact parity), and `irS`/`irS2`/`irT`/`irT2` are just
`tnsc / age` and its cohort median -- a downstream ratio, not a clock, and `irT` is cohort-level
so it has no per-sample column in an `n x k` score matrix. That leaves `tnsc` (full model, what we
ship) and `tnsc2` (the beta0 = 0 approximation). The catalog carries one recipe with
`model: "full"` and one fixture, so `tnsc2` has no clock id and no golden; adding it would mean a
new master row upstream in `methylCIPHER-meta` first. Revisit only if the approximation is
requested -- the branch already has `delta`/`beta0` in hand, so it is a one-line variant.

**Why the author's ground-state alert is not ported.** The header prose says the function warns
when many beta values fall below the ground state; the shipped code never does, and the EPIC
cohort fixture has 1148 such values while still matching at max_abs_diff 0. Porting a warning the
oracle does not emit would be noise on ordinary data.

**On `computation_type = "reference_code_required"`.** That tag means "the meta recipe does not
fully specify the math, read `code_ref`" -- it is not a licence for an interpreter. It routes to a
named branch like every other non-linear clock. Two catalog clocks still carry it and remain
unsupported: `MiAge` and `CellDRIFT`.

---

## 2026-07-22 -- one hash, not four: `file_sha256` and the runtime drift re-hash are gone

**Decision.** The external-pack layer keeps exactly two integrity mechanisms, and neither costs
anything at score time:

- **`payload_hash`** (maintainer-side, build) -> pack filename + release tag. This *is* the
  "don't re-upload unchanged weights" guard, via the content-addressed name plus `gh_upload.py`'s
  `get_release(tag)` + asset-name skip. Kept.
- **qs2's `validate_checksum`** on read -> transfer integrity and bit rot. Kept; free.

Removed:

1. **`file_sha256` entirely** -- `file_sha256_of()`, the field on the shipped registry row, and the
   `sha256` key in the `gh_upload.py` stdin manifest. It was dead in both directions: no R code read
   it out of `mc_provenance`, and `gh_upload.py` only ever used `group_id`/`tag`/`path`/`name` --
   the `sha256` it was handed was never referenced. The 2026-07-21 entry already called it a
   "candidate for later removal"; this completes that.
2. **The runtime `payload_hash` drift check** (`mc_payload_hash()` + the WARN in `mc_read_pack()`,
   both deleted; `load_mc_assets()` now calls `qs2::qs_read()` directly). This reverses the
   2026-07-21 "the recompute is cheap belt-and-braces" line, which measurement did not support:
   re-serializing and hashing the PCClocks payload cost **0.204s against 0.203s for the entire
   `qs_read`**, i.e. it doubled every load. What it bought was near-zero -- the filename is
   `<group>-<payload_hash>.qs2`, so a name match already asserts the content address, and the only
   scenario it caught was a deliberately mis-named cache file.
3. **`payload_hash` as a registry-row field** -- redundant with `file`/`release_tag`, which embed it.

**Why now.** The package is pre-alpha and moving fast; layered hashes and warn-only checks are
friction on the path that has to stay quick, and they were guarding against failure modes that
either cannot happen (content-addressed names) or are already covered (qs2 checksum). Correctness is
the parity fixtures' job, not a hash's.

**Cost.** Changes the shape of `mc_provenance$external_assets` rows. Nothing reads the dropped
fields, so the currently shipped `sysdata.rda` stays loadable until the next `sync()`.

**Also swept in the same pass** (all confirmed unreferenced across `R/` and `tests/`):
`mc_set_cache_dir()` (the `methylCIPHER.cache_dir` option is still honored by `mc_cache_dir()`),
`FIXTURE_FIELDS` in sync.R (`prune_fixture()` hand-builds its stub and never used it), and
`R/generics.R` -- whose `rbind.methylCIPHER` was an empty stub registered on class `methylCIPHER`
while records are actually classed `methylCIPHER_result`, so it could never have dispatched.

---

## 2026-07-22 -- external packs carry weights only; scorer tags name the branch they run

**Decision.** Three related cleanups, all found by tracing sync.R -> `calc_clocks()` end to end.

1. **`catalog` and `group` are dropped from the external pack payload**
   (`stable_external_payload()` in [sync.R](../data-raw/sync.R)). This completes the 2026-07-21
   "external scoring panels come from the pack, not the catalog" split, which only did half the
   job: `drop_external_probe_cpgs()` stripped the resolved panels from the shipped `mc_catalog`,
   but `bundle$catalog <- catalog$clocks[...]` put a *second*, unstripped copy back inside the
   pack. Measured on the PCClocks pack: `pack$catalog` was **75.5 MB of an 89.9 MB in-memory
   pack**, because all 14 member entries each carried a full 78,464-element `probe_sets$cpgs`
   vector that is `identical()` to the `pack$cpgs` sitting next to it. The genuine metadata is
   71.8 KB. Nothing in `R/` or `tests/` ever read `pack$catalog` or `pack$group` -- the runtime
   reads `mc_catalog`/`mc_groups`, per the "accessors are the executable schema" invariant.
   Same reasoning as the 2026-07-17 F1 entry (write-only fields leaking into `payload_hash`).

2. **`payload_hash` gets closer to "moves only when the weights move."** Dropping the embedded
   catalog removes the largest source of hash churn that is not a weight: previously any catalog
   metadata edit (a `pmid`, a fixture stub, a recipe `note`) re-minted all three release tags.
   `schema_version` and `encoding_version` remain in the payload and still churn it; that is
   intentional (a layout change *must* re-mint) but it means a schema bump alone re-uploads
   unchanged weights.

3. **Scorer tags name their branch.** `score_type()` returned `"linear"` for PCClocks/PCBrainAge
   members and `"linear"`/`"systemsage"` for SystemsAge members, but `calc_clocks()` discarded
   the value for every external clock and dispatched on `clock_is_external()` + `group_id`
   instead -- so the tag documented a branch that never ran. Tags are now `pack_linear` and
   `pack_systemsage`, `is_pack_scored()` is the single predicate behind both `pack_groups_needed()`
   and the `is_pack` split, and `score_pack_group()` switches on the tag. No behavior change:
   SystemsAge organ sub-clocks were already routed to `score_systemsage_group()` (which splits
   organs from composites internally) despite being tagged `"linear"`.

**Cost.** (1) changes `payload_hash` for all three groups, so the next `sync(upload = TRUE)`
re-mints and re-uploads them once. Worth it: it is also a large cut to the download and to
resident memory per pack.

**Also.** Renamed sync's `resolve_cpgs()` -> `pack_canonical_cpgs()` (it only exists because packs
store one CpG-name vector for a set of order-aligned numeric tensors); the name now belongs solely
to the runtime present/absent splitter in `R/resolve_inputs.R`. Deleted the empty `impute_DNAm()`
stub and renamed `R/impute_DNAm.R` -> `R/missingness.R` to match what the file actually holds
(`scan_missing_cpgs()`, `build_partial_cache()`).

---

## 2026-07-22 -- External scoring panels live in the pack, not the bundled catalog

**Decision.** `clock_scoring_cpgs(id, packs)` returns `pack$cpgs` for an external clock and never
consults the catalog; `sync()` correspondingly stops writing those CpGs into `mc_catalog`
(`drop_external_probe_cpgs()`). Panels for bundled clocks are unchanged. `clock_panels()` and
`needed_cpgs_union()` gained a `packs` argument, and `calc_clocks()` now resolves packs *before*
panels, since the panel is in the pack.

**Why -- correctness first.** The panel was stored twice: in `mc_catalog` and in the pack. Verified
`identical()` for all 28 external clocks today, but nothing enforced it. Packs are content-addressed
and independently re-uploadable, so a weights update that shifted the panel would silently disagree
with the bundled copy, and the failure mode is misaligned scoring rather than an error. One source
kills that class of bug. It also matches the existing convention -- `clock_coefs(id, packs)` and
`clock_impute_ref(id, packs)` already read external data from the pack; the panel was the outlier.

**Why it is also the big perf lever.** External panels are 3 distinct vectors but **92%** of all
panel elements (561,491 of 610,619; PCClocks 78,464 x 14 members, SystemsAge 125,175 x 13,
PCBrainAge 357,852). Dropping them takes `mc_catalog` from 1.79 MB / 1.26s to **0.09 MB / 0.021s**.
Separately, all 14 PCClocks members now receive the *same SEXP* from `pack$cpgs`, so `dedup_panels()`
short-circuits on pointer identity and there is no per-clock `unique()` over 78k strings:
`clock_panels()` 0.112s -> **0.008s**.

**Why not a keyed panel store in sysdata.** The alternative was splitting panels into a separate
lazy object (or per-group objects) so the catalog stays cheap while pack-free `sim_DNAm("PCClocks")`
keeps working. Rejected: it preserves the duplication that is the actual hazard, and buys a dev
convenience that is not worth it -- an external clock cannot be *scored* without its pack anyway, so
a pack-free simulated panel only ever exercises machinery.

**Consequence.** `sim_DNAm()` takes `assets` / `ask` and requires the pack for external clocks;
with a cached pack the old call still works. `test-score-external.R` now mints synthetic panels for
its in-memory packs instead of borrowing the real ones (faster, and it was circular once the pack
owns the panel). The code change and the data change are independent: until `sync()` is re-run the
bundled copy is simply ignored, so there is no flag day.

**Not yet landed.** The 0.09 MB / 0.021s catalog needs a `sync()` re-run to regenerate
`R/sysdata.rda`; that is maintainer-side and was not done here.

---

## 2026-07-22 -- sysdata.rda ships xz, not use_data()'s bzip2 default

**Decision.** `sync()` passes `compress = "xz"` to `usethis::use_data()`, and the committed
`R/sysdata.rda` was recompressed in place to match (content verified `identical()` per object, so
no `sync()` re-run was needed).

**Why.** `use_data(internal = TRUE)` defaults to `compress = "bzip2"`, which is the worst of both
worlds for this blob. Measured on the same five objects: bzip2 11.89 MB / 4.36s load, gzip
13.40 MB / 2.06s, **xz 2.31 MB / 1.85s**, uncompressed 55.22 MB / 1.34s. xz wins on size *and*
load; nothing trades off. Relevant to the CRAN target too -- 11.89 MB blew the 5 MB tarball
guidance on its own, 2.31 MB does not.

**Why this is the lever for the dev loop.** Every `load_all()` / `test()` / `check()` / `document()`
pays this deserialization eagerly (pkgload `load()`s `sysdata.rda` directly); `document()` in
particular is ~all sysdata, since parsing all 15 `R/*.R` is 0.015s.

**Sharing SEXPs does nothing; restructuring does.** Making the 114 stored probe-set vectors share
one SEXP per distinct panel was measured and is a no-op -- 1.79 MB either way, load 2.14s -> 2.01s
-- because R's serializer writes each vector out in full and only ref-caches the CHARSXPs.
*Structurally* splitting the panels out of `mc_catalog` into a keyed store is a different thing and
does work: catalog metadata alone is 0.01 MB / **0.009s**, the 84-panel store is 1.76 MB / **0.385s**,
versus 1.79 MB / 2.05s monolithic. The win is element count -- 3.1M nested elements collapse to
590k, so the string-cache does 5x less work. (Supersedes the first version of this entry, which
concluded from the SEXP-sharing result that catalog dedupe was a dead end. It is not.)

---

## 2026-07-22 -- CpG set math: dedupe identical panels, not special-case external groups

**Decision.** `resolve_cpgs()` runs the present/absent split **once per distinct CpG panel**, not
once per clock. `clock_panels()` fetches each clock's scoring + norm panels a single time and maps
identical vectors onto a shared set (`dedup_panels()`, linear scan on `identical()`); member clocks
then share the resulting `present`/`absent` vectors by reference. `calc_clocks()` calls
`clock_panels()` once and feeds both `panels_union()` (missingness scan) and `resolve_cpgs()`.

**Why.** Every member of an external pack group carries the *same* panel: PCClocks is 14 clocks x
one 78,464-CpG panel, SystemsAge 13 x 125,175. The old code hashed that panel once per clock, and
did the whole fetch twice -- `needed_cpgs_union()` and `resolve_cpgs()` each called
`clock_scoring_cpgs()` for all 14. Measured on `sim_DNAm("PCClocks")`: `resolve_cpgs` 0.333s ->
0.022s, whole set-math path 0.486s -> 0.140s. It also drops 14 redundant copies of a 78k character
vector out of `cpg_list`.

**Why not scope it to PCClocks/SystemsAge.** The obvious fix is a group-keyed branch for external
packs, but that buys a new special case for a property that is not actually about externality --
it is about panels being equal. Deduping on `identical()` needs no group knowledge, keeps the
closed branch set closed, and stays correct if a bundled family ever shares a panel. Bundled groups
pay only the dedupe scan, which is a no-op for them (GrimAge: 12 clocks / 11 distinct panels, all
under 1030 CpGs); `clocks = "all"` is 115 clocks / 85 distinct panels and the scan short-circuits
on length mismatch.

**Verified.** Deduped per-clock sets are `identical()` to the old per-clock math, and cohort parity
runs green (198 passed / 0 failed), so the shared vectors reach the pack scorers unchanged.

**Not addressed here.** The remaining ~1.1s of a cold `calc_clocks("PCClocks")` is `load_mc_assets()`
reading the qs2 pack from the runtime cache; passing `assets =` an already-loaded pack drops the call
to ~0.3s and steady-state scoring is ~0.09s. Pack I/O is a separate question from the set math.

---

## 2026-07-22 -- Missing Age/Female: propagate NA, warn once in check_pheno

**Decision.** An `NA` in a required covariate (`Age`, `Female`) is a **legal input**. It propagates
through the linear algebra to an `NA` score for that sample and nothing else: the row is never
dropped, no error is thrown, and every other sample scores normally. `check_pheno()` emits **one**
warning naming each affected covariate column and its NA sample count.

**Why warn upstream, not per clock.** Propagation already worked (`X %*% cox` in `score_grimage()`,
`cov_mat %*% cov_coefs` in `linear_predictor()`), but it was *silent* -- a user with a few missing
ages got an `NA` column with no pointer back to `pheno`. Warning at the single validation site is
succinct and clock-count-invariant: one message regardless of how many clocks consume the covariate,
emitted before any scoring. A per-clock or per-scorer warning would fire N times for the same root
cause and would have to be threaded through every branch.

**Why not error or drop.** Dropping rows silently changes `nrow(scores)` and breaks id alignment
against the caller's pheno; erroring makes one bad row block a whole cohort. `NA` in, `NA` out is
the least surprising contract and keeps the result rectangular.

**Scoping.** The count is taken over rows that survive the id-join (`match(sample_id, pheno[[ID]])`),
so an `NA` on a cohort row that is not being scored does not warn. Under positional alignment all
rows count, since all are scored.

**Consequence.** Surrogates inherit the `NA` only where they declare the covariate: with `Age`
missing, `DNAmPAI1` / `DNAmLeptin` / `DNAmlogCRP` (no covariates) still score finite while the
Age-consuming surrogates and both GrimAge composites go `NA`. Coverage counters are unaffected --
they track CpGs, not covariates.

**Sites.** `check_pheno()` / `warn_missing_covariates()` in `R/resolve_inputs.R`; call site in
`calc_clocks()` now passes `sample_id`. Tests in `tests/testthat/test-score-grimage.R`.

---

## 2026-07-21 -- Dunedin is a bundled QN family (not an external adapter); Hybrid missing data

**Decision (routing).** DunedinPoAm38 and DunedinPACE ship **bundled** (`cpg_coefficient` /
`linear`, group `Dunedin`), not as external `computation_type=wrapper` adapters as the early
plan assumed. Both route to a shared `score_dunedin()` branch. PoAm38 is plain linear; PACE
quantile-normalizes a 20000-CpG gold panel to `gold_standard_means` (via `betanorm::quantile_norm`,
bit-exact with the author's `preprocessCore::normalize.quantiles.use.target`) and then scores
linearly on its 173 model CpGs (a subset of the panel). The gold panel ships as the
`gold_standard_means` tensor; its names are the QN background (there is no inlined
`quantile_normalization_background` probe_set in the compiled catalog), so `clock_norm_cpgs()`
falls back to those names and `dunedin_gold_means()` reads the tensor. `present_needed_union`
now also unions `norm_present`, so partial-NA gold-panel CpGs get the shared cohort cache.

**Why bundled + betanorm.** Both clocks reproduce the author fixtures to machine epsilon on the
EPIC cohort (PoAm38 4.4e-16; PACE 1.8e-15). `betanorm` stays a **Suggests** soft dep: it is a
GitHub-only (`hhp94/betanorm`) compiled package, so a CRAN target cannot Import it; PACE errors
cleanly and sim-smoke / value goldens skip when it is absent. Reimplementing preprocessCore's QN
in pure R was rejected -- betanorm already matches it exactly and is the maintainer's own package.

**Decision (missing data -- Hybrid).** The author's `PACEProjector`/`PoAmProjector` cascade uses a
single `proportionOfProbesRequired` threshold to (a) NA a whole clock under-covered on its model
(PACE also gold) panel, (b) NA a sample observing too few panel cells, (c) impute lightly-missing
present probes to the cohort mean, (d) replace heavily-missing present probes with vendor means,
and PACE fills fully-absent gold probes with `gold_standard_means` before QN. We adopt a **Hybrid**:
keep methylCIPHER-native imputation (partial NA -> shared cohort cache; fully-absent -> vendor
means (PoAm) / gold fill (PACE)), and add only the two coarse **NA-gates** (a) and (b). We drop the
per-probe threshold split (d) and the author's quirk that leaves heavily-missing non-model gold
probes as NA. The gates reuse the existing `min_coverage` argument (default 0.8) -- it *is* the
author's proportion-of-probes-required -- so no new per-clock argument was added; `min_coverage = 0`
disables warning and gating together. The parity cohort has full coverage, so the divergent
degraded-coverage rules are exercised only by always-on value goldens, not parity.

**Decision (error handling -- deliberately not defensive).** `score_dunedin()` mirrors the author
projectors' stance: **fill or drop to NA, never stop.** Insufficient coverage yields NA (whole clock
or per sample); a missing scoring CpG is filled. The author's only `stop()` is a non-numeric-matrix
check, and ours is likewise a single one -- `betanorm` missing -- placed at the **point of use**
(inside the scoring block) so an under-covered PACE call still returns NA without needing the
package, exactly as the author reaches `preprocessCore` only after its NA-gate. Three earlier guards
were removed as dead defensiveness: "model CpGs absent from the gold panel" and "absent CpG lacks a
fill value" can never fire (model CpGs are a subset of the gold panel; `gold_standard_means` and
`model_means` are complete by construction), and a `cov_ratio()` divide-by-zero helper guarded a
panel that is never empty. Catalog-integrity assertions belong to sync, not to a scorer.

**Why not a general per-clock arg mechanism / `...`.** Only one override (the coverage threshold)
exists and it already had a home (`min_coverage`), so a `clock_args` dict or a `...`-intersection
dispatch would be speculative. A `...` that silently swallows a mistyped argument is a correctness
hazard in a scoring package (a typo'd threshold changes results with no signal) and is the same
arg-walker anti-pattern the engine forbids -- the dead `...` in `calc_clocks()` was removed. If a
general mechanism is ever needed (>= 2-3 real params), it will be an explicit `clock_args` list that
errors on an unknown clock id or arg name, never `...`.

**Parity.** PACE's catalog `fixture$parity_policy` is still `skipped` (a meta-repo value that the
portable path could not previously recompute). Now that betanorm makes it bit-exact, flipping it to
`exact` is a meta.json + re-`sync()` step (maintainer-side); until then PACE's numeric correctness
rests on the always-on value goldens (allowed while parity is skip-listed). PoAm38 (policy `exact`)
is picked up automatically by `parity_targets()` and passes.

**Also fixed here.** `test-fixtures-parity.R` wrapped its parity loop in
`if (skip_if_no_cohort() | skip_if_no_pack(clock_id))`, referencing an undefined `clock_id` at file
scope -- it errored on any enabled parity run. Moved both skips inside each `test_that` (matching the
PhysAge block), which is what the 2026-07-21 parity-gating entry below intended.

## 2026-07-21 -- External clock tests smoke-only; check_DNAm orientation keys on ^cg

**Decision (external tests).** `test-score-external.R` no longer re-derives scoring values. Every
external member (PCClocks, PCBrainAge, SystemsAge, incl. the SystemsAge systems_PCA composite) has
a passing `exact` cohort-parity fixture, so parity owns the goldens; this file keeps only always-on
coverage parity cannot give: the in-memory `assets=pack` closed-set path runs and returns the right
shape (smoke), the wrong-pack error, the PCClocks no-family-expansion subset contract, and the
vendor-fill coverage flags (`score_imputed_full`). Completes the partial removal in the "test
altitude" entry below -- the composite golden had survived in the tree; now all external in-test
goldens are gone. Trade-off: external numeric correctness now rests entirely on the (CRAN-skipped,
pack-gated) parity tier, mirroring how bundled clocks lean on sim-smoke + parity. Shared per-group
packs are built once at file scope (they were rebuilt byte-identically per test).

**Decision (check_DNAm).** The orientation guard drops the `nrow > ncol` dimensional heuristic --
it false-positived whenever samples outnumber a small CpG panel (every parity fixture) -- and keys
solely on the `^cg` prefix: warn when no column looks like a CpG, naming the transposed case when
the rows do. A genuinely transposed matrix is also caught downstream by `warn_low_coverage`.

## 2026-07-21 -- Cohort parity tier gated behind METHYLCIPHER_PARITY, not just file.exists()

**Decision.** `test-fixtures-parity.R` now skips unless `METHYLCIPHER_PARITY=1` is set (in
addition to the cohort duckdb being staged). The env flag gates both the file-scoped duckdb
connection and each test's `skip_if_no_cohort()`. A dev-only, build-ignored helper
`test_parity()` (`R/dev-utils.R`) sets the flag via `withr::with_envvar()` for one run.

**Why.** The parity tier is ~34s of a ~37s suite (one duckdb query + a real `calc_clocks()`
per fixtured clock), so on a machine with the cohort staged every `devtools::test()` paid it
on each change. Default runs are now ~4s; the science gate runs on demand via `test_parity()`
(or any env set). No CI workflow exists yet; whatever runs the gate must export the flag.
`R/dev-utils.R` is git-tracked but `.Rbuildignore`d, so it is available after `load_all()`
yet never ships to CRAN.

## 2026-07-21 -- Batched pack scorers for PCClocks / SystemsAge; no general grouping layer

**Decision.** External groups whose members share one large CpG panel (PCClocks, SystemsAge) are
now scored per *group* in a single shared `DNAm[, present]` subset + one matmul over all requested
columns (`score_pack_group()` in `R/score_pack.R`), instead of routing each member through
`linear_score()`. `resolve_clocks()` is unchanged -- no family expansion, no reverse-resolve, no
member-input restriction; requesting one member returns one member, and the batch just collapses
whatever pack members are in the plan into one matmul.

**Why.** The dominating cost of scoring a pack member is the panel-wide subset, not the matmul.
Per clock it is repeated once per member: measured **13.7x** for 14 PC clocks (78k-wide panel);
the SystemsAge composite repeats it ~24 times across its age/organs/systems sub-steps. Batching
pays the subset once. It reproduces the exact imputation contract (partial-NA cohort cache; absent
-> vendor-mean offset from `pack$impute`), so parity holds within fp tolerance -- guarded by the
in-test goldens (multi-member PC with covariates + `anti.trafo`; the SystemsAge systems_PCA
composite re-derived) and cohort parity.

**Why not a general `calc_group` layer** (fold *all* same-recipe linear clocks into one matmul).
Considered and rejected. The win only exists where the panel is huge *and shared*: the bundled
clocks' scoring CpGs overlap just **1.6x** (sum 33,249 vs union 20,724), so batching all 86 saves
~26ms on a ~60ms operation while the dense union matmul does ~50x more (zero-filled) flops -- and
it would fold per-clock vendor refs (bundled refs differ per clock, unlike the packs' single shared
`pack$impute`), mean denominators, covariates, output transforms, and coverage into a batched
engine: the recipe-ish per-clock logic the "one engine + closed branch set" invariant keeps out.
Not worth it. Batching earns its complexity only for the two shared-panel packs.

**Why GrimAge is excluded** though it is a family. Its surrogates have disjoint CpG sets (no shared
panel -> no subset win), and `"GrimAge"` as a group token already expands to its members. Member ->
family expansion was dropped entirely: it made output unwieldy (one member request -> 14 columns)
and required reverse-resolve, while buying nothing the group token doesn't already give.

## 2026-07-21 -- Test suite deliberately loosened; "test altitude" policy set

**Decision.** The suite had grown too tight for a fast-moving pre-alpha: it pinned exact
error-message strings, maintainer-side plumbing shapes (asset filenames, release tags, download
URLs, the 4-layer cache-dir precedence order), per-clock internal dispatch-tag tables
(`clock_reduction()`, `score_type()`), and re-derived full recipes in-test. All of these break on
refactors that change no observable behavior. Loosened them so tests assert on `calc_clocks()`
output and the parity science gate only; the policy now lives in CLAUDE.md "Test altitude".

- Error/warning/message assertions drop the regex -> bare `expect_error/warning/message`; the
  behavior after the throw (writes nothing, still returns the payload, etc.) is still asserted.
- Dropped the per-clock reduction/route tag tables; kept the closed-set "every clock maps to a
  known tag" guard. Routing is proven through `calc_clocks()` output.
- Dropped filename/release-tag/URL format asserts and the cache-dir precedence spec; kept the
  download *behaviors* (verify, no scratch, warn-not-stop on hash drift, closed-set, never
  auto-delete).
- Removed the SystemsAge composite recipe re-derivation -- cohort parity passes it at ~1e-11, so
  the fixture owns the golden. **Kept** the DNAmFitAge KDM re-derivation: its parity is skip-listed,
  so the in-test math is currently its only numeric gate. Revisit when that parity lands.

**Added.** `test-sim-smoke.R` -- `expect_no_error` over every bundled, supported clock via
`sim_DNAm()` + `calc_clocks()`. This is the always-on crash net the Testing section had described
but no test implemented. Net: ~256 lines of brittle assertions removed, ~84 crash-smoke cases
added; full suite 363 pass / 0 fail.

**Not done.** S3 record-verb tests (`as.matrix` / `[` / `cbind` / ...) -- the generics are not
implemented yet, so there is nothing to test; add them when the methods land.

---

## 2026-07-21 -- SystemsAge orchestrator landed (external.md step 4); external scoring complete

**Decision.** Wired the SystemsAge family (Sehgal 2024, 13 members) -- the last external.md gap.
All 13 members now pass cohort parity at exact tolerance (max_abs_diff ~1e-11, gate 1e-6), so the
three external groups (PCClocks, PCBrainAge, SystemsAge) are fully scored.

1. **Organ sub-clocks go through the shared linear engine, not the orchestrator.** The 11
   organ/system members (Blood..MusculoSkeletal) are `cpg_coefficient`/`linear` -- literally
   `intercept + sum(coef*beta)` with vendor-mean fill. `score_type()` now returns `"linear"` for
   them; `clock_coefs()` gained a one-line SystemsAge branch sourcing the column from `pack$organs`
   (PCClocks/PCBrainAge use `pack$coefficient_matrix`; SystemsAge's pack carries `$organs` +
   `$systems` + `$age`). This maximizes reuse of the tested engine (coverage, transform, reduction,
   impute) -- the alternative of routing all 13 through the orchestrator would re-implement
   `linear_score`'s vendor accounting 11 extra times for no gain. Rejected the earlier
   "hold the organ members back with the composite as one slice" framing: the members are
   independent (the composite recomputes its own raw system vectors), so each of the 13 dispatches
   on its own.

2. **The two composites (`Age_prediction`, `SystemsAge`) are one named branch,
   `score_systemsage()`.** They are not linear in the CpGs, so they get a family orchestrator (like
   GrimAge/FitAge/PhysAge), NOT a recipe walker. The branch hard-codes the pipeline *shape* --
   age-linear front -> quadratic; and for the composite: 11 raw system predictors + poly-scaled age
   -> center/scale -> systems_PCA project -> linear head -- and reads every *constant* from the
   catalog recipe via targeted accessors (`systemsage_step/_age_intercept/_poly/_raw_intercepts/
   _stack_order/_final_intercept`), never hand-copied into R source. Every linear stage reuses the
   shared `linear_predictor()` kernel via a small `sa_linpred()` (present cols + partial-cohort
   cache + vendor offset for absent).

3. **The systems_PCA tensor tree stays in the pack `$tensors`, consumed as-is.** `encode_systemsage`
   already leaves the four small tensors (center[12], scale[12], model[12], rotation[12x12]) in
   `$tensors` keyed by their `weights/` paths; `systemsage_pca()` resolves them by component name
   (from the catalog) and aligns rows/cols to the recipe stack order. Deliberately did **not** bump
   the pack encoding to promote them into named `$pca` fields: the accessor is small, the pack is
   already built/verified, and a re-encode would churn the payload_hash and re-staging for no
   scoring benefit. Column-order care: the pack `$organs`/`$systems` matrices are alphabetical while
   the systems_PCA order is the recipe stack order (Blood, Brain, **Inflammation, Heart**, ...,
   **Metabolic, Lung**, ..., Age_prediction) -- everything is indexed by name, never position.

4. **`Age_prediction` vs the composite's `ap_scaled` use different poly coefs, on purpose.** The
   standalone `Age_prediction` clock folds the author's `transformation_coefs` affine + /12 into its
   quadratic; the composite feeds the un-transformed age (just /12) into systems_PCA, so its
   `ap_scaled` poly differs. Both come straight from their own catalog recipe step, so the scorer
   never conflates them.

## 2026-07-21 -- external scoring landed for PCClocks + PCBrainAge (external.md steps 1-5); SystemsAge deferred

**Decision.** Wired external.md steps 1-5 for the two plain-linear external groups (PCClocks,
PCBrainAge): accessors read external coef/impute from the loaded pack, `load_mc_assets` is the
single loader, `calc_clocks()` resolves packs upfront and scores them on the shared linear engine.
SystemsAge (organ members and its component-matrices composite) stays `"unsupported"` -- its
orchestrator is its own later slice.

1. **`external_pack()` -> `load_mc_assets(groups, assets = NULL, ask = TRUE)`** (renamed, verb-named)
   in [`R/methylCIPHER_data.R`](../R/methylCIPHER_data.R). Returns a **named list of packs keyed by
   `group_id`** (even for one group), not a single pack. It is the sole consent/download/read site.
   The single-pack read/validate/drift-warn body is now the internal helper `mc_read_pack()`;
   `mc_download()` split into pure `mc_fetch()` (stage -> qs2-validate -> atomic rename, no prompt)
   plus `mc_consent()` (one **batched** prompt for the union of missing packs, refusing
   non-interactively). `mc_data_download()` reuses both.

2. **Closed-vs-open `assets`.** `assets = NULL` -> default cache dir, missing packs consent-downloaded
   (open set). `assets` **explicitly provided** -> **closed set, never downloads**, a coverage gap is
   fatal. `assets` accepts a cache-dir path **or** loaded object(s) -- a bare pack, a list of packs,
   or a path all canonicalize (via `mc_canonicalize_assets()`) to the named-list registry; objects
   key by their own `$group_id`. This reverses the old semantics where `assets` was only a cache-dir
   path that downloads landed into. Data-layer tests were rewritten to drive open-mode downloads via
   `options(methylCIPHER.cache_dir=)` and to cover the closed path/object modes.

3. **Accessors take the registry on the external path.** `clock_coefs(id, packs = NULL)` and
   `clock_impute_ref(id, packs = NULL)` gained an external branch: for an `external_group` clock they
   pull the coef column from `packs[[group_id]]$coefficient_matrix[, id]` (named by `$cpgs`) and the
   vendor ref from `$impute`, via the new `clock_pack()` helper (errors if the group's pack is not in
   the registry). The bundled path is unchanged; `packs` defaults to `NULL`, so every existing
   bundled caller (`linear_score`, `score_physage`, tests) is untouched. `linear_score()` gained
   `packs` and forwards it.

4. **Upfront resolution is gated on `score_type() != "unsupported"`, not merely "is external".**
   `calc_clocks()` collects the external groups whose clocks route to a pack-consuming scorer and
   calls `load_mc_assets()` once before the pure loop. `assets`/`ask` were added to the
   `calc_clocks()` signature.

5. **Dispatch flip (step 3), narrowly.** `score_type()` now routes external + `cpg_coefficient` /
   {`linear`,`linear_transformed`} -> `"linear"`, **except** the SystemsAge group, which stays
   `"unsupported"`. So PCClocks + PCBrainAge score (step 4 for them is pure fall-through to
   `linear_score`); SystemsAge's cpg_coefficient organ members are held back with the composite so
   the group is enabled as one orchestrated slice, and the gating in point 4 means SystemsAge never
   downloads-then-errors.

6. **Tests (step 5).** `test-score-external.R` drives `calc_clocks()` end-to-end on PCBrainAge with
   an **in-memory pack** (closed set -> no disk, no network): full-coverage golden, vendor-fill
   golden (absent CpGs from `$impute`), and the closed-set "pack absent" error. Real cohort parity
   is now **pack-gated** too: `skip_if_no_pack()` skips an external clock's parity unless its pack is
   cached, and open-mode `calc_clocks()` then reads the cached pack without a download prompt.

**Why hold SystemsAge back.** Its composite is `component_matrices`/`wrapper` (organs -> systems ->
age -> composite over the pack's `$organs`/`$systems`/`$age` matrices), a genuinely new branch;
enabling only its linear organ members would half-support the group and download-then-error on the
headline composite. It ships as one slice.

---

## 2026-07-21 -- external asset layer simplified: qs2 checksum is the integrity guard, payload_hash drift WARNs, no memoise

**Decision.** Rewrote [`R/methylCIPHER_data.R`](../R/methylCIPHER_data.R) from ~479 to ~230 lines.
The module now reads as one flow: resolve cache dir -> find expected file -> consent-gated download
if absent -> `qs2::qs_read(validate_checksum = TRUE)` -> WARN-only content-hash check -> return
pack. `external_pack(group_id, assets = NULL, ask = TRUE)` is the single runtime entry; see
detail-plan 9.4.

**What changed and why (reverses two earlier stances).**

1. **Transfer integrity is qs2's own `validate_checksum`, not a parallel `file_sha256` recompute.**
   This reverses the 2026-07-20 "kept, deliberately ... runtime download-integrity `file_sha256`"
   line. qs2 already stores and verifies a checksum on read; recomputing a separate sha256 over the
   file was redundant. `file_sha256` stays in `mc_provenance` as a maintainer-side record but is no
   longer read by R (candidate for later removal). The `.part`-stage -> qs2-validate -> atomic
   rename gives the same "a corrupt transfer never counts as cached" guarantee the old sha256 +
   staging dance did, in three lines.

2. **`payload_hash` mismatch WARNs, never errors.** The old `external_pack()` hard-errored on four
   structural identity fields (`group_id`/`encoding`/`encoding_version`/`n_cpgs`). Replaced by one
   `mc_payload_hash(pack)` vs provenance `payload_hash` compare that only `warning()`s. The filename
   already embeds `payload_hash`, so a name match is strong; the recompute is cheap belt-and-braces
   and a drift is a "your scores may not match this version" signal, not a reason to refuse to
   score. (`mc_payload_hash` mirrors sync's `payload_hash_of`; the saved object *is* the stripped
   stable payload, so the hashes line up with no re-derivation.)

3. **No memoise -- confirmed by benchmark, not assumed.** Cold `qs_read` of the three packs is
   ~0.10s (PCClocks 13.5 MB), ~0.13s (PCBrainAge 8.6 MB), ~0.20s (SystemsAge 30.7 MB), and
   `validate_checksum` is free within noise. So there is no session-cache/global-env tier and no
   `clear_mc_cache()`-clears-a-memo semantics; `calc_clocks()` loads each needed group's pack
   once per call and passes it down. This also retires the stale detail-plan 9.4 "session cache ->
   bundled -> R_user_dir -> download" order (external packs are never "bundled").

4. **`ask = TRUE` default, not `ask = interactive()`.** `interactive()` as the default silently
   auto-downloads in scripts/CI (evaluates to `FALSE` = "consent given"). `ask = TRUE` prompts when
   interactive, **refuses** non-interactively, and `ask = FALSE` is the explicit-consent signal --
   matching the "no download without consent" invariant.

**Cut as over-engineering** (the module treated our own compiled `mc_provenance` as a hostile wire
contract): `mc_validate_row` field-by-field validation, `mc_resolve_groups` multi-group ceremony,
`mc_cache_status` data.frame, `mc_release_repo` slug regex, `mc_chr1`/`mc_num1`, the elaborate
`mc_confirm`, and the separate sha256 verify + `.part-<pid>` machinery.

**Cache-dir model.** `assets` arg (per-call) > session option `methylCIPHER.cache_dir` (via
`mc_set_cache_dir()`, a session-consistent setter) > `METHYLCIPHER_CACHE_DIR` (.Renviron) >
`mc_default_cache_dir()` (`R_user_dir`).

**Kept.** `clear_mc_cache()` as a must-have **stub** -- it reports what is cached and never
auto-unlinks; the interactive delete flow is deferred by design. A thin `mc_data_download()`
pre-fetch verb for offline/CI staging.

**Sites.** `R/methylCIPHER_data.R` (rewrite); `tests/testthat/{test-methylCIPHER_data.R,
helper-external-data.R}`; `dev/detail-plan.md` 1.1, 9.1, 9.4. Full suite after: 0 failures (the 13
parity warnings are the pre-existing per-clock-subset orientation notice, unrelated).

## 2026-07-21 -- no roxygen yet; plain `#` comments are the only in-source docs pre-alpha

**Decision.** The package carries **no roxygen** during the rewrite. Do not author roxygen
blocks and do not run `devtools::document()`. In-source documentation is short `#` comments
only (the existing "Comments" rule: 1-2 sentences on *what*, not *why*).

**Why.** The API surface is still moving; regenerating `NAMESPACE` / `man/*.Rd` on every change
is churn with no reader yet, and half-written roxygen would rot against the code. `man/*.Rd` and
the live `NAMESPACE` remain the hand-managed pre-rewrite leftovers noted below until roxygen is
switched on.

**Trigger.** Turning roxygen on is a **human-decided override** tied to the alpha release --
there is **no** automatic condition (no version tag, no milestone gate). Claude must not enable
it on its own initiative.

**Supersedes.** The earlier CLAUDE.md guidance to "run `devtools::document()` after any
export/doc change" and to treat `NAMESPACE` / `man/*.Rd` as roxygen-generated. Those lines were
rewritten to the "No roxygen yet" invariant.

## 2026-07-21 -- maintainer handoff to Hung Pham; license is BSD-3 (not the anticipated GPL-2)

**Decision.** Two DESCRIPTION-level facts recorded now that they are set in package metadata.

1. **Authorship / maintainer.** Hung Pham (ORCID 0000-0002-8271-9355) is now `aut` + `cre`
   (maintainer). Kyra L. Thrush stays `aut` (prior author). Copyright holder stays the **Albert
   Higgins-Chen Lab, Yale University** (`cph`) -- Hung takes authorship / maintenance, not copyright.
   DESCRIPTION moves from the malformed plain-text `Author:` / `Maintainer:` pair (which held a
   literal `person(...)` string) to a proper `Authors@R` block, clearing that `checking DESCRIPTION
   meta-information` NOTE.

2. **License is `BSD_3_clause`, superseding the GPL-2 anticipation.** The 2026-07-17 third-party
   weights entry assumed a future GPL-2 finalization; the package instead ships `BSD_3_clause + file
   LICENSE` (LICENSE stub + LICENSE.md, copyright the Lab). **Caveat unchanged:** BSD-3 is more
   permissive than GPL-2 but does **not** dissolve the redistribution question for bundled
   third-party clock weights carrying NC / ND / research-only / unspecified terms. That
   pre-submission gate (chase grants, lean on the facts stance, or move to the external tier) still
   applies whatever license the package itself ships under; only the CRAN GPL-2-sublicensing framing
   from 2026-07-17 is moot.

**Sites.** `DESCRIPTION` (`Authors@R`, `License`); `LICENSE` / `LICENSE.md` (copyright unchanged).

---

## 2026-07-21 -- plan/code reconciliation: sync §9 lockfile story, dead `impute_DNAm.R` path

**Decision.** Documentation-only reconciliation of `dev/migration-plan.md` and `dev/detail-plan.md`
with what `data-raw/sync.R` actually does after the 2026-07-20 prune. No behavior change.

1. **`manifest_key` is gone from the plans.** The 2026-07-20 prune removed the catalog build-skip
   cache (`manifest_key` / `build_key` / `data-raw/lockfile.json`), but both plans still described it
   as current truth in six places (migration §Source-of-truth + §Packaging; detail §7, §9.1, §9.2,
   §14 diagram). Rewritten to the current reality: the full `sysdata.rda` rebuild is uncached (~2s),
   and the only remaining sync-side identity is the external packs' content-address `payload_hash`
   (filename / release tag / reupload skip) plus runtime `file_sha256` in `mc_provenance`.
2. **The lockfile survived, narrower.** The 2026-07-20 entry below overstated the removal: a
   gitignored `data-raw/assets/lockfile.rds` (+ `read/write_lockfile`, `lockfile_hit`, `sync(force=)`)
   remains, but its only job now is to skip **external qs2 pack rebuild** when the meta
   `source_git_sha` and the staged packs are unchanged -- not catalog-change detection. detail §9.1 /
   §9.2 now describe this instead of the old `manifest.json -> manifest_key -> no-op` flow.
3. **Dead path fixed.** The plans (and the 2026-07-18 entry below) referenced `R/impute_DNA.R`; the
   file is `R/impute_DNAm.R`. Corrected in detail §2.3a; the 2026-07-18 entry is left as-was
   (append-only).

**Why documentation-only.** The plans state current truth; these were plan<->code drifts left by the
2026-07-20 prune, not new design. DECISIONS stays append-only, so this entry records the drift rather
than silently rewriting history.

**Sites.** `dev/migration-plan.md` (Source of truth, Packaging); `dev/detail-plan.md` §7, §9.1, §9.2,
§9.5, §14, §2.3a.

---

## 2026-07-20 -- `data-raw/sync.R` pruned to two requirements: consume upstream, don't re-upload unchanged assets

**Decision.** Cut `data-raw/sync.R` from 2172 to ~1350 lines (~40%). The script's only jobs are (1)
pull the scoring contract + tensors from methylCIPHER-meta into `R/sysdata.rda` + `inst/` + the three
external qs2 packs, and (2) not re-upload an external pack whose stable weights didn't change.
Everything not serving those two was removed.

**Removed.** (a) The **build-skip cache** entirely -- `read/write_lockfile`, `manifest_key`,
`build_key`, `missing_outputs`, `data-raw/lockfile.json`, the `prior_assets` reuse path, and the
`force`/`dry_run` modes. The full build takes ~2s, so caching it bought nothing but cache-invalidation
edge cases (which `missing_outputs` existed to patch). (b) **Trust-upstream input validation** --
manifest field/dup checks, `papers.csv` uniqueness/column checks, `assert_repo_rel` path-escape
guards, `verification_status` threading, sha-format and post-checkout existence assertions. methylCIPHER-meta
owns provenance integrity; re-policing it here is redundant. (c) The `manifest_generated_at_sha`
parent-commit lag-stamp. (d) ~250 lines of rationale-essay comments -- the "why" lives here in
DECISIONS.md and `dev/*.md`.

**Kept, deliberately.** The scoring-correctness guards, because a silent failure there produces a
confident wrong number (the class of bug the whole design fights): `align_double`/`cbind_aligned`
alignment checks, the missing-tensor stop in `build_group_bundles`, `resolve_cpgs` probe-set-mismatch,
`materialize_probe_set`'s empty/duplicate-CpG stops, and `read_tensor_csv`'s duplicate-key stop. Also
kept the runtime **download-integrity `file_sha256`** in the shipped `mc_provenance` registry (protects
the consumer, not the maintainer) and the git subprocess error surfacing.

**Why no re-upload without the lockfile.** Idempotency was never the lockfile's job. The content-address
(`payload_hash` -> filename -> release tag) plus `gh_upload.py`'s remote skip (`get_release(tag)` +
"asset name already present") together guarantee it: unchanged weights -> same hash -> tag+asset already
on the release -> skipped. Upstream moving without touching the three external groups leaves their
hashes untouched (pin fields are stripped by `stable_external_payload`). Verified: a full `sync()`
rebuilds all outputs in ~2s and the three packs keep stable hashes across runs.

---

## 2026-07-20 -- external asset content-address hashes over a version-pinned serialization, not `rlang::hash`

**Decision.** `payload_hash_of()` in `data-raw/sync.R` now hashes `serialize(payload, connection =
NULL, version = 2L, xdr = TRUE)` with `digest::digest(algo = "sha256")`, replacing `rlang::hash(payload)`.
`EXTERNAL_ENCODING_VERSION` bumped 2 -> 3; `rlang` dropped from the sync `library()` preflight (it
was the only use). This is the content-address for the three external qs2 release assets (SystemsAge,
PCClocks, PCBrainAge) -- it sets the filename, the GitHub release tag, and the reuse/skip key.

**Why.** `rlang::hash()` hashes R's *serialization* of the object, so it inherits ALTREP compaction
and serialization-version drift: an R/rlang upgrade could flip the hash with **zero content change**
and force a spurious rehash + reupload of a blob that is meant to be frozen. `version = 2L` predates
ALTREP (the object is written fully expanded) and `xdr = TRUE` fixes endianness, so the byte stream
depends only on content. `rlang::hash` exposes no knob to pin its serialization, so it could not be
fixed in place. The new form also matches the `digest::digest(..., serialize = FALSE)` idiom already
used by `manifest_key()` and `build_key()` -- `rlang::hash` was the lone outlier.

**Why the one-time reupload is acceptable.** Switching the hash function is itself a genuine
content-address change: the three current blobs rehash once on the next `sync(upload = TRUE)`, then
never move again unless their group's own contents change. The package is private and pre-release, so
there are no pinned installs in the wild to strand, and old release tags are never deleted anyway.

**Scope note.** Only SystemsAge/PCClocks/PCBrainAge are content-addressed qs2 assets. Clocks added to
ship groups (GrimAge, Dunedin, ...) live in `mc_bundles` inside `R/sysdata.rda` and never touch these
hashes -- so in normal development the three external hashes are effectively constants, which is
exactly why removing the R-version drift vector matters.

---

## 2026-07-20 -- dependency clocks are returned as columns (reverses "compute-only intermediates"); `sim_DNAm` mirrors the compute plan; a prepare-time coverage floor

**Decision.** Three coupled changes, all traceable to one bug: `calc_clocks(sim_DNAm("DNAmFitAge"),
"DNAmFitAge")` returned a confident, meaningless number.

1. **Auto-added dependency clocks are RETURNED as score columns.** This **reverses** the
   2026-07-17 rule that they are compute-only intermediates. `calc_clocks()` assembles
   `results[output_ids]` where `output_ids = c(clock_ids, setdiff(clock_sequence, clock_ids))` —
   requested first in request order, deps after in compute order — and
   `$provenance$requested` / `$dependencies` partition `$provenance$clocks`.
2. **`sim_DNAm()` resolves through `resolve_clocks_sequence(resolve_clocks(clocks))`**, not
   `resolve_clocks()` alone.
3. **A coverage floor at prepare time.** `warn_low_coverage(cpg_list, min_coverage = 0.8)` in
   `R/resolve_inputs.R`, called from `calc_clocks()` right after `resolve_cpgs()`.

**The bug that forced all three.** `sim_DNAm` resolved namespace only, so `"DNAmFitAge"` produced
the group's 627 CpGs while the compute plan needs 1643 — the 1016 missing belong to `GrimAgeV1` and
its 8 surrogates, pulled in as cross-group deps. Those surrogates are `imputation_policy = omit`, so
at 0 present CpGs they collapse to intercept-only, `GrimAgeV1` is built on them, and `DNAmFitAge`'s
KDM `DNAmGrimAge` term becomes a constant. It did not error. It returned 37–67 where the correct
panel gives 90–127. Three independent layers were each silent, hence three fixes.

**Why return deps rather than propagate their coverage into the composite's row.** The considered
alternative was to keep deps hidden and fold their coverage into `DNAmFitAge`'s tier-1 row (or add a
`dep_coverage` field). Returning the columns is strictly better: the deps *were* computed, each
already carries its own complete coverage row, and stacking them costs nothing. Folding would also
make `score_needed` stop matching the clock's own panel. Consequence: `score_fitage_composite()` and
`score_grimage()` need no change at all — the transparency problem was in assembly, not in them.

**Why not keep deps hidden for a stable column set.** That was the original argument, and it loses
to the failure above: `calc_clocks(x, "DNAmFitAge")` returning 16 columns is surprising once; a
composite silently resting on a collapsed input is dangerous every time. Callers who want a fixed
set subset downstream, which is cheap and explicit. Note the column set was never truly stable
anyway — it already varied with group-vs-clock token expansion.

**Why the floor warns instead of erroring.** A partial panel is often legitimate: `vendor_mean`
clocks fill, `omit` clocks degrade smoothly, and a user may knowingly score a 27k array. Silence is
what is indefensible. 0.8 is a judgment call, not a derived threshold — hence `min_coverage` is an
argument, with `0` to silence.

**Why the floor runs over the plan, not the request.** That is the whole point: the collapsed clock
was `GrimAgeV1`, which nobody requested. Over the request alone the motivating bug stays invisible.

**Why it is not inside `resolve_cpgs()`.** `resolve_cpgs()` is pure set math over the catalog (no
numerics, no `pheno`, no side effects) — 2026-07-18 established that purity, and scope creep into it
was already flagged once as a tell. The floor is a separate consumer of the skeleton it returns.

**Why the floor is sample-invariant only.** It reads `score_present / score_needed`, so it catches
panel/array mismatch, never per-sample NA concentration. The per-sample analogue (a threshold read
off the tier-2 `n × k` matrix, §4.2) is deliberately **left unbuilt** — record-and-report stays the
default there, and no per-sample gate ships until there is a concrete case for one.

**Sites.** `R/sim_DNAm.R` (`clock_cpgs` call site), `R/calc_clocks.R` (`output_ids`,
`min_coverage` arg, `construct_methylCIPHER(requested_ids)` + provenance fields),
`R/resolve_inputs.R` (`warn_low_coverage`). Plan text: detail-plan §1.3, §2.1, §4.1, §7, §10.1.

---

## 2026-07-18 -- partial-NA is a shared cohort cache; all-NA columns reclassify absent; empty rows throw; coverage gains a per-sample tier

**Decision.** Four coupled rules for the missingness front end, settling how partial-NA (cohort)
imputation is organized and how degenerate inputs are handled. Numeric machinery in
`R/impute_DNA.R` (`scan_missing_cpgs`, `build_partial_cache`); pure set math in `R/resolve_inputs.R`
(`resolve_cpgs`). Order in `calc_clocks()`: `scan_missing_cpgs → resolve_cpgs → build_partial_cache`.

1. **Partial-NA fill is a shared, precomputed cache — not per-scorer.** Component families reuse the
   same CpGs, so per-clock cohort-mean imputation would recompute the identical column mean once per
   clock. Build it **once** over `intersect(present_needed_union, partial_na_cols)` via
   `slideimp::mean_imp_col(DNAm[, cache_cpgs])`; every scorer reads the same filled column. This is
   also a **consistency** guarantee (shared CpG ⇒ identical imputed value across clocks), not only
   speed. Subset-columns-first is load-bearing: `mean_imp_col` returns an input-width matrix, so
   running it on the full panel allocates a second `n × p` copy (the §3 spike); narrowing bounds the
   cache to the partial-NA needed probes (empty when clean).
2. **All-NA present column → reclassify ABSENT.** Zero observed values ⇒ no cohort mean exists ⇒
   informationally identical to a probe off the panel. Removed from `usable_cols` via `setdiff` (a
   name-set op, **not** a `DNAm[, keep]` copy) so it routes to the vendor/drop path and is counted in
   `n_cpg_score_miss`. `mean_imp_col`'s own "0 observed → unchanged" branch therefore never fires in
   our path.
3. **Empty sample (all-NA row) → HARD ERROR (no permissive dial).** The reason it is an error, not a
   warning: cohort-mean imputation would otherwise fill its every cached cell with the column mean and
   emit a plausible-looking **fabricated** score instead of an honest `NA`. (An all-NA row is
   `isnan`-skipped in every column mean, so it does not bias others — the damage is its own fake
   score, but that suffices to refuse it.) Unlike rowname-less DNAm (§5.1 `allow_positional_ids`),
   there is no opt-in; empty rows are always fatal.
4. **Coverage gains a tier-2 per-sample matrix.** Tier 1 (`summary()`) stays the per-clock aggregate;
   tier 2 is an `n × k` per-clock × per-sample missingness matrix (`sample_coverage(x)`) for a QC
   heatmap. Motivation: a sample fine *globally* can be ~100% missing *for one clock* if its gaps
   concentrate in that clock's set — and that is **recorded, never thrown**. Each scorer emits a
   length-`n` row-wise NA count over its own present scoring CpGs, computed on the **raw** subset
   before the cache fills it (post-fill reads ~0). It is the row-wise decomposition of tier 1's
   `score_imputed_partial`, so no new mechanism.

**Gate cascade (both cache and tier 2).** Global `anyNA(DNAm)` is the monotone loose upper bound —
clean betas skip everything. Per clock, `intersect(score_present, partial_na_cols)` is a **free** set
op reusing the single prepare-side `mat_miss(col = TRUE)` result; a local `anyNA(subset)` gate is
rejected — on a clean subset it full-scans (nothing to short-circuit on), costing the same as
`mat_miss` while telling you less. `mat_miss(col = FALSE)` then runs only over a clock's NA-bearing
columns.

**Boundary (why `resolve_cpgs` stays pure).** `resolve_cpgs` does **no** numerics and takes **no**
`pheno` (an earlier draft signature `resolve_cpgs(DNAm, pheno, clock_ids)` was the scope-creep tell —
`pheno` is only needed by sex-keyed *vendor* fill, which is per-clock and downstream). It maps
`usable_cols` + the compute plan to the per-clock present/absent skeleton + `present_needed_union`.
Partial (cohort, shared cache) and absent (vendor, per-clock, needs `Female`) stay two separate
inputs to one gather; they are never merged upstream (the §2.3 never-cross invariant).

**Rejected.** (a) Per-clock partial-NA fill (recomputes shared column means; loses the cross-clock
consistency guarantee). (b) A union submatrix / full-panel `mean_imp_col` (the §3 memory spike; the
bounded cache is a narrow, evidenced exception to §3's blanket "no union submatrix" and §12's "no
copy-opt before evidence" — the evidence is shared component CpGs, and it is empty in the clean
case). (c) Warning (not error) on empty rows (lets imputation fabricate a cohort-average score). (d)
Baking absent counts into the tier-2 vector (absent is sample-invariant — keep it as tier 1's scalar,
combine at render).

**Sites.** `R/impute_DNA.R` (`scan_missing_cpgs`, `build_partial_cache`), `R/resolve_inputs.R`
(`resolve_cpgs`), `R/calc_clocks.R` (wiring), `R/score.R` (scorer emits tier-2 vector; consumes
cache), `dev/detail-plan.md` §2.3a, §4.1/§4.2. Uses `slideimp::{mat_miss, mean_imp_col}` (already in
Imports).

---

## 2026-07-17 -- clock resolution: maximal group expansion + a `resolve_clocks_sequence` dep topo-plan

**Decision.** Two settled rules for turning user tokens into a compute worklist, both in
`R/resolve_inputs.R`:

1. **`resolve_clocks()` — maximal expansion on a name clash.** A token resolves by fixed
   precedence, MOST clock ids first: `"all"` -> every clock; a `group_id` -> all members; a bare
   `clock_id` -> itself. When a token is **both** a `group_id` and a `clock_id` (e.g. `"DNAmFitAge"`
   names both the 7-member group and one composite member), the **group wins** -- a group request
   returns as many member ids as possible. One rule covers every shape: a group_id that is no
   clock_id (`GrimAge`, `SystemsAge`, `PCClocks`), a group_id that is also a member id
   (`DNAmFitAge`), and a singleton group (`group_id == its sole clock_id`, e.g. `MiAge`); the last
   two coincide at one member. This was previously an *accident* of `c()` ordering + `[`
   first-match; now it is explicit (`resolve_one()` checks group before clock) and documented.
2. **`resolve_clocks_sequence()` — the dependency plan step.** This settles the "orchestrator
   nesting vs. a resolve/plan step that topo-sorts on `depends_on_clocks`" question that
   detail-plan §2.1 left open: it is a **topo-sort plan step**. It takes `resolve_clocks()` output
   and returns the **transitive closure** over `depends_on_clocks` (read via the new
   `clock_depends_on()` accessor), topologically ordered so every clock comes AFTER its deps, with a
   cycle guard. Only 3 clocks declare any dep (`GrimAgeV1`, `GrimAgeV2`, `DNAmFitAge`), so like the
   covariate union, deps are pulled in **on demand** -- `SystemsAge` and its organ/system members
   declare none and pass through unchanged.

**Coupled fix.** Auto-added deps (in the plan, not in the requested set) are compute-only
intermediates; the caller returns score columns for the requested set only. Because a dep can need
a covariate the request does not, `calc_clocks()`'s covariate union (and the `check_DNAm_extra`
worklist) now run over the **plan**, not the raw `clock_ids`: the single `DNAmFitAge` clock declares
only `Female`, but its `GrimAgeV1` / `DNAmVO2max` deps also need `Age`. (In current data this is a
robustness fix -- the `"DNAmFitAge"` *token* resolves to the group, which already includes an
Age-needing member -- but it is reachable for the single clock id and if the catalog shifts.)

**Why.** `resolve_clocks` stays a pure token->id map; dependency logic lives in one named plan step
rather than smeared into each pack orchestrator, mirroring how the covariate union is computed once
over the worklist. Maximal expansion is the least-surprising reading of a group name and is what the
code already did.

**Rejected.** (a) Clock-beats-group (minimal expansion) on a clash -- makes `DNAmFitAge` return one
column when the user named the group. (b) Pure reorder without auto-adding missing deps -- leaves
`resolve_clocks_sequence` unable to actually satisfy a `DNAmFitAge`-only request. (c) Handling deps
inside the GrimAge/FitAge orchestrators only -- keeps the covariate-union coupling invisible and
duplicates topo logic per pack.

**Sites.** `R/resolve_inputs.R` (`resolve_clocks`, `resolve_clocks_sequence`), `R/accessors.R`
(`clock_depends_on`), `R/calc_clocks.R` (plan step + union over plan), `dev/detail-plan.md` §2.1.

---

## 2026-07-17 -- rowname-less DNAm: positional ids, refuse at cbind (not a front-door error)

**Decision.** `rownames(DNAm)` is no longer a hard requirement at the front door. `calc_clocks()`
handles identity **inline** before `check_DNAm`: when rownames are absent it stamps positional ids
`sample1..sampleN` (gated on `allow_positional_ids`). Such a call is flagged
`$provenance$positional_ids = TRUE`, and `cbind.methylCIPHER` gains a **gate 0** that throws on any
positional record. `check_DNAm` keeps `null.ok = FALSE` (an unconditional invariant, since ids exist
by the time it runs). `calc_clocks(..., allow_positional_ids = FALSE)` restores the strict
hard-error behavior. (Stamping is silent -- no warning; the `cbind` gate 0 is the guard.)

**Why.** The prior contract hard-errored on rowname-less DNAm to prevent a `cbind` footgun (two
matrices with `seq_len` ids binding as if the same samples). But that put maximal defensiveness at
the front door, which some workflows don't want -- they legitimately score by row order and never
`cbind`. The footgun only bites at the *bind* step, so guard it *there*: manufacture ids at the
boundary (keeping the "unique id everywhere" invariant every downstream consumer, incl. `cbind`,
relies on), warn, and quarantine positional records from `cbind` rather than corrupting it. Real
rownames remain preferred and are what enable binding + `Reduce(cbind, ...)`.

**pheno alignment follows the mode** (`resolve_pheno()`). When DNAm has real rownames, pheno is
**id-joined**: `rownames(DNAm) ⊆ pheno[[pheno_id]]` (superset ok, subset + reorder). When DNAm is
positional, there is no id to match on, so pheno is **row-order joined**: `pheno_id` ignored,
`nrow(pheno)` must equal `nrow(DNAm)` **exactly** (taller = ambiguous, shorter = short), row *i* is
sample *i*. This is the only way covariate-needing clocks can pair Age/Female to rowname-less
samples. Both modes keep the record scores-only; the aligned pheno feeds scorers, it is not
appended.

**Rejected.** (a) Flipping `check_DNAm`'s `null.ok = TRUE` and letting a positional `sample_id`
flow unlabeled -- reopens the silent-bind footgun and leaks positional semantics into
provenance/cbind. (b) Row-order join with a **superset** pheno (`nrow(pheno) >= nrow(DNAm)`) -- in
positional mode a taller pheno gives no way to know which rows are the samples, so require exact
equality. (c) Auto-appending pheno onto the scores (`cbind(scores, pheno)`) -- the record stays
scores-only regardless of mode; joining is `augment()`'s / `as.data.frame()`'s job.

**Sites.** `R/calc_clocks.R` (`allow_positional_ids` arg + inline positional-id stamping),
`R/resolve_inputs.R` (`resolve_pheno`, `check_DNAm` / `check_pheno`), `R/generics.R` (identity
contract + cbind gate 0), `dev/detail-plan.md` §1.3/§5.1/§7/§7.1, `dev/migration-plan.md`
(Public API).

---

## 2026-07-17 -- drop `covers` / `shared` from the shipped catalog (runtime-catalog trim)

**Decision.** `covers` and `shared` are no longer shipped on `mc_catalog` entries. They stay in
`FIELD_REGISTRY` (so the build-time resolver can still read them), then a new
`CATALOG_BUILD_ONLY_FIELDS` strip in `build_sysdata()` removes them from every clock entry
immediately after `resolve_group_scoring_probe_sets()`.

**Why.** Both duplicate `mc_groups`, and are dead weight at runtime: `covers` (composite
membership) mirrors `mc_groups$members` and is only ever consumed by the build-time CpG resolver
(`resolve_scoring_cpgs()` tier 5); once `probe_sets[role=scoring]` is materialized nothing reads
it again. `shared` (clock-level shared-tensor name list) mirrors `mc_groups$shared_tensors` and
is redundant with `imputation$ref` + the group bundle -- it was write-only (pruned in, never read
by sync's logic). Removing them sharpens the invariant "an `mc_catalog` entry is one clock's
runtime scoring contract" and drops two drift-prone copies of group facts. (~662 -> ~647 KB.)

**Ordering (load-bearing).** The strip runs *after* resolution, not by dropping the fields from
`FIELD_REGISTRY`: `resolve_scoring_cpgs()` reads `entry$covers` from the *pruned* catalog entry,
so deleting it pre-resolution would silently empty composite scoring sets (GrimAgeV1's 1030-CpG
union resolves through its `covers` list). Verified post-rebuild: GrimAgeV1 still resolves 1030
scoring CpGs with `covers`/`shared` absent from all 113 entries.

**Kept.** `components` and `recipe` stay -- they are the pack scorers' runtime instruction set,
not group duplication. Membership at runtime comes from `mc_groups$members`, never a resurrected
`covers`. Added `clock_covariates_coefs()` so the just-confirmed uniform `$covariates` coef map is
read through the accessor layer.

**Sites.** `data-raw/sync.R` (`CATALOG_BUILD_ONLY_FIELDS`, `build_sysdata`); `R/accessors.R`
(`clock_covariates_coefs`, `clock_coefs` now routes on `weights_format`).

---

## 2026-07-17 -- GrimAge V1/V2: do not symmetrize upstream; components emit V2 precision

**Decision.** Leave `weights/GrimAge/v1` and `v2` as-is upstream. Do **not** restructure
GrimAgeV1 to carry an inline component set mirroring V2. The GrimAge pack orchestrator
(`score_grimage`) absorbs the version asymmetry.

**Why.** The asymmetry encodes a real published difference, not a modeling accident: V1's 8
surrogates *are* the exposed standalone component clocks (DNAmADM, DNAmB2M, ...), so V1 reuses
them via `depends_on_clocks` with no duplication; V2's 8 surrogates differ by ~1e-4 to 1e-9,
were never published as standalone clocks, and so live as private `v2/_internal` components.
Forcing symmetry would either duplicate the published surrogate coefficients into GrimAgeV1's
meta (a one-writer violation) or mint ~1e-4-apart near-duplicate clock_ids that pollute
`list_clocks()`. Downstream gains nothing: `score_grimage` is a hand-written orchestrator
(detail-plan sec 2.5), explicitly not a generic recipe path, so it reads V1's `depends_on_clocks`
and V2's `components` and does the right thing per version.

**Component columns come from the V2 path** (detail-plan sec 2.5): when a user requests the
component id `DNAmADM`, the pack emits the V2-precision value from `v2/_internal/DNAmADM`; the
standalone V1-precision `DNAmADM` clock_id stays available by explicit name. The `_internal`
tensors do not become addressable clock_ids.

**Rejected.** (a) Inline components on V1 for uniformity -- duplicates published coefs. (b)
Expose V2 `_internal` surrogates as clock_ids -- near-duplicate namespace pollution. If the
`component_matrices` format doing double duty (surrogate-by-ref vs surrogate-by-inline-component)
ever bothers a schema consumer, that is a `weights_extraction.md` doc clarification upstream, not
a data restructure.

**Sites.** Upstream `weights/GrimAge/` (unchanged); future `score_grimage` in `R/`.

---

## 2026-07-17 -- sync identity keys no longer bleed into shipped `sysdata` (F1/F2)

**Decision.** Enforce sec 9.1 literally -- the sync-internal identity keys live **only** in
`data-raw/lockfile.json`, never in `R/sysdata.rda`:

- **F1.** `bundle_hash` and `out_sha256` are no longer written onto `mc_catalog` entries
  ([sync.R](../data-raw/sync.R) `build_catalog`). They were write-only (no reader in `R/`, in the
  accessors, or in sync itself -- `manifest_key()` reads `man$clocks` directly) and were leaking
  into the external `payload_hash` via `bundle$catalog`, so a manifest-only hash change could
  re-mint a release with identical scoring content. The upstream-freshness reasoning is kept as a
  comment; the values are not recomputed or stored.
- **F2.** `manifest_key` is no longer put in `mc_provenance`, and `catalog$manifest_key` is no
  longer set. It reaches the lockfile via `mkey -> write_lockfile()` as before. `source_git_sha`
  and `manifest_generated_at_sha` **stay** in `mc_provenance` as honest build provenance (only
  written on rebuild, so they accurately name the bytes -- not covered by the lockfile-only
  clause). The `external_assets` registry stays too (runtime download integrity, sec 9.4).

Also completes the `verification_status` code follow-up flagged in the reconciliation entry
below: it is dropped from both `mc_index` and `mc_catalog` (validated maintainer-side only).

**Why.** These fields answered exactly one question -- "has the catalog changed since last
sync?" -- which is the lockfile's job. Shipping them in `sysdata` contradicted sec 9.1, bloated
the catalog, and coupled the content-addressed release hash to manifest bookkeeping.

**Rejected.** Rewording sec 9.1 to bless identity keys in `sysdata`. Nothing consumes them there,
so removal is strictly cleaner than relaxing the rule.

**Sites.** `data-raw/sync.R` (`build_catalog`, `build_sysdata` `mc_provenance`, `sync`);
`dev/detail-plan.md` sec 9.1. A future `names(mc_index)` structural test must match the trimmed
schema (no `verification_status`).

---

## 2026-07-17 -- third-party weight licenses: audit + CRAN-vs-GitHub contingency

**Context.** A GPL-2 CRAN package that **bundles** third-party clock coefficients may only
redistribute weights whose own license permits GPL-2 sublicensing. CRAN additionally rejects
non-commercial, no-derivatives, research-use-only, and unlicensed ("all rights reserved")
content outright. Audited `weights/**/*.meta.json`'s `license` field (already kept via sync's
`FIELD_REGISTRY`, so this is queryable from the catalog).

**Findings.**

- **Clean (GPL-2-redistributable / CRAN-OK):** MIT (CellPopAge, Knight, PCBrainAge); GPL-2(+)
  (Bohlin96/251, DunedinPoAm38); GPL-3 (StocClocks P/H/Z); GPL-3 / CC-BY-4.0 (EpiTOC2,
  HypoClock); CC BY (Horvath1, RepliTali(Norm), RetroAge 450K/EPICv2); `open-redistributable`
  (HRSInChPhenoAge, all PCClocks, all SystemsAge, all SenescenceAge).
- **Blocked / unclear (NOT clearly redistributable; CRAN blockers):**
  - `non-commercial`: IntrinClock 370/380, the entire **PhysAge** family (DNAmCRP, DNAmHbA1c,
    DNAmDHEAS, DNAmHDL, DNAmPeakflow, DNAmPhysAge(_years), DNAmWHR, DNAmPulsePr,
    DNAmCystatinC_PhysAge).
  - `CC BY-NC-ND 4.0`: Mayne (NC **and** ND).
  - `research-use-only`: DunedinPACE.
  - `unspecified`: entire **GrimAge** family (all DNAm* components + GrimAgeV1/V2), MiAge.
  - `public-github-unspecified` (= no grant = default all-rights-reserved): DNAmFitAge family,
    CellDRIFT, DNAmClockCortical, Zhang2019, PedBE.
  - `journal-supp` (published as supplementary; license usually unstated): Hannum, PhenoAge,
    Horvath2, Lin, Weidner, McCartney family, Lee* placental, EpiTOC, DNAmTL, DNAmStress,
    VidalBralo, ZhangMortality, DNAmFI_Li, DNAmIC, CausalityAge (Caus/Dam/Adapt).

**Reasoning.** The technical "can't ship blobs needing rev-eng" concern (DunedinPACE's historic
obfuscation) is **orthogonal** to license: even plaintext, freely downloadable coefficients can
carry NC/ND/unspecified terms that forbid GPL-2 redistribution. A counter-argument -- raw
regression coefficients may be uncopyrightable facts (Feist), making the labels unenforceable on
the numbers -- exists but CRAN will not adjudicate it; a clean grant is required regardless.
Moving to a GitHub release **sidesteps CRAN policy** but does not by itself cure the underlying
redistribution grant (it leans on the facts argument + academic precedent, e.g. the current
methylCIPHER already redistributes many of these).

**Decision (contingency, not yet locked).** Treat license as a **pre-submission gate**: before
any CRAN attempt, resolve each non-clean clock -- chase the author for an explicit grant, rely on
the facts stance, or move it to the external / user-fetched tier -- so the CRAN-bundled set is
clean-only. If that leaves too little to be worth CRAN, ship the whole package as a **GitHub
release** instead (accepted fallback). License compatibility, not the stale `DESCRIPTION` string,
is the real CRAN risk.

**Sites.** `data-raw/methylCIPHER-meta/weights/**/*.meta.json` (`license`); `dev/detail-plan.md`
sec 9.3 / sec 10.2 (packaging + CRAN checklist); `DESCRIPTION` (License, when finalized to GPL-2).

---

## 2026-07-17 -- plan/code reconciliation: bibliography input, verification_status is dev-only, sim_DNAm smoke, self-contained `%||%`

**Decision.** Four corrections aligning the plans (and one code hygiene fix) with what
`data-raw/sync.R` and `R/sim_DNAm.R` actually do:

1. **Sync inputs include `bibliography/`.** The "inputs R may read" contract now lists
   `bibliography/{papers.csv,clocks.bib}` alongside `manifest.json` + `weights/**`. Sync joins
   `pmid -> bib_key` from `papers.csv` in memory and vendors `clocks.bib` to `inst/`; `papers.csv`
   is never shipped. (`control/`, `papers/`, `scripts/` remain unread.)
2. **`verification_status` is a maintainer/dev signal, not user surface.** It records why a clock
   is deliberately `skipped` (e.g. Horvath online variants) for the human running sync. It is
   **not** exposed by `list_clocks` and should not be canonicalized into `mc_index` / `mc_catalog`
   -- it stays on the maintainer side (manifest / lockfile). Sync still validates the value at
   build time. **Code follow-up (not yet done):** drop `verification_status` from `build_index`
   ([sync.R:1366](../data-raw/sync.R)) and the `mc_catalog` entry ([sync.R:766](../data-raw/sync.R));
   nothing in `R/` consumes it, and it needs a re-sync to take effect.
3. **`sim_DNAm` smoke tier documented.** The always-run test tier now names the shipped
   `sim_DNAm()` `expect_no_error` machinery check over every bundled clock (random uniforms on the
   catalog's scoring probe sets; no cohort, no golden, no drift). Distinct from the hand-authored
   engine-arithmetic units and from the cohort-gated parity gold.
4. **Self-contained `%||%` in sync.** Restored the small `%||%` helper (was commented out, leaving
   the script silently dependent on base R >= 4.4 shipping its own). Audited all ~55 call sites:
   every LHS is freshly built or from `fromJSON(simplifyVector = FALSE)` (NULL, never scalar NA),
   so base's null-only `%||%` and the restored NA-aware one are behaviourally identical on shipped
   output; the restore only removes the hidden version dependency.

**Why.** These were plan<->code drifts I would rather kill than annotate -- the plans state current
truth only. None reverse a design; they correct descriptions that had fallen behind the
implementation.

**Sites.** `dev/migration-plan.md` (Source of truth, Testing); `dev/detail-plan.md` sec 5, sec 8,
sec 9.5, sec 10.1; `data-raw/sync.R` (`%||%`; `build_index`/`build_catalog` follow-up);
`R/sim_DNAm.R`.

---

## 2026-07-17 -- runtime is `linear_score` + a finite branch set (walker fully dropped)

**Decision.** No recipe interpreter, not even as a reference/fallback. Runtime is one shared
`linear_score()` plus a small, closed set of branches dispatched on the catalog pair
`(weights_format, computation_type)`: linear, linear + pre-transform (Zhang `sample_scale`),
sex-split, family orchestrators (GrimAge, SystemsAge), `external_package`, `custom` (MiAge).

**Rejected.** Keeping a "14-verb recipe walker as reference/fallback for remaining
`component_matrices` clocks" (the language that lived in detail-plan sec 2.2 / sec 11). It was
Schrodinger's deliverable -- simultaneously rejected and planned -- and a third scoring path
that would have to agree with the engine and the packs on every fixture.

**Why.** Every remaining clock resolves to either the linear engine or a named pack. A walker
buys a third code path and more fixture surface for zero clocks it alone can score.

**Sites.** `dev/migration-plan.md` (Runtime), `dev/detail-plan.md` sec 2.2, sec 11. Extends
the 2026-07-15 "engines + packs" entry (which already rejected the walker *as the only
runtime*); this removes it entirely.

---

## 2026-07-17 -- result object inherits `list`, not `matrix`

**Decision.** `calc_clocks()` returns an S3 record over `list`, class `"methylCIPHER"`:
`list(scores = <n x k double>, coverage = <per-column>, provenance = <clocks, covariates_used,
batch_set_id>)`. Verbs are methods: `as.matrix` (naked scores / escape hatch), `as.data.frame`,
`[` (subsets scores **and** the matching coverage/provenance), `cbind` (binds columns, checks
`batch_set_id`), `augment` (imported from the `generics` generic, not a home-grown one),
`statistics`, `codebook`, `citation`, `print`.

**Rejected.** `class = c("methylCIPHER", "matrix")` carrying coverage/provenance in attributes.
Base R drops class and attributes on the first `x[, "Horvath"]`, `t()`, or arithmetic, so a
user's first subset silently discards provenance. A bare matrix subclass with load-bearing
attrs is a known R foot-gun.

**Sites.** `dev/detail-plan.md` sec 1.3, sec 4, sec 7.

---

## 2026-07-17 -- imputation source split: partial = cohort, absent = vendor; never cross

**Decision.** Two missing-kinds, two sources, non-interchangeable:

- **Partial NA** (probe present in `DNAm`, `NA` in some samples): impute from the **current
  cohort** (mean over observed samples for that probe). Cohort-dependent by definition -- that
  is what imputation is -- and it is a **fallback**: the contract is that users pre-handle their
  NAs; we only avoid crashing.
- **Completely absent** (probe not in `DNAm` at all): no cohort data exists to borrow, so fill
  from the **vendored** ref the clock specifies (mean / median / sex-wise) or drop per policy.

Crossing them is the bug: vendored-filling a partial-NA injects a foreign cohort's data into the
user's samples; cohort-filling a fully-absent probe is impossible.

**Rejected.** An earlier framing that partial-NA should use the vendor ref "to stay
batch-independent." Batch-independence is not a goal here -- imputation borrows within-cohort
information on purpose. `mean` is merely the current fallback statistic, not a fixed rule.

**Sites.** `dev/detail-plan.md` sec 2.3, sec 4 (coverage splits cohort-borrowed
`score_imputed_partial` from vendor-borrowed `score_imputed_full`), sec 6.

---

## 2026-07-17 -- no SHA / pin provenance in the R package (correctness gated by fixtures)

**Decision.** The package carries **no** commit pin as product identity. `manifest_key` /
`bundle_hash` / `out_sha256` in `data-raw/sync.R` are **sync-internal skip keys only** (answer
"has the catalog changed since last sync?") and never reach `mc_provenance` or a result attr.
Result provenance is: clock ids scored, covariates actually used, coverage, and batch-set id
when applicable. Correctness is proven by **fixtures**, not by a pin.

**Rejected.** Stamping results with the resolved checkout SHA (former detail-plan sec 7 / sec
9.1). Because a no-op sync (`manifest_key` unchanged) leaves `sysdata.rda` lagging the checkout,
a SHA stamp could disagree with the data actually bundled -- provenance that can lie. Dropped
rather than patched.

**Sites.** `dev/detail-plan.md` sec 7, sec 9; `dev/migration-plan.md` (Source of truth,
Packaging).

---

## 2026-07-17 -- fixtures: always-run engine units + cohort-gated parity (no shipped golden)

**Decision.** Two test tiers doing **different jobs**:

- **Engine / machinery unit tests -- always run, no meta dependency.** Hand-authored toy inputs
  (a few probes x a few samples, coefficients chosen so the expected dot product is written by
  hand in the test). They verify `linear_score` arithmetic, impute accounting, accessors,
  coverage math, and result-class methods. Golden values tie to nothing upstream -> zero drift.
- **Parity fixtures -- the single clock-golden source, cohort-gated.** The upstream golden
  fixtures vs `fixtures/cohort_EPIC/beta.duckdb`, skipped via `file.exists()` when the cohort is
  not staged under `data-raw/methylCIPHER-meta`. CI may stage/regenerate the cohort and run the
  full gate; CRAN skips it.

**Rejected.** Committing a small slice of the real cohort's golden values into `tests/` to make
parity run on CRAN. That is a second copy of upstream golden values -> it drifts. A 5-sample
slice also cannot validate a `correlation`-policy clock anyway; only `exact` policy is reliable
at small n, and those do not need a bespoke shipped fixture once engine units cover the wiring.

**Sites.** `dev/detail-plan.md` sec 10.

---

## 2026-07-17 -- sysdata schema lives in accessors + a structural test, not a schema.md

**Decision.** No hand-written `schema.md`. The sysdata shape is defined by `build_index()` /
`build_catalog()` in `data-raw/sync.R`, and the **accessor layer is the executable schema**:
`get_clock()`, `clock_scoring_cpgs()`, `clock_norm_cpgs()`, `clock_impute()`, `clock_coefs()`,
`clock_group_bundle()`. Their `#` comments state the fields/types (roxygen deferred to alpha, per
the 2026-07-21 entry); a `testthat` test asserts `names(mc_index)` and a sample `mc_catalog`
entry's structure, so a shape change breaks a test and forces the comment update. `calc_clocks` code consumes only accessors -- never raw
`entry$components[[i]]$file`. Covariate requirements are one flattened catalog field
(`covariates_required`, computed by sync's `extract_covariates`), read once at runtime, never
re-derived from three meta keys in the scoring path.

**Rejected.** A prose `schema.md` for future agents -- the exact prose<->code drift class the
meta repo's `meta_schema.py` exists to eliminate. If a browsable schema is ever wanted, generate
it from the objects; do not hand-maintain it.

**Sites.** `dev/detail-plan.md` sec 5, sec 8; `data-raw/sync.R` (`build_index`,
`extract_covariates`).

---

## 2026-07-15 -- engines + packs (not pure recipe VM; not 113 calc*)

**Decision.** Runtime is a shared **linear engine**, small **pre-transforms** (e.g. Zhang
`sample_scale` via full-panel moments then coef subset), and a few **family orchestrators**
(GrimAge, SystemsAge), plus external/custom one-offs. Public scorer is `calc_clocks`.
Meta/fixtures remain the scientific contract; orchestrators are optimized implementations.

**Rejected.** (1) Pure “one recipe interpreter walks all 113 clocks” as the only runtime.
(2) Restoring per-clock unexported `calc*` bodies as the engine. (3) Global probe intersect
in `calc_clocks` before scoring. (4) Full-panel size warning on every call — localize to
clocks that need full-probe-panel stats (Zhang today).

**Also.** Pass raw beta into scorers; small subset copies inside scorers are fine. Coverage
is recorded at score time; `statistics()` only formats attrs. GrimAge UX prefers pack /
components, not V1/V2 as primary ids (see detail plan).

**Sites.** `dev/migration-plan.md` (overview), `dev/detail-plan.md` (full design).

---

## 2026-07-15 -- archive legacy `R/*.R` outside the package tree

**Decision.** Snapshot every pre-rewrite `R/*.R` into one reference file
`dev/legacy/R-pre-rewrite.R` with `# ---- <filename> ----` separators, remove those sources
from `R/`, keep `R/sysdata.rda`. Rebuild via `uv run --no-project python dev/legacy/archive_R.py`
(simple: for each `R/*.R`, emit separator, append body). Also snapshot
`NAMESPACE` → `dev/legacy/NAMESPACE-pre-rewrite`.

**Why not wipe without archive.** Legacy calculators still encode edge-case math, imputation
quirks, and API shapes worth grepping while porting verbs. Pure delete forces `git show`
archaeology for every question.

**Why not leave sources in `R/` (even as one mega-file).** `load_all()` / install load all
`R/*.R`. Legacy in `R/` keeps the old export surface and collides with the rewrite.

**Why one file under `dev/legacy/`.** Frozen reference for `rg`; separators recover original
filenames. Not multi-file because we are not editing the legacy tree further.

**How to use.**

```text
rg calcHannum dev/legacy/R-pre-rewrite.R
# rebuild archive from git history only if R/*.R still exist:
# uv run --no-project python dev/legacy/archive_R.py
```

**What stays in `R/` after archive.** `R/sysdata.rda` plus new rewrite sources as they land.

**Related leftovers.** `man/*.Rd` and live `NAMESPACE` still describe old exports until
roxygen is regenerated. Do not reintroduce deleted `R/` sources just for old `.Rd` files.

**Sites.** `dev/legacy/archive_R.py`, `dev/legacy/R-pre-rewrite.R`,
`dev/legacy/NAMESPACE-pre-rewrite`, emptied `R/*.R`, `R/sysdata.rda` retained.
