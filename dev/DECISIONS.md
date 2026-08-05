# Decisions log (package rewrite)

Append-only, date-stamped. Records **why** we chose a design -- the "we tried X / other
maintainers will ask this" history that should not bloat an operational doc.

**Scope:** R package rewrite, packaging, API, and local maintainer workflow. Upstream
metadata contract decisions live in `data-raw/methylCIPHER-meta/control/DECISIONS.md`.

Newest first. Add an entry when a decision reverses a prior approach or is likely to be
second-guessed; do not restate rules already stated in `CLAUDE.md`.

**Archive:** entries before 2026-07-30 live in `dev/DECISIONS.old.md` (unchanged full log).
Older dated citations in `CLAUDE.md` resolve there. Do not restate that history here.

---

## 2026-08-04 -- source-tree-only tests are build-ignored, not CRAN-skipped

The first `devtools::check()` after the trim came back `1 error | 1 warning | 0 notes` in 51s (it
did not hang -- see the tier entry below for why). Both findings had the same root cause and the
same fix, and both falsify something asserted earlier the same day.

**The error.** `test-source-hygiene.R` failed twice in the installed package. The prediction was
that its scans would "pass vacuously over an empty file list" and that `skip_on_cran()` made that
an honest skip. **Both halves were wrong.** `scan_sources()` is `unlist(Map(f, R_PARSED, ...))`,
and `unlist(list())` is `NULL`, not `character(0)` -- so `expect_equal(NULL, character(0))` is a
*failure*, not a vacuous pass. And the gate never fired, because **`devtools::check()` sets
`NOT_CRAN=true`**. A `skip_on_cran()` is inert in exactly the check most likely to be run locally.

**The warning.** `checking for unstated dependencies in 'tests'` flagged `duckdb`, from
`test-fixtures-parity.R`. This one **cannot be fixed by any runtime skip**: the check is a *static
scan of the shipped sources*, so a `duckdb::` behind `if (parity_on)` still counts. `duckdb` and
`DunedinPACE` are deliberately undeclared -- the parity tier is maintainer-gated and CI installs
them itself -- so the only way to satisfy the scan without declaring a dependency the package does
not have is for the file not to be in the tarball. (`DBI` is in `Suggests` and was never at issue.)

**The rule that falls out.** A test whose subject is the **source tree** rather than the package's
behaviour does not ship. `test-fixtures-parity.R` reads `data-raw/methylCIPHER-meta/` (already
build-ignored) and `test-source-hygiene.R` reads `R/*.R` (an installed package ships `R/*.rdb`);
neither can ever work from a tarball, under any `NOT_CRAN` value. Both are now in `.Rbuildignore`,
joining the existing `^R/dev-utils\.R$` -- which is where `test_parity()` lives, so the precedent
was already set and simply had not been followed through to the test files.

Verified by building the tarball and listing it: 28 test files ship, neither of these two among
them. `skip_on_cran()` was then **removed** from the hygiene tests -- with the file build-ignored
the gate is dead weight, and worse, it implies the file ships. They now run on every
`devtools::test()`, which is the dev loop they exist for.

**What this does not change:** `skip_on_cran()` stays the right tool for the ~86 internal tests that
*do* ship and *can* run installed. The distinction is whether the test can execute at all outside
the source tree, not whether it is internal.

---

## 2026-08-04 -- Horvath1 is admitted to parity under a measured BMIQ snapshot

The third tolerance regime that 2026-07-29 declined to decide. Deciding it retires the last
hand-authored clock golden in the always-on tier.

**The reading being tested.** The `horvath` block is skipped wholesale because the oracle filled
completely-absent probes server-side with an unpublished constant, so the residual tracks the
absent-probe count. `Horvath1` was flagged as the exception -- the one clock the oracle BMIQ'd --
whose gap should therefore be a *normalization* gap, fixable by turning normalization on rather
than by widening a tolerance. That was an argument, never a measurement. It is now measured, and
it holds:

| cohort | scoring panel | `normalize` off | `normalize` on |
|---|---|---|---|
| `cohort_450K` | 353/353, zero absent | max_abs **7.72** | max_abs **0.114**, max_rel 1.93e-03 |
| `cohort_EPICv1` | 334/353, 19 absent | max_abs 6.53 | max_abs **3.96**, max_rel 1.71e-01 |

So on the complete panel, normalizing closes 98.5% of the gap and what is left is a
BMIQ-implementation difference. On the 19-probe-short panel it closes almost nothing, because
there the *fill* gap dominates -- exactly the original diagnosis, now with the two effects
separated by data rather than asserted.

**Why the split is a guard, not a second tolerance.** Admitting `cohort_EPICv1` at 4.0 years would
be the vacuous bound 2026-07-25 warned about. The test instead **skips unless the cohort leaves
zero scoring probes absent**, which is the precise condition under which the oracle's undisclosed
input cannot contaminate the comparison. That predicate is computed from the loaded matrix, so no
cohort is named anywhere and a restaged or extended cohort re-decides itself.

**Why a snapshot rather than an agreement target.** 1.9e-03 relative is far too loose to be an
agreement claim and is not offered as one; `PARITY_REL_TOL` stays 1e-10 for everything that is a
real gate. What this pins is *that the residual does not move*, which is the only thing worth
asserting against an oracle whose own pipeline we cannot reproduce. Verified deterministic --
bit-identical `max_abs`, `max_rel` and score checksum across three runs -- so the ceiling sits just
above the measurement instead of being padded, and drift will actually trip it. A pair that clears
the absent-probe guard with no entry in `HORVATH_NORM_TOL` **fails**: a newly admissible pair needs
its residual measured, never defaulted to a neighbour's.

**Membership is derived.** `is_normalized_horvath()` is horvath-online **and**
`clock_norm_scheme() %in% NORM_SCHEMES`. Today that is `Horvath1` alone -- of the 15 horvath-online
clocks, 13 declare `none` and `Horvath2` declares `noob`, which is not expressible here. (CLAUDE.md
previously said "14 of the 15 declare `scheme = none`"; corrected in the same pass.)

**What it replaced.** The `test-normalize.R` BMIQ golden, which composed
`betanorm::bmiq_calibration()` with the linear score by hand. **This is a real trade and worth
naming**: the deleted test was exact and ran on any machine with `betanorm`; the new one is
maintainer-gated and bounded at 1e-03. What is bought is that the gate is now against the
**oracle** rather than against the same library we call, so it can catch a wiring error the
hand-composed golden shared. The *record* half of the old block did not move -- `provenance$
normalized`, `cov$normalizes`, the `sample_miss$norm` column are things parity never looks at, so
they stay in `test-normalize.R` under a name that says so.

Standing parity state goes 264 -> 266 blocks.

---

## 2026-08-04 -- the always-on suite is cut to ~800, and CRAN sees only the exported surface

The suite had reached ~1284 expectations across 27 always-on files and was no longer maintainable:
a routine change touched a dozen files, and the cost of that was not buying proportional safety
because **the parity tier already proves the arithmetic**. Trimmed to 801 local / 391 under
`NOT_CRAN=false`, 0 failures either way, with the parity file untouched.

**The two cuts, and why they are different.** Deletion was applied to tests that could not fail
for a reason anyone would act on: goldens parity already owns (every catalog clock declares a
fixture, so a broken scoring branch fails parity before it fails a unit test), assertions about
maintainer-side plumbing that "Test altitude" already banned, message-wording pins, and loops
restating one invariant over many clocks. Gating with `skip_on_cran()` was applied to what
survives but is **internal** -- the kernel contracts, the catalog accessors, the chunking path,
the asset transfer mechanics. Those still run in dev and CI; they just stop being CRAN's problem.

**What CRAN actually runs now, stated plainly so nobody misreads a green check:** the smoke tier
in the default configuration, the front-door refusals, and the record contract. **No numeric
gate.** That is not a regression -- parity was already `MC_PARITY`-gated and CRAN-skipped, so CRAN
never proved a score. The change is that this is now visible instead of implied by a suite that
looked comprehensive. A CRAN green says the package loads, refuses correctly, and returns a
well-formed `mc_result`; it says nothing about the numbers.

**But `skip_on_cran()` is narrower than it sounds, and the trim should not be read as "these tests
now only run locally".** `NOT_CRAN` is unset on r-hub and on a GitHub Actions `R-CMD-check`, so the
gated tier runs on both -- across platforms, which is where a Windows-encoding or a long-double
difference would actually surface. What the flag buys is CRAN's own machines not paying for a tier
that cannot tell them anything. `devtools::check()` is the opposite case: it sets `NOT_CRAN=true`
and runs the whole suite, so a local check and a tarball check are not the same run.

**Addendum, same day, from the first `devtools::check()` after the trim.** Two files needed
`.Rbuildignore`, not `skip_on_cran()` -- see the entry below.

**Four goldens were kept against the "parity owns it" rule**, because parity is structurally blind
to them and the blindness is not incidental:

- **Alias routing** (`test-score-fitage.R`, `DNAmGrip_wAge`). Fixtures are declared on the 14
  routed *members*, never on the 7 aliases -- so *which sex's model scored which sample* has no
  parity coverage at all.
- **DunedinPACE quantile normalization** (`test-score-dunedin.R`). Since the reference golden moved
  to the parity tier earlier today, this is the only always-on proof normalization is applied.
- **PhysAge mean-divisor fill offset** (`test-score-physage.R`). Fill landing inside vs outside the
  divisor is a silently wrong number on a degraded panel; parity scores clean panels.
- **Wang mixed-request domain isolation** (`test-score-wang.R`). Parity scores one clock per call,
  so `sample_scale` contamination between two clocks in one request cannot appear there.

The BMIQ golden in `test-normalize.R` was kept for the same class of reason, and then **superseded
within the day** -- see the next entry.

**`test-sim-smoke.R` is untouched and ungated.** It is the only tier running `calc_clocks()` in the
default configuration and the only caller of `sim_DNAm()`, which is exactly what a CRAN machine
should be exercising. **`test-source-hygiene.R` is gated, and not because it is internal**: an
installed package ships `R/*.rdb`, not `R/*.R`, so `list.files(..., "\\.R$")` returns nothing there
and both scans would pass vacuously over an empty file list. `skip_on_cran()` converts a fake pass
into an honest skip.

**The gate goes inside the block, never at file level.** testthat runs top-level code at collection,
so a file-level `skip_on_cran()` reads as one skipped file rather than N skipped tests and hides how
much is off. One call, first line of each gated `test_that`; there is no top-level `skip_on_cran()`
in any file. Under `NOT_CRAN=false` the suite reports **89 skipped blocks** spread across 24 files
(the parity tier being one of them), which is how to check the gating is per-block rather than
per-file: a file-level gate would have shown 24.

---

## 2026-08-04 -- the parity tier gates its generator, and the suite runs silent

Two separate complaints about the same thing: the default `devtools::test()` was unreadable.

**The parity wall.** `run_parity_target()` opened with `skip_if_no_cohort()`, so with the tier off
all 263 generated targets ran far enough to skip. testthat prints skips grouped by reason **with a
location per skip**, so one sentence came back as ~90 wrapped lines of
`test-fixtures-parity.R:339:9`, repeated. The reason is identical every time and the locations are
all the same line, so the block carried exactly one bit of information and buried the run summary
under it.

The fix is to gate the **generator** rather than the generated test: `staged_cohorts` (the flag AND
a live duckdb connection) drives `parity_targets()`, the PhysAge loop, the census test and the
Dunedin golden, and a tier that cannot run emits **one** `test_that` saying so, plus one per
unstaged cohort. Verified generation is otherwise untouched -- with both cohorts staged the file
still produces 264 blocks, the same 146/28/30/56 split across `core`/`fitage`/`horvath`/`packs`.

**Why not just quieten the reporter.** Because the skips were never informative individually. The
counterpressure is real and was weighed: the block count is how a parity run is read (CLAUDE.md's
"264 blocks / 32 skip", checked against each other before the pass number), and a generator gate
makes that count depend on what is staged. That is the right dependency -- a test that had no cohort
to read was never a test -- but it means the standing figure is now "with both cohorts staged", and
the guard against a *dropped fixture* has to be the census test, which is already there and already
ungated by cohort. `skip_if_parity_off()` and `skip_if_no_cohort()` are gone; nothing else used them.

**The message noise.** Unrelated in mechanism, same symptom. Three sources, all of them the package
working correctly:

- `mc_spec()`'s full-panel note (`say_full_panel_clocks()`) fires once per spec build, so any test
  touching `Zhang2019EN` printed three lines. It is emitted from `mc_spec()`, **not** only from
  `score_cohort()` -- a test that never calls `calc_clocks()` still triggers it.
- `say_pending()` on a multi-batch `rbind` in `test-bind.R`, at two sites that were not asserting on
  it (the one that asserts uses `expect_message()` and stays).
- `utils::download.file(quiet = FALSE)` in `mc_fetch()`. This one is environment-dependent and was
  the confusing one: with `method = "auto"` a `file://` URL takes the `internal` method and says
  nothing, but an interactive session with `options(download.file.method = "libcurl")` -- what
  RStudio/Positron set -- prints `trying URL` / `Content type` / a progress bar. Those go to
  **stderr as raw text, not as conditions**, so `suppressMessages()` cannot touch them; only
  `capture.output(type = "message")` can. `test-mc_data.R` has a local `quietly()` that does both.
  The package side is unchanged: a progress bar on a 300 MB pack download is worth having.

Same file's one warning came from `download.file()` warning on its way to the failed status that
`mc_fetch()` turns into the abort under test; the test now suppresses it, because it restates the
abort in worse words and is not a second assertion.

**Result: 0 fail, 0 warn, 2 skip on a default run, with no stray output.** Parity was not run.

---

## 2026-08-04 -- `codebook()` is reinstated, and it is blocked upstream

**Reverses the 2026-07-31 decision** that kept D3 out of the finalizer family. That entry rejected
`codebook` on three grounds: it touches no result, it reads a `bib_key` that does not exist, and it
is a third view of `list_clocks()` / `clock_cpgs()`. The first is true and is not disqualifying, the
second named the wrong field, and the third is what changed.

`codebook()` returns `data.frame(clock_id, description)` and dispatches like `cite_clocks()`. It is
not a third view of the catalog: `description` is a sentence per clock saying what the score means,
which is the one column `list_clocks()` does not carry and cannot be derived from anything the
package holds. That is the whole justification, and it stands or falls on the column existing.

**Blocked upstream, and that is where the work is.** `description` is not verified in
`methylCIPHER-meta` across the 137 clocks. The method itself is small. Sourcing and checking one
description per clock is not, and it is upstream work rather than package work. **Do not build
against a partially populated field** -- a `codebook()` that returns `NA` for most of the catalog is
worse than no method, because it looks like a defect in the package.

## 2026-08-04 -- `duckdb` and `DunedinPACE` leave Suggests: a dep is declared for code, not for tests

Both were in `Suggests`, and `DunedinPACE` also needed `danbelsky/DunedinPACE` in `Remotes:`. Both
are gone; `Remotes:` now carries `hhp94/betanorm` alone.

**The line is where the *package* reads the dependency, not whether a test does.** `betanorm`
appears in `R/`: `require_betanorm()` in `R/score_normalized.R` gates every normalizing branch, so a
user who passes `normalize =` needs it at runtime, and it stays declared with its `Remotes:` entry.
Neither `duckdb` nor `DunedinPACE` appears anywhere in `R/`. They are read only by
`tests/testthat/test-fixtures-parity.R` and `tests/testthat/test-score-dunedin.R`, and both tiers
are already gated for other reasons -- parity needs `MC_PARITY=1` plus a staged duckdb cohort no
user has, and the Dunedin golden needs a GitHub-only package. Declaring them made every installing
user resolve a dependency for a tier they cannot run, and in `DunedinPACE`'s case made a CRAN
submission carry a `Remotes:` entry for a package CRAN does not have.

**So the reference golden moved to the tier that matches it.** The degraded-coverage test against
`DunedinPACE::PACEProjector()` was the only third-party-dependent test in the always-on value-golden
tier, where an undeclared dep means it skips silently on every machine forever. It now lives in
`tests/testthat/test-fixtures-parity.R`, behind the tier flag and **with no
`skip_if_not_installed()` on the reference at all**. That is the point: parity only ever runs on a
maintainer machine, which has both `duckdb` and `DunedinPACE`, so a skip there cannot protect
anybody and can only hide a silent non-run. If the reference is missing the test errors, which is
the correct signal on the one machine that runs it. It needs no duckdb and no staged cohort -- it
builds its own holed panel -- so the file's guards split: `skip_if_parity_off()` is the flag alone,
and `skip_if_no_cohort()` calls it and then demands a connection.

The other four tests in `test-score-dunedin.R` are in-package goldens with no third-party dep and
stay in the always-on tier. Measured after the move: always-on 1271 pass / 0 fail, parity 264 blocks
/ 707 pass / 32 skip / 0 fail (was 263 / 699 / 32 / 0 -- one block and eight expectations, which is
exactly the moved test).

**What this costs, stated plainly.** The degraded-coverage path is now checked only when a
maintainer runs parity, and no CI job gets it from DESCRIPTION. Accepted: the alternative is a
CRAN-facing DESCRIPTION advertising deps for tiers CRAN never runs. **A CI job that means to run
that tier must install `duckdb` and `DunedinPACE` itself.**

## 2026-08-04 -- `payload_hash_of()` stops hashing the serialize header, re-addressing every pack once

`payload_hash_of()` hashed `serialize(payload, version = 2L, xdr = TRUE)` as a raw stream.
`serialize()` writes a 14 byte header whose bytes 7 to 10 carry the **writer's R version**, so the
content-address moved on an R upgrade with not one coefficient changed. That is the one thing the
identity key exists to prevent: `payload_hash` sets the pack filename and the release tag, and its
whole job is to make re-upload of unchanged weights a no-op. It is now
`digest::digest(payload, algo = "sha256", serializeVersion = 2L)`, whose `skip = "auto"` drops that
header.

**The fix is not free: it re-addressed all four packs once.** `pcbrainage`, `pcclocks`, `systemsage`
and `zhang2019` each have a new hash, so `R/sysdata.rda` declares four filenames and four release
tags that did not exist before. The assets were re-uploaded under the new tags before this landed,
so no download breaks. Every existing user's cached packs become **superseded** on upgrade, exactly
as any `payload_hash` move does -- `list_mc_assets()` reports them and `clear_mc_assets()` reclaims
them, which is designed behaviour and not a special case. A one-time cost, paid to stop paying it on
every R release.

## 2026-08-04 -- The prose files split from one rule set into two, and stop being hard-wrapped

`dev/WRITING.md` section 10 governed `vignettes/*.Rmd` and `README` with a single rule set. Writing
the README exposed three things wrong with it. The detail is in that file. What matters here is
why each changed.

**The build split is the load-bearing correction.** The old section said an evaluated chunk "must
run offline, with no asset", derived from `R CMD check` building vignettes with no network. That
premise is true of vignettes and **false of `README.Rmd`, which is `.Rbuildignore`d** and is
rendered by the maintainer alone. Nothing on CRAN executes it. Applied as written the rule would
have forbidden the README that now exists, which scores `SystemsAge` against a staged asset and
shows the coverage record doing real work. The cost is that rendering `README.md` needs the assets
staged, so a collaborator with an empty assets directory cannot rebuild it. Accepted: the
alternative is a README that cannot show the package's most interesting output. `README.md` itself
is **not** ignored and does ship, so what lands in it is still held to the rules.

**Prose is no longer hard-wrapped in those two files.** One paragraph is one source line and the
editor soft-wraps. A hard-wrapped paragraph reflows on any edit, so changing one word rewrites five
lines and the real change is lost in the diff. Code blocks stay at 80 columns, because that text is
read as code and nothing reflows it. Scoped to `README.Rmd` and `vignettes/*.Rmd`. **`dev/` docs,
including `WRITING.md` itself, stay hard-wrapped.**

**ASCII binds the author, not the run.** The rule was stated as though it covered the whole file,
which is unachievable: cli emits its own symbols and captured chunk output carries them. Suppressing
them with `options(cli.unicode = FALSE)` was tried in `README.Rmd` and then removed, because it made
the rendered output disagree with the reader's own console, which is the one thing pasted output
exists to demonstrate. So: never hand-write a non-ASCII character, and leave cli's alone. Pandoc is
the third case and needed a fix rather than a rule, since its `smart` extension manufactures curly
apostrophes out of ASCII input. `README.Rmd` disables it with `md_extensions: -smart`.

**`CLAUDE.md` gains an invariant making the read mandatory before any user-facing text is edited**,
widened from the roxygen-only wording it had. The roxygen bullet no longer restates it.

---

## 2026-08-04 -- The assets dir stays `"cache"`, and the CRAN claim that justified it was false

Two separate things, and the first is why this entry exists at all. The rule does not change. Its
justification was wrong, and it was wrong in the direction that forbids a legal option.

**The policy claim was false.** `CLAUDE.md` said the dir must "never" be `which = "data"`, and
`DECISIONS.old.md` (2026-07-24, section 5) called such a change "a policy violation". Re-fetched
from <https://cran.r-project.org/web/packages/policies.html> on 2026-08-04, verbatim and unchanged:

> For R version 4.0 or later (hence a version dependency is required or only conditional use is
> possible), packages may store user-specific data, configuration and cache files in their
> respective user directories obtained from `tools::R_user_dir()`, provided that by default sizes
> are kept as small as possible and the contents are actively managed (including removing outdated
> material).

All three of data, configuration and cache are permitted. Two consequences. The "actively managed"
clause attaches to all three equally, so it is **not** evidence for `"cache"` over `"data"` --
`clear_mc_assets()` is required either way and stays exactly as it is. And "sizes kept as small as
possible" is satisfied by the consent design, not by the directory: nothing downloads unprompted,
so the default size is 0 bytes everywhere. Neither condition discriminates. The archive is not
edited, so this entry is the correction.

**The sizes that framed the question were wrong by roughly 9x.** The corpus is 43.0 MB total
(SystemsAge 22.5, PCClocks 8.9, PCBrainAge 6.7, Zhang2019 5.0), not the "several hundred megabytes"
the question was posed against. That figure came from the hand-written example output in
`vignettes/assets.Rmd`, which shows about 402 MB and does not match what `list_mc_assets()` returns.
Recorded because it moved the argument: at 43 MB a re-download is cheap, which weakens the
durability case that was the only case for `"data"`.

**The decision: keep `"cache"`, on platform merits.** Both are legal, so this is a trade-off, and it
is genuine rather than a clear win.

- Windows is where the hard failure lives. `"data"` is `%APPDATA%`, the roaming half of the profile:
  under a roaming profile it is copied across the LAN at every logon and logoff, and under Folder
  Redirection it lives on a network share so every read during scoring is SMB traffic. Roaming
  quotas around 100 MB are common and 43 MB spends a large slice of one. `"cache"` is
  `%LOCALAPPDATA%`, which never roams and which Storage Sense and Disk Cleanup do not touch.
- macOS is where the soft failure lives. `~/Library/Caches` is Apple's discardable location, sits
  in Time Machine's default exclusions, and is what third-party cleaners target. But there is no
  scheduled purge, and the cost when it happens is one consented re-download of at most 43 MB.
- **Severity is asymmetric, and that decides it.** The `"cache"` cost is consented, one-shot and
  user-initiated. The `"data"` cost on institutional Windows is unconsented, recurring, and paid by
  users who may never score an external clock again.
- A re-download is already routine anyway. Filenames are content-addressed, so every sync that moves
  a `payload_hash` orphans the old pack and every user re-fetches that group. A purge is the same
  event the normal lifecycle already produces; `"data"` would only make it rarer, not novel.

**Durability is `MC_ASSETS_DIR` / `set_mc_assets_dir()`, and the vignette must sell them on that.**
It currently sells them on control. This is the third option from the to-do item, taken deliberately:
the genuinely offline case -- air-gapped node, no-egress cluster -- is already told to use
`ext_data = <path>`, so it is not on the default dir at all. The default's durability matters least
exactly where the offline risk is worst.

**Precedent, measured rather than recalled.** On a machine with both installed,
`BiocFileCache::getBFCOption("CACHE")` and `ExperimentHub::getExperimentHubOption("CACHE")` both
resolve under `R_user_dir(..., "cache")`, and ExperimentHub was holding 1.9 GB there -- 45x this
package's whole corpus. Those are the two general-purpose data caches of the Bioconductor stack,
which is the stack this package's users already run. R's own `R_LIBS_USER` is under `%LOCALAPPDATA%`
on Windows too.

**`torch` is the one CRAN precedent that does otherwise, and it does not transfer.** `inst_path()`
reads `TORCH_HOME`, else falls back to `system.file("", package = "torch")` -- it writes into its own
installed package directory and never calls `R_user_dir`. That is dictated by dynamic linking:
libtorch is shared libraries the loader must find next to the package's compiled code. Our payload
is a `qs2` file read at runtime, with no linking and no loader. Copying torch here would put 43 MB
somewhere `update.packages()` destroys, which is worse durability than the `~/Library/Caches` risk
that opened the question, and it breaks outright on a read-only system library.

**The migration question is dissolved rather than answered.** It only existed if the directory moved.
Nothing moves, no user's assets go invisible, and no detection-or-migrate code gets written.

**A platform-conditional default was considered and rejected.** It buys the macOS property at the
cost of a default path that differs by OS for reasons the user cannot see, plus a branch in code
that has none, to protect 43 MB that re-downloads on request.

**The runtime download path stays as it is**, and the properties it already has are the reason: an
unauthenticated CDN GET rather than a **GitHub API** call, whose anonymous limit of 60 requests/hour
**per IP** is shared across an institutional NAT or an HPC login node; a declared `release_tag` +
`file` pointer resolved by `mc_asset_url()` rather than a list-the-release-then-find-the-file
search, per the standing accessor rule; and `mc_fetch()`'s staged `.part` download,
`validate_checksum` and atomic rename. Any replacement has to keep all four.

---

## 2026-08-03 -- The first vignette, and why almost none of it runs

`vignettes/assets.Rmd` is the package's first vignette. It documents the assets directory, the
four-source resolution order, and the download / load / clear lifecycle.

**Every block that downloads, deletes, prompts, or prints a per-machine path is `eval = FALSE`,
with its output written out by hand.** `R CMD check` builds vignettes, and the standing invariant
is no network at install, build, check or CRAN test. A vignette about downloading is exactly the
document most likely to break that, so the policy is stated in `dev/WRITING.md` section 10 rather
than left to the next author to infer. Exactly one block evaluates, a `list_clocks()` call over
the shipped catalog, which is offline and deterministic.

This is the `@examplesIf interactive()` policy in the other syntax. The two must not drift: an
example and a vignette chunk that show the same call under different guards would be a contract
with itself.

**`vignettes/` and `README` join the R1 to R8 scope.** They were outside it only because neither
existed. Markdown invites a chattier register than a `@details` block, so R2 and R3 are the rules
that slip there, and section 10 says so explicitly.

`DESCRIPTION` gains `VignetteBuilder: knitr` and `knitr` / `rmarkdown` under Suggests. Both are
build-time only and neither reaches `Imports`.

---

## 2026-08-03 -- The finalizer set is derived, and the two batch counts are cross-checked

An independent audit found `score_associations()` calling `finalized()` while the `CLAUDE.md`
invariant named only three finalizers and its own manual said nothing. Three changes came out of
it, and the first is the one that matters.

**The finalizer set is now derived from a mechanical test, not enumerated.** A finalizer is any
exit that takes an `mc_result` and returns something that is not one. That is checkable against any
signature, so the set cannot drift. It drifted twice under the enumeration: `as.matrix()` was
outside it until earlier the same day, and `score_associations()` was re-finalizing silently. The
test also settles `rbind` without a special case -- it returns an `mc_result`, so it is not a
finalizer, which is the same conclusion the old recursion argument reached the long way round.

**The two batch counts are cross-checked rather than chosen between.** `provenance[[mc_batch_id]]`
and `names(per_clock)` are independent derivations of "how many batches", and four sites read the
`per_clock` one while the four exit frames read the provenance one. They agree today. `n_batches()`
(`R/mc_result.R`) now derives the count from provenance -- authoritative, because it is the vector
that fills the column -- and `stop()`s if `per_clock` disagrees. Every finalizer and both coverage
frames route through it. **A record whose counts disagree is malformed, so the right behaviour is
to stop, not to pick one.** Picking one is how a disagreement becomes wrong numbers instead of an
error.

**`finalized()` calls it before the `pending` test, deliberately.** The obvious form,
`if (length(pending) && n_batches(x) > 1L)`, short-circuits: a record with no `pending` -- the
common case -- never evaluates the guard, so `as.data.frame()` and `as.matrix()` accepted a
malformed record while both coverage frames rejected it. Measured, not reasoned about: the first
version of the test failed on exactly those two exits. The check is unconditional now.

**`set_mc_assets_dir()` keeps returning the previous override, `NULL` included. The bug was the
documentation.** The audit flagged `@returns A string.` against a function that returns
`getOption("mc.assets_dir")`, which is `NULL` in the default state. That was first "fixed" by
returning the resolved directory instead, and then reversed the same day, because the resolved
value is the wrong thing for a setter to hand back.

**A setter returns whatever restores the previous state exactly.** `setwd()`, `options()`, `par()`
and Python's `ContextVar.set()` all work this way, and the property they share is that "was unset"
survives the round trip as a distinguishable state. Resolving it away conflates *no override,
falling through to `MC_ASSETS_DIR`* with *pinned to that path*, so restoring the return value
silently shadows the environment variable. It also duplicates the getter: "which directory is in
effect" is exactly what `get_mc_assets_dir()` answers, and it is always a string.

So the split is the ordinary one. **The getter returns the effective value, always a path. The
setter returns the previous override, `NULL` when there was none.** A test now pins the case that
distinguishes them: with `MC_ASSETS_DIR` set and no option, set-then-restore must leave the
environment variable back in charge.

---

## 2026-08-03 -- The `@seealso` groups, and why one generic stopped being four topics

`dev/WRITING.md` section 6 banned `@seealso` outright until the groups could be decided once with
the whole surface in view. That pass is done, and the rule is now "the groups are closed" rather
than "never write one".

**The mechanism is a set, not a per-topic judgement.** A topic's `@seealso` is the **union of the
groups it belongs to, minus itself**. Symmetry is then true by construction, including where a
topic is in two groups, which is the case that would otherwise rot. Five groups, 17 tagged topics,
58 links. `list_mc_assets` is in both discovery and assets and carries the union at eight links --
the widest list in the manual, and intended: it genuinely answers both questions.

**`calc_clocks` points nowhere on purpose.** It is the entry point, every relevant topic already
reaches it through the inherited `x` param text ("The value returned by `[calc_clocks()]`"), and a
hub that links to everything is a table of contents in the wrong place. `sim_DNAm`, `predict_sex`
and the three `print` methods are likewise untagged by decision, not by oversight.

**Four `cite_clocks` topics became one, and that replaced a group rather than joining one.**
`cite_clocks`, `.character`, `.mc_result` and `.default` are all `(x, ...)`, so `@rdname` merges
them with no param collision, and `?cite_clocks` now shows all four usage lines. Co-location beats
cross-reference here: the methods of one generic are the same verb, and a reader wants them on one
page, not linked from four.

**It also closed the footgun `WRITING.md` section 5 calls the one live risk in the donor scheme.**
The donor's `x` is an `mc_result`, so `cite_clocks.character` had to be kept from inheriting it --
inheritance matches on name alone and yields confidently wrong text rather than an error. Merged,
`x` is written once, locally, covering both accepted forms. `lint_roxygen()` went from 6 rows to 0
as a side effect, because the merge retired `cite_clocks.default`'s untyped `Any object.` /
`Nothing.` pair and three duplicated `@returns`.

**The `mc_result` methods were deliberately not merged.** `rbind.mc_result` cannot join them: a
merged topic has one `@param ...` slot, and its `...` means "two or more mc_result objects" while
the others mean "not used" -- a direct collision with the three fixed `...` sentences in
`WRITING.md` section 4. `print.mc_sim` cannot merge into `sim_DNAm` either, because `n` is the
sample count in one and the rows to display in the other. Same name, one slot, two meanings.

**A closed set needs an instrument, so `lint_seealso()` joins `R/dev-utils.R`.** It catches the
two failures nothing else does: a link whose topic does not exist, which is an `R CMD check`
WARNING and therefore invisible until the check that is deferred to section 3 of the to-do finally
runs, and a one-way link, which no tool has ever caught because Rd has no notion of a reciprocal
link. Both were verified against injected faults, not just against a clean tree.

---

## 2026-08-03 -- A frame column keyed on a record fact, and `all_columns` as the escape hatch

`clocks_coverage()` returned 17 columns and `list_clocks()` returned 10. Both now return a narrow
default and take `all_columns = FALSE`. The two use **different mechanisms**, and that difference is
the decision.

**`clocks_coverage()` keys on a declared record fact.** Nine columns are always there --
`clock_id`, `group_id`, `policy`, and the six `score_*` counts. Four appear only where they say
something: `role` when the record holds a routing target, `normalizes` plus the five `norm_*`
counts when a clock normalizes, `missing_cpgs` when a CpG is absent, and `mc_batch_id` under its
own existing rule. Its own `@examples` block dropped from 17 columns to 9, and the nine it lost
were all-`FALSE`, all-zero, or empty.

**The norm block is not "usually" empty, it is structurally empty.** Three of 137 catalog clocks
normalize -- `DunedinPACE` on quantile, two on `bmiq`, which is opt-in and off by default. So six
columns and a flag were dead weight in nearly every run the package will ever serve.

**The precedent that decides the mechanism is `samples_coverage()`, which was already conditional
in the row direction.** A `panel = "norm"` row exists only when the clock normalizes. Keying the
`clocks_coverage()` columns on the same fact makes the two frames agree about when the norm axis
exists, instead of one frame being conditional and its neighbour always-on. `samples_coverage()`
itself is untouched and gets no `all_columns`: it has nothing to drop, and a no-op argument for
symmetry is worse than the asymmetry. Its `panel` column stays even when constant, because
dropping it would change the row key from (sample, clock, panel) to (sample, clock).

**`normalizes` travels with the block it describes rather than staying always-on.** `CLAUDE.md`
calls it the one declared panel fact and forbids re-deriving it from `norm_needed`. That still
holds: the column is present whenever the answer is non-trivial, and where it is absent the fact
is simply that nothing normalized. No reader is pushed toward the derivation the rule bans.

**`list_clocks()` cannot use the same mechanism, so it gets a fixed set.** Its content depends on
the user's filter, not on a declared fact, and a generic "drop the constant columns" rule would
delete `group_id` from `list_clocks(pattern = "^Horvath")`. The default is `clock_id`, `group_id`,
`request_as`, `covariates`, `external`, `tags`. `callable` goes because it is **exactly**
`request_as == clock_id` -- measured, both split the same 14 of 137 rows, and there is now a test
asserting the identity rather than a comment claiming it. `group_size` goes because the frame
already carries what it counts, `batch_dependent` and `normalize` because they are set on 2 and 3
of 137 rows.

**The flag does not touch `mc_batch_id`, deliberately.** Folding it in would give one uniform story,
but `as.data.frame()` and `calc_accel()` have no such flag, so `all_columns = TRUE` on a
single-batch record would make the coverage frames carry a join key the other two exits do not.
That reopens the four-exits-together rule of 2026-08-03 to solve a cosmetic problem. The batch
label keeps its own gate.

**The cost is accepted, not overlooked.** A conditional column means `cov[["norm_needed"]]` returns
`NULL` rather than erroring, which is the silent-`NULL` hazard this package bans `$` over. This is
the second and third such axis after `mc_batch_id`. `all_columns = TRUE` exists precisely so code
that names a column directly has a fixed schema to name it in, and both `@details` sections say so.

---

## 2026-08-03 -- An empty `groups` selects nothing, and `"all"` is the only way to say all

`mc_resolve_groups(NULL)` and `mc_resolve_groups(character(0))` used to return **every** external
group. They now return `character(0)`, and `load_mc_assets()` routes through the same function
instead of normalizing `groups` itself.

**The two readings were already in conflict, in one file.** `pack_groups_needed()` returns
`character(0)` on any run that requests no external clock, which is the ordinary case, so
`load_mc_assets(character(0))` had to mean "load nothing" and did. `mc_resolve_groups(character(0))`
meant "every group". The same value meant opposite things eight functions apart, and the obvious
future cleanup -- route `load_mc_assets()` through `mc_resolve_groups()` -- would have made every
`calc_clocks("Hannum")` call try to download all four asset groups. That cleanup is now done, and
it is safe because the semantics were unified first.

**The failure modes are asymmetric, which decides the direction.** These verbs download hundreds of
megabytes and delete files. `clear_mc_assets(x)` where `x` came out of a filter that matched nothing
deleted the entire assets directory under the old rule. Empty-means-nothing fails as a no-op that
the caller notices immediately and harmlessly; empty-means-everything fails expensively and
silently. Nothing is lost, because the default is already `"all"` -- a caller who wants every group
omits the argument or spells it.

`NULL` is folded into the same rule rather than kept as "unspecified". With a default of `"all"`,
`NULL` never arrives from omission, so it is always either an explicit choice or a computed value
that came out empty, and the computed case is the one worth protecting.

**The validation is `checkmate::assert_subset()`, not `match.arg(several.ok = TRUE)`.** This was
tried both ways. Measured, `match.arg(several.ok = TRUE)` on these choices:

- `character(0)` errors, which is the one property it has that we wanted;
- `NULL` returns `choices[1L]` silently, which here is `"all"` -- exactly the mass-download
  behavior this entry removes, reintroduced through a back door;
- `c("SystemsAge", "Nope")` returns `"SystemsAge"`. It only errors when **every** element fails,
  so a typo beside a valid group is dropped without a word;
- it does not deduplicate.

So it gives one of the four properties and silently breaks two. `assert_subset()` gives exact
matching, an error naming both the bad element and the valid set, and `empty.ok` as an explicit
flag, and it is what `list_clocks(tag =)` already uses for the same shape of argument.

Partial matching is dropped with it, deliberately. The group set grows with each sync, so an
abbreviation that resolves today can become ambiguous later and break calling code that worked.
That is the same reasoning behind the blanket `$` ban in `CLAUDE.md`: a convenience that resolves
silently to something the caller did not name.

## 2026-08-03 -- `as.matrix()` is a finalizer, and a finalizer is defined by the exit

`as.matrix.mc_result()` was `x[["scores"]]` and nothing else. It now calls `finalized()` like
`as.data.frame()` and `calc_accel()` do.

Found while documenting the method, not while testing it: writing the `@details` paragraph meant
claiming what the function does with `pending`, and the claim was not true of this one.

**The rule was drawn at the wrong place.** The old set was "the two exits that return a frame",
which is a fact about the return type and has nothing to do with why re-finalizing is necessary.
The reason is that the caller is leaving the `mc_result` structure: `pending` lives in
`$provenance`, a bare matrix has nowhere to carry it, and the value cannot be recovered from. That
is equally true of `as.matrix()`. So the definition is now **any exit that leaves the S3
structure**, and the frame-vs-matrix distinction is not part of it.

The bug this closes is small but real: on a multi-batch record holding a non-empty `pending`,
`as.matrix(res)` and `as.data.frame(res)` returned **different numbers** for the same clock, and
neither said anything. Both now reduce over every sample and both announce it under
`say_pending()`'s unchanged guard.

`rbind` is still excluded, for its own unchanged reason: it recurses under `do.call()`, and it
hands back an `mc_result` rather than leaving the structure at all.

## 2026-08-03 -- The public message surface: audience, not transport, and one English

Four agents audited all 70 user-visible message sites against a proposed rule set. 11 were clean.
The rules and the per-site evidence were in `dev/cli-audit.md` and `dev/cli-audit/` (local only).
`dev/cli-audit.md` was retired on 2026-08-03, once the rules it argued for became `CLAUDE.md`
invariants and `dev/WRITING.md`. Only the four per-site part files remain.
`CLAUDE.md` now carries R1 to R7; this entry records the reversals and the things that will be
second-guessed.

**The cli keep-set was an enumeration and is now a rule.** cli was "front-door only", named file by
file. That put `resolve_normalize()`'s four aborts on the `stop()` side, even though every one of
them is about a `calc_clocks()` argument the caller typed. Those four could not carry markup, so
they could not satisfy R6, so a whole class of user-facing messages was structurally exempt from
the rules. That is a two-tier system, not a rule. The line is now **user-choosable input against
package defect**, which leaves the old keep set almost unchanged in practice and fixes the one
place where the transport was deciding the register.

Worth recording because it looks like a loss: the plain `stop()` text was **better** than its cli
neighbours on R2, with no first person and no "please". The voice problem was never caused by cli.
The bullet grammar invites it.

**`--` and `;` are banned in public-facing text, and that is not a reversal of the ASCII section.**
The ban is scoped to what a user reads -- cli message text now, roxygen prose later. Comments,
dev-facing `stop()` text and `dev/` docs keep `--`, which the ASCII section still requires. The
reason is accessibility: the maintainer has dyslexia that makes `--` and em-dashes require a guess
at the intended meaning, where a single `-` does not. Evidence that the ban costs nothing: all 14
banned characters across the audit were sentence boundaries in disguise, and not one rewrite needed
the ` - ` allowance that was left open for them.

**"Scoring continues" is deleted, not reworded.** Several warnings spent a whole bullet saying the
run had not stopped, which is what a warning already means. The premise behind keeping it was also
wrong: we cannot promise a transposed `DNAm` will fail at the coverage gate, and scoring can fail
for reasons that have nothing to do with the diagnosis. So the messages describe the problem and
name an instrument instead. One casualty is honest: `coverage_gates.R`'s marginal warning carried a
real fact, that more of the panel gets filled by imputation, and restating it accurately needs a
per-clock policy lookup (`vendor_mean` fills, every other policy drops). It was dropped rather than
stated approximately.

**No "... and N more" tail.** The cap was already 10 in four hand-rolled spellings with two
different tails. Less is more: the true total is already in the lead line, so the tail was
restating a number the reader has. The one exception is the interactive delete-consent manifest,
which counts the remainder on its own line -- there the list is what is being consented to. This
narrows `CLAUDE.md`'s older promise that the delete prompt "lists every file", knowingly: the
counts above the list stay exact, and an unbounded render is the failure the cap exists to prevent.

**`sprintf` never feeds cli, and interpolating is not enough on its own.** cli parses message
elements and bullets as templates, so a built string carrying a `{` is read as syntax --
`mc_manifest_bullets()` genuinely aborted with `Could not evaluate cli {} expression` on a brace in
a file name. Building with `format_inline()` fixes the *inputs* but not the output, which still
goes back through cli as a template, so `bullets()` escapes braces at the one door every bullet
passes through. The tempting shortcut, `cli_vec(vec-trunc)`, was rejected: it truncates the display
while still handing cli the whole vector, so it looks like the fix and does not meet the
requirement.

**`capped_bullets()` caps before it formats.** The old order formatted every element then capped,
which is why capping and markup were mutually exclusive: a finished `sprintf` line cannot carry
markup, and marking up an unbounded vector is what the cap forbids. Reversing the order dissolves
the conflict and is also strictly less work.

**`gate_label()` marks up the token, not the whole label.** It renders `"DNAmFitAge" (female
model)`, so the quotes sit on the part the caller can type back. The audit proposed embedding
`sprintf("{.val %s}", id)`, which is exactly the `sprintf`-into-cli mixing banned above; it uses
`format_inline()` instead. Two assertions in `test-coverage-gate.R` tested the exact string shape,
which the altitude rule discourages, and now assert what the label must and must not name.

**Two class names got `{.cls}` and nothing else changed at those sites.** Worth noting because it
is the one rule that reaches sites the other six leave alone: R6 caught `class(x)[[1L]]` rendering
bare in prose at two sites that were otherwise clean.

## 2026-08-03 -- The Age units gate is per row, not distributional, and it lives in `check_pheno()`

`check_pheno()` bounded `Age` only by `assert_numeric(finite = TRUE)`, which accepts an age in
months, weeks or days without comment. It now also warns per row, on both sides:
`AGE_MAX_YEARS <- 122` (the verified human maximum) and `AGE_MIN_YEARS <- -2`.

**Per row, not a cohort statistic.** The queued item floated keying on the distribution -- a median
far past any plausible age -- as the way to separate "wrong units" from "one unusual subject". That
is the weaker test and was rejected. A cohort statistic can only fire once **most** of the cohort is
wrong, so it is blind to the row-wise failure: a degenerate upstream `ifelse`, or a merge of two
files that coded age differently. That case leaves every other row looking fine, which is exactly
when a human will not catch it by eye. The per-row test also catches everything the distributional
test would have -- if the median is past 122 then rows are too -- so it strictly dominates. It is a
`warn`, so the usual objection to a sensitive test (it stops good runs) does not apply.

**Both bounds are at the edge of possibility, not the edge of typical**, because a units error is an
order-of-magnitude error: a 50-year-old is 600 in months and 18000 in days. So `122` never doubts a
real centenarian, and `-2` leaves the legitimate pre-birth convention (`-0.5`, `-1`) alone while
still catching gestational age in weeks (`-40` to `0`). The two sides warn **independently**,
mirroring `check_col_values()`'s `min_val` / `max_val` -- one pheno can carry both.

**It lives in `check_pheno()`, and the alternative would have covered half the surface.** The
obvious placement is after `resolve_pheno()`, which is where the canonical id vector and the
narrowed covariates first exist together. But `calc_accel()` **never calls `resolve_pheno()`** -- it
merges `data =` over the record's `$pheno` in `merge_accel_data()` and does its own id-join. It
does call `check_pheno()`, with `extra_columns = vars` (which includes `"Age"` whenever
`type = "diff"` or the formula names it) and `sample_id` in hand. So `check_pheno()` is already the
one place both entry points meet, and `warn_missing_covariates()` already establishes the pattern
there: narrow to the join-surviving rows, then warn per column. Both warners now share
`joined_rows()` rather than each re-deriving the join -- one fewer hand-rolled left join, not one
more (cf. the `collapse` question).

**Gating on `extra_columns` means the check fires exactly when the age is consumed.** If no
requested clock requires `Age`, `resolve_pheno()` drops it from `$pheno` and nothing warns -- the
package does not audit columns it does not read. If that same record later reaches
`calc_accel(type = "diff")`, the age must arrive via `data =`, `vars` contains `"Age"`, and the
gate fires there instead. There is no path where an age is used unchecked.

**The flagged ids are returned, not stored.** `warn_age_units()` hands back the offending ids on the
`sample_id` axis, so it is testable and could be threaded into `$provenance` later. It deliberately
is not: the flag is a pure function of `$pheno$Age`, which the record already carries, so storing it
duplicates derivable state -- the same reasoning that refused a `below_min` column on the coverage
frames. What would change the answer is a case where the age is **not** recoverable from the record,
and the one candidate (`calc_accel(data =)`) has nowhere to store it anyway: `calc_accel()` returns
a frame, not a record.

Suite: 0 failed / 0 error / 0 warning / 264 skipped / 1242 passed (+8, five new tests).
`R CMD check` not run. Parity not run.

---

## 2026-08-03 -- Assertions live at the boundary, and `.var.name` is filled only where the deparse lies

One pass over `R/`: 32 `checkmate` calls across 12 files -> 30 across 9. `coverage_gates.R`,
`missingness.R` and `score_cohort.R` now use `checkmate` not at all, which is the point -- all
three were asserting values that had already crossed the front door.

**What forced the pass was a measured message, not a principle.**
`calc_clocks(min_clocks_coverage = "a")` reported

```
Assertion on 'threshold' failed: Must be of type 'number', not 'character'.
```

because that floor was asserted **nowhere** at the front door -- `calc_clocks()` asserted only
`min_samples_coverage`, and the floor's sole validation was `assert_number(threshold, ...)` three
frames down in `check_coverage()`, reached through `mc_cohort()`. So the assertion named a variable
the caller never typed, while the cli messages in that same function correctly said
`{.arg min_clocks_coverage}`.

**The tempting fix is the wrong one.** Setting `.var.name = "min_clocks_coverage"` inside
`check_coverage()` makes the message right for today's one caller and a lie for any other, and it
buries a missing boundary check under a cosmetic patch. The floor is now asserted in
`calc_clocks()` beside its sibling, and `check_coverage()` asserts nothing -- at which point the
deparse is correct for free. **That is the general shape: a wrong variable name in a `checkmate`
message is usually evidence the check is in the wrong frame.** Look there first.

**So `.var.name` is filled only where relocation cannot help.** `checkmate` derives the name by
deparsing the expression, so at a boundary the default is already the caller's own word and a
hand-written string is pure staleness risk -- it survives the next rename; the deparse does not
need to. The rule: fill it **iff the deparsed expression does not name something the caller can
locate in their own call.** This is deliberately *not* "iff it is not a bare symbol" --
`assert_character(colnames(DNAm))` deparses to `colnames(DNAm)`, which any caller can find, and
filling it would be noise. The one site in this pass that qualified was `check_pheno()`'s
`assert_character(pheno[[ID]], ...)`, where `ID` is an internal parameter name; it now carries
`.var.name = paste0("pheno$", ID)`, which names the *actual column* and so beats both the internal
name and a generic `"pheno_id"`. Precedents already in the tree: `predict_sex.R`'s
`"pheno$Female"` and `missingness.R`'s `sprintf("moment_sets[[%s]]", who)`.

**`check_moment_sets()` went the other way -- checkmate to bare `stop()`.** Its input is
catalog-derived (`resolve_moment_domains()`), so a failure is a package bug and a `checkmate`
message aimed at a user is the wrong register. It is *not* deleted, because it is the guard between
a bad index and an out-of-bounds kernel read; the bounds, NA, integer and arity checks are all
still there, just as `stop()` with greppable text. Its existing tests assert *that* it errors and
nothing about wording, so they carried over untouched -- which is the altitude rule paying for
itself.

**`check_pheno()`'s missing-id-column refusal became cli.** `assert_choice(ID, names(pheno))`
printed `Assertion on 'ID' failed: Must be element of set {'zz'}, but is 'ID'` -- where `'ID'` is
simultaneously the fake variable name and the value. Pheno structure at the `calc_clocks` front
door is already on the cli keep-list, so this was never a `.var.name` question.

**`check_mc_result()` is a front-door refusal.** This was the open question section 4.1 flagged: it
is not an S3 method, so the keep-list did not name it. All six call sites (`rbind.mc_result`,
`refinalize_clocks`, `as.data.frame.mc_result`, `calc_accel`, `clocks_coverage`,
`samples_coverage`, `score_associations`) are exported verbs receiving a user-supplied first
argument, which is the keep-list's own pattern. It is now `cli_abort` and reports the class it got.

**The suite did not move: 0 failed / 0 error / 0 warning / 264 skipped / 1234 passed, identical to
the pre-pass baseline.** That is the check on the pass, not a coincidence -- tests assert *that* an
error fires rather than its wording, so relocating and rewording checks is invisible to them. A
failure here would have meant a test was too tight, not that the pass was wrong.

`R CMD check` not run (maintainer-only). Parity not run.

---

## 2026-08-03 -- `predict_sex()` compares against the recorded sex; the join is by id

`predict_sex()` already took `DNAm` and `pheno` on one call over one matrix and passed `pheno`
through purely for the id column. It now also reads `Female` off it and returns two companion
columns beside the declared `predicted_sex`: `recorded_sex` and `sex_mismatch`.

**The recorded sex is read from the caller's `pheno`, not from the record, and that is forced.**
`resolve_pheno()` narrows `$pheno` to the id column plus the covariates the run *required*, and the
two `DNAmSex_Wang` members declare none -- so `Female` never reaches `$pheno` on a `predict_sex()`
run. There is no version of this that reads the record.

**So it is a left join, and the join key is the id column.** Not row order, and not row names --
`match(out[[pheno_id]], pheno[[pheno_id]])`, with `pheno_id` read back from
`$provenance$pheno_id` rather than re-derived, because `...` may have carried a non-default one
into `calc_clocks()`. The match is total and one-to-one for a reason worth writing down rather than
re-checking: `calc_clocks()` has already run, and `check_pheno()` refuses a duplicated or NA id
column while `resolve_pheno()` refuses a pheno that does not cover every `DNAm` row. An `NA` from
the match is therefore a package bug and is raised as one, not handled.

**`Female` is validated here because nothing else does.** It is not a required covariate for these
clocks, so `check_pheno()`'s `assert_integerish(lower = 0, upper = 1)` never fires on it. The same
assertion runs in `recorded_from_female()`, so a factor or an `"M"`/`"F"` column is refused with
the message it would have got from `calc_clocks()`, rather than being guessed at.

**Only an unambiguous binary disagreement is flagged.** The rule table emits `47,XXY` and `45,XO`,
so a flat mismatch would fire on biology; those are shown against the record and never flagged, as
are an unscored sample (NA call) and an unrecorded one (NA `Female`). `BINARY_CALLS` is checked
against the declared `karyotype_calls(kc)` on every run, so an upstream rename fails loudly instead
of silently never flagging. The flag is an invitation to investigate -- either side can be the
wrong one -- and `say_mismatch()` says so in as many words.

**Not done, and not to be re-proposed: auto-resolving sex inside `calc_clocks()`** when a requested
clock needs `Female`. (a) Sex-chromosome probes are routinely filtered out -- `cohort_450K` has
none, which is why both `DNAmSex_Wang_*@cohort_450K` sit in `KNOWN_PARITY_GAPS` -- so an implicit
check would be unavailable on an unpredictable fraction of matrices, and an explicit surface can
refuse where an implicit one can only shrug. (b) It would put the sex panel under the coverage gate
on runs that do not need it, turning a working `DNAmFitAge` call into a hard `ratio == 0` stop on
any sex-filtered matrix -- avoidable only by exempting one clock from the gate, a special case in a
system whose pitch is that routing is total. (c) It is quality-of-life, not a correctness guard.
The residual gap -- the user who never calls `predict_sex` -- belongs in the docs for the
sex-requiring clocks, not in a runtime hint that would fire on every GrimAge and FitAge run.

---

## 2026-08-03 -- `col_stats()` carries the observed range, not two booleans

`any_lt0` / `any_gt1` are gone. The kernel now carries `min_val` / `max_val`, **seeded at 0.0 and
1.0** so only an out-of-range value ever moves them, plus `min_col` / `max_col` -- the 1-based
position within `cols`, the same convention `overflow_col` already uses, so R can name the probe.
The old flags are derived (`min_val < 0`, `max_val > 1`) and nothing is lost.

Cost is unchanged: the same two comparisons on the same branch, two assignments on the rare side,
no extra pass and no new pass structure. The `else if` stays correct because `min_val <= 0 <= 1 <=
max_val` holds by construction, so the two branches remain exclusive. The kernel is serial -- no
OpenMP anywhere in `src/` -- so the reference-accumulation pattern the booleans used carries over
with no reduction; the four fields are one struct only to keep the signatures readable.

**Only the panel columns are scanned, and that is enough.** Scale is a whole-matrix property, so
the panel is a valid sample of it. Pass 2 (`sweep_moments_remaining`) still tracks nothing, exactly
as before.

The point of the change is the message. The `< 0` warning already named M-values and gave the
conversion; the `> 1` one only said "double-check the scale". Both now report the observed extreme
(at `signif(, 4)` -- a diagnostic, not a value to compute with) and the column it came from, and
above **50** the `> 1` warning names percent methylation and gives `DNAm / 100`. 50 is not a
guess at a boundary: a units error is an order-of-magnitude error, so anything between 1 and 50 is
not a scale story and gets the honest "check the scale" instead.

---

## 2026-08-03 -- The coverage floors go into `$provenance`, batch-wise

`min_clocks_coverage` and `min_samples_coverage` were arguments that reached no result field, so a
saved record could not say what floor it was scored under and the `check_row_coverage()` warning
was unreproducible from the record. Both now sit in `$provenance`, **keyed by batch label** like
`$coverage$per_clock`.

**Batch-wise, not record-wise, and `rbind` reconciles nothing** -- the same line every other bind
gate follows: record what batching forced, refuse what the caller chose differently. A differing
floor across batches is not an error and does not throw. It cannot be: the raw per-sample
`coverage` column survives the bind, so which samples in which batch sat low stays answerable
whatever the floors were.

**The exits finalize it, the way the finalizers resolve `pending`.** `samples_coverage()` takes the
**most restrictive** floor across the bound batches (`max`) and re-warns under it. That is what
makes a post-`rbind` `sc[["coverage"]] < threshold` filter well defined at all -- without it there
is no single number to compare against.

**No `below_min` logical column.** The cell axis already exists -- one row per (sample, clock,
panel) carrying `coverage` -- and a conditional second column would mean different things on
different rows after a bind, which is precisely the failure the batch-keyed storage avoids.

Asymmetry between the two is deliberate and follows from what each gate does. `min_clocks_coverage`
aborts, so a record's existence already proves it passed and nothing reads it back; it is stored
for the record's own account of itself. `min_samples_coverage` only warns, so the raw per-sample
numbers are what matter -- and those were already carried.

Resolved while scoping, recorded so it is not re-derived: NA in `$scores` has exactly four sources
-- NA covariates, unknown sex on a routed alias (both recoverable from the retained `$pheno`), BMIQ
unfit samples and Wang domain failures (both in `$provenance$scoring_failures`). None come from the
beta matrix. Unverified corner: MiAge L-BFGS-B non-convergence.

---

## 2026-08-03 -- The test-suite trim is the last step before public alpha, not the next one

The queued audit of the suite (assert output not wiring, let parity own the goldens it already
owns, delete what a no-behavior refactor would break) is **deferred until immediately before
public alpha**, and it is the last piece of work before it. It was previously slotted first, on
the reasoning that a faster suite makes every later iteration faster.

**Why it moves to last.** The suite is the thing that has to be stable when `R CMD check` starts
running, and check does not run yet -- it is maintainer-on-demand precisely because the suite is
bloated. So the trim and the first real check run are one piece of work, and doing the trim now
buys a faster suite for a period in which the code it tests is still moving. Two things in
particular are still unsettled and both rewrite test surface: **roxygen prose does not exist**
(the exported surface still carries a placeholder block, and turning prose on is its own pending
decision), and the validation/message pass will relocate and reword checks across `R/`. Trimming
against a surface that is about to move means auditing the same files twice.

**What this is not.** Not a reversal of the altitude rules -- they stand and bind on every test
written meanwhile, which is what keeps the eventual trim from growing. Not permission to add
tests loosely on the theory that a cleanup is coming. And not a change to the standing
prohibition on running `R CMD check`: it stays maintainer-on-demand until the trim lands, which
is the point at which running it is expected to become routine.

**Ordering consequence.** The message/validation pass (validate at the boundary, then tighten the
bounds, then tone) is a deliberate forcing function for the trim: those passes reword and
relocate messages, and by the altitude rule tests assert *that* an error fires rather than its
wording, so anything that breaks under them was too tight and is already on the trim's list.
Doing them first means arriving at the trim with the list half-written.

---

## 2026-08-03 -- `check_DNAm()` diagnoses shape and replicate probes; it never dedups

Three changes to `check_DNAm()`, and one standing refusal.

**Orientation moved ahead of the matrix refusal, and now runs on a data.frame too.** `dim()`,
`colnames()` and `rownames()` all work on a data.frame, so the orientation and replicate checks are
computed before `is.data.frame()` throws. A data.frame caller therefore gets the actual problem --
and the throw itself adapts, offering `t(as.matrix(DNAm))` rather than `as.matrix(DNAm)` when the
probe ids are in the rows. Previously they got only "must be a matrix" and had to discover the
orientation trap on the next call.

**`nrow > ncol` is evidence, not a verdict.** The obvious cheap orientation test false-positives on
a real run: 1000 samples x 353 CpGs, one clock over a large cohort, is correctly oriented and
taller than it is wide. So probe ids in the rows is the decisive signal, and the dimension ratio
only warns when the columns *also* fail to look like probe ids. Pinned by a test.

**EPICv2/MSA replicate suffixes are now named.** Both arrays suffix every probe with its address
(`cg00002033_TC11`) and ship several rows per CpG -- 8523 duplicated stems on MSA (up to 8 copies),
5225 on EPICv2 (up to 10). Panels are declared on the unsuffixed id, so every such column matches
nothing and is silently vendor-filled or dropped as absent. `clocks_coverage()` already surfaced
this as elevated `score_absent`; what was missing was the diagnosis.

`PROBE_REPLICATE_SUFFIX` is `_[BT][CO][0-9]+$`. Measured against the manifests rather than guessed:
it matches **all** 937055 EPICv2 and **all** 281806 MSA ids, and **zero** EPICv1 (865918) or 450K
(485577) ids, which carry no underscore at all. `[0-9]+` and not `{2}` because EPICv2 has exactly
one five-character suffix, `cg06373096_TC110`. The package's own panels are unaffected: 3 of 452499
panel ids contain an underscore (`ch.13.39564907R_II_R_O_37491` and two siblings, from the
Retroelement clocks) and none match.

**The scans are a bounded stride sample, which replaces the old `ncol < 2e5` guard.** That guard
was sound for orientation -- 2e5 columns cannot be samples -- but it points the wrong way for this
check, because EPICv2 is 937k columns wide, so any width threshold disables the replicate warning
on exactly the array that needs it. Lowering the threshold makes it strictly worse. Sampling 2000
ids is decisive rather than probabilistic here, because the suffix is a whole-array property: every
EPICv2/MSA id carries one, no EPICv1/450K id does. Cost is then constant in the matrix width. (For
the record the full scan was affordable anyway -- 0.17s of `grepl` at 937k -- so this is about
keeping a per-call check free, not about rescuing an expensive one.)

**Standing refusal: `calc_clocks()` will not collapse replicate probes.** Do not re-propose it. (a)
It needs an external array manifest, which brings manifest versioning fragility and a tail of edge
cases (chr0 probes among them) into a package that otherwise ships a closed, self-contained
contract and reaches no network. (b) Doing it means allocating a second full copy of a very large
matrix inside the scoring call -- the same objection that sank `coerce_dnam()`'s
`as.matrix(data.frame)` in PR #3 sec 3.4. The collapse is the user's, done outside `calc_clocks()`,
where they can choose the manifest version and the aggregation rule. The package's job here is to
tell them it is needed, which is what the warning now does.

---

## 2026-08-03 -- One beta entry point, and no pre-flight surface

Written down because it was already true and nobody had said it: **`calc_clocks()` is the only
public surface that reads a beta matrix.** Verified over the exported surface -- only
`calc_clocks` and `predict_sex` take a `DNAm` argument, `predict_sex` reads it exclusively through
`calc_clocks`, `sim_DNAm` generates rather than reads, and every remaining export takes either
catalog arguments or an `mc_result`.

**Why state it now.** It settles a class of proposals in one line instead of re-arguing each one.
PR #3's `build_coverage_table()` (per-clock coverage on a bare matrix), the `report(DNAm)` arm, and
any future "dry run" or preflight helper are all the same shape: a **second beta reader**. The
objection is not duplication, it is decoupling -- a second reader is handed its own matrix, so its
verdict can be about a different object than the one that eventually gets scored, and nothing in
the package can detect the substitution. `predict_sex()` already shows the shape that is fine:
it is a `calc_clocks()` call whose output feeds a later call, which is composition, not a
pre-check.

**Nothing is bought by the second reader.** Scoring is a matmul over a matrix that is already
resident; the costly parts (materializing the betas, loading packs) are paid by any caller either
way. And both coverage gates are arguments, so `calc_clocks(min_clocks_coverage = 0,
min_samples_coverage = 0)` is already the full-report dry run, with the real numbers rather than a
prediction of them.

**Why the pre-flight model feels obligatory anyway, and why it does not transfer.** It is inherited
from upstream, where it is correct: an ENmix or minfi pipeline threads one object through its
steps, caches an expensive IDAT parse at the front, and has genuine cross-sample and cross-probe
stages where dropping a sample changes what follows. Within a coupled pipeline there is no
decoupling hazard -- the object *is* the state. Neither condition holds here. Users will still
arrive with the habit; the answer is that the model is right where they learned it and does not
apply to a calculator.

**The premise, kept visible.** This rests on scoring being cheap enough that running it is not a
commitment. A chunked or streaming path over an on-disk store weakens exactly that, which is the
one change that would put a preflight surface back on the table -- and it would need an explicit
decision at that point, not an inherited one. See `dev/to-do.md` sec 7, where the chunked front end
is still an open question.

**Not in scope:** internals under `calc_clocks()` obviously touch the matrix (`check_DNAm()`,
`col_stats()`, the scan and fill machinery). The rule is about the public surface and where a
matrix enters the package, not about who may hold a pointer to it.

---

## 2026-08-03 -- Working trees get native line endings; only the index is pinned to LF

`.gitattributes` said `* text=auto eol=lf` -- LF in the repo **and in every working tree, on every
platform**. The second half is now dropped: `* text=auto`, so the index is still always LF and the
working tree follows `core.eol` (CRLF on Windows).

**Why.** `Rcpp::compileAttributes()` -- which `devtools::document()` and `load_all()` both call --
generates `R/RcppExports.R` and `src/RcppExports.cpp` inside its compiled
`.Call("compileAttributes", ...)`, writing through a text-mode `std::ofstream`. That is CRLF on
Windows and there is no R-level knob for it. Under `eol=lf` those two files came back ` M` after
every document run, forever, on every Windows machine. `git diff` showed nothing (it normalizes
before comparing) while `git status` showed the modification, which is the confusing part: git's
index refresh short-circuits on a size mismatch and never reaches the content comparison, so the
CRLF file reads as changed even though its filtered hash matches the blob exactly.

**Why not just live with it.** It is cosmetic -- `git add` normalizes, so no CRLF can reach the
index -- but it trained every Windows collaborator to ignore a dirty `git status`, which is the one
signal that has to stay trustworthy. Cost of the fix is one attributes line; cost of the noise is
paid on every commit.

**What stays pinned to LF**, because these are read on a platform other than the one that checked
them out: `src/Makevars` and `src/Makevars.win` (they ship inside the source tarball, and a stray CR
lands in a GNU make variable value), `.Rbuildignore` (`readLines()` does not strip a CR off it
outside Windows, so the ignore patterns silently stop matching), plus `*.sh` and `*.py` on the same
reasoning. Everything else is native.

**Do not "fix" a Windows worktree by re-normalizing the whole tree.** Existing LF files stay clean
because their recorded stat still matches; only a file some tool *rewrites* can churn, which today
is exactly the two generated ones. The one-time repair is `rm` + `git checkout --` on those two --
a plain `git checkout` alone is a no-op, since git skips writing an entry it already considers
up to date.

---

## 2026-08-02 -- `score_associations()` ships as a disposable advisory, and it is the one sanctioned `cor()`

**This is a deliberate, scoped carve-out of "Correlation is never a numeric gate."** The invariant
stands, and its subject is unchanged: an **agreement gate** -- a test that asks "do we match the
oracle" -- may never be a correlation, because `cor()` is offset- and scale-invariant and cannot
tell "correct" from "uniformly wrong". Nothing about parity, unit tests or any numeric bound moves.

`score_associations()` is a different use. The correlation is not an agreement statistic against a
reference implementation; it **is the estimand** -- the cohort's score-age association, which is the
quantity the reference table reports. It gates nothing: it returns a frame, stops no scoring, filters
no clock, and emits no verdict. Read the invariant as binding on gates, not on reported statistics.

**The reference intervals are wide, and that is measured, not suspected.** 31 of the 69 `age_r`
prediction intervals include zero, median width 0.71 on a 2.0 scale, so for ~45% of clocks "no age
association" falls inside the expected range. This was quantified before the function was written.
Shipping a v1 on that basis is a deliberate call: the flags are a starting point for a user who wants
one, not a measurement anyone should rely on.

**The function reports and says nothing else.** It emits no message and returns only the frame. A
caveat on every call is noise in programmatic use, and the appropriate place for a limitation of this
kind is documentation. The cost is accepted: a caller who reads neither this entry nor
`dev/PR3-respond.md` gets two booleans without context. There is still no PASS/WARN/FAIL and there
must not be one -- that was the specific thing declined in PR #3 (see `dev/PR3-respond.md` sec 3.6).

**The reference table is taken as given, and we differ from the author on its construction.**
`inst/extdata/clock_reference.csv` and `data-raw/build_clock_reference.R` are shipped unmodified,
authored by `dsborrus`, whose meta-analytic work is the only artifact in PR #3 that could not be
sourced elsewhere -- the prediction intervals use the Higgins-Thompson-Spiegelhalter form, pooling
is DerSimonian-Laird, correlations are pooled on the Fisher-z scale, and `MIN_N = 20` is a
reasonable floor. The disagreement is about what the two-stage design can support downstream, not
about the execution. Stage one fits a separate `lm` per (dataset, clock) pair; stage two pools those
estimates. A pool of that shape summarizes datasets rather than samples, so the individual level is
not recoverable from it, and the only comparison available at runtime is a user-cohort statistic
against a reference-cohort statistic -- which inherits the user's study design, and is why the
intervals have to be as wide as they are to stay honest.

Our position is that this wants a one-stage hierarchical model instead (`score ~ age + sex + tissue +
(1 | dataset)` fit on pooled individual-level rows), shipping fixed effects + covariance + tau^2 +
residual sigma so a `predict` step gives a closed-form per-sample interval and covariates can be
marginalized when absent. That version has a known null -- 5% of samples outside a 95% interval --
so no threshold has to be chosen anywhere. It needs the per-sample tables behind
`HigginsChenLab/TranslAGE-workflows` and is a separate project, which is why the v1 here ships on
the aggregate table rather than waiting.

**Disposability is the design constraint, so it is written down as a contract.** The whole feature is
`R/score_associations.R`, `inst/extdata/clock_reference.csv`, `data-raw/build_clock_reference.R`, and
one `export(score_associations)`. Delete those four things and nothing else in the package changes.
It reaches the record through exactly one internal, `finalized()` -- the same re-finalize hook
`as.data.frame()` and `calc_accel()` use, so a multi-batch record with pending cross-sample columns
is not correlated stale -- plus the documented `$scores` / `$pheno` fields. It adds no dependency,
no S3 method, no catalog field, and no coupling to routing, coverage or the pack machinery. Keep it
that way: when the hierarchical version lands, this should be removable in one commit.

## 2026-08-02 -- Wang scores as a matmul plus two scalars, and the scalars stay runtime

**`score_DNAmSex_Wang()` no longer materializes the z-scored or centred matrix.**
`sum_j ((x_ij - m_i)/s_i - c_j) * r_j` is algebraically `(x_i . r - m_i * sum(r)) / s_i -
sum(c * r)`, so the projection is one `n x p` matmul against a vector plus two reductions over
`p`. The two `n x p` temporaries -- the z-score and the `rep(center, each = n)` broadcast -- and
the `rep()` itself all vanish. Same rearrangement `score_Zhang2019()` already used, so this is the
established form here, not a new trick.

**Measured** (ChrX, p = 4047, real tensors): n=500 30x, n=2000 33x (0.99s -> 0.03s over 10 reps),
and 247 MB of transient allocation per call gone at n=2000. The win grows with n because the old
form's cost is `O(n*p)` elementwise work while the new one's is the matmul BLAS was going to do
anyway.

**It reorders the summation, so it is not bit-identical.** Old vs new on the shipped tensors, over
100%/90%/50% coverage and both members: worst 1.1e-12 absolute on scores of scale ~230, worst
1.5e-11 relative (that relative is inflated by samples whose score is near zero, where a 1e-12
absolute miss is large in ratio). Against the 1e-10 abs / 1e-10 rel `core` parity gates that is
~100x and ~7x of margin respectively. **Re-run parity before trusting this**: EPICv1 previously
passed at 5.5e-12 abs (ChrX), so expect ~6.5e-12; the relative axis is the tighter one and is the
one to read.

**The two scalars are computed at scoring time, from `present`, and must stay there.** The
tempting next step is precomputing `sum(center * rotation)` upstream as an intercept -- it looks
like a catalog constant. It is not: both Wang members declare `imputation: omit`, so an absent CpG
is *dropped*, and the correct constant is over `present`, not the declared panel. Measured drift of
`sum(center*rotation)` on ChrX as coverage falls: 0.37 at 1% absent, 1.9 at 5%, 3.8 at 10%, against
a score scale of ~+/-137 -- i.e. 0.3% to 2.8%, applied as a **uniform offset to every sample**. That
is orders of magnitude past the parity gate and is exactly the failure a correlation check cannot
see, which is why the gate is a bounded per-element difference. `sum(r)` in the z-score half is
coverage-dependent for the same reason and drifts faster (0.55 at 1% absent).

**So the contrast with the SystemsAge sync precompute is the point:** that one survived review
because its arithmetic does not depend on which CpGs the user measured. Wang looks like the same
shape and is not. A precomputed constant here would be correct only at exactly 100% coverage and
silently wrong everywhere else.

---

## 2026-08-02 -- Wang parity: EPICv1 passes, cohort_450K is a declared gap

**Parity was run** (maintainer-authorized). **FAIL 0 / SKIP 32 / PASS 699**, 61s. Re-measured
standing state below; the figures in `CLAUDE.md` were stale and are updated. `R CMD check` not run.

**The first run failed 8 targets, all Wang, and the cause was ours.** `run_parity_target()` carried
its own projection -- `panels_union(clock_panels(seq_ids, packs))` -- so it fetched the 4047 / 284
panels and none of the 442533-probe ref. Every EPICv1 sample scored `NA`. This is the **same defect
as the `clock_cpgs()` one two entries below, in a second copy of the same logic**, which is the
argument for the fix that landed: `sequence_cpgs()` (`R/clock_cpgs.R`) is now the one union, called
by both, with an always-on test in `test-score-wang.R` asserting they agree. A private copy of "what
must I measure" is the bug; deleting the copy is the fix.

**EPICv1 then passed with room:** ChrX 5.5e-12 abs / 2.3e-13 rel, ChrY 9.4e-13 / 9.2e-14, against
`core`'s 1e-10 on both axes. First evidence the branch is right against something other than a
re-implementation of it.

**Free finding, since the cohort carries `Female`:** the karyotype call agrees with reported sex
**71/71** on EPICv1 (43 Female, 28 Male, no aneuploid calls). Recorded as an observation, **not
wired as a test** -- see the reasoning in `dev/moment-domains-plan.md` sec 12: parity gates "we match
the oracle", and per-sample agreement with a binary `Female` cannot represent `47,XXY` / `45,XO`, so
any true aneuploid would be a permanent red. What it does settle cheaply is inversion: a swapped X/Y
binding scores 0/71, not 71/71.

**cohort_450K is a genuine gap and is now declared** -- the first two entries in the deliberately
empty `KNOWN_PARITY_GAPS`. The deposited 450K matrix carries **no sex-chromosome probes** (473034 of
485512), so both panels are 0% present; upstream's own fixture declares all 4047 / 284 missing and
expects `0`, which is what the author's code returns from summing an empty panel. We refuse instead:
`check_coverage()` treats `ratio == 0` under a non-`vendor_mean` policy as an unconditional stop,
independent of `min_clocks_coverage`. **That rule is right and was not relaxed to make a fixture
pass** -- a `0` here is not a small score, it is the `Female` quadrant of the sign map, so scoring it
would be a meaningless number wearing a real answer's clothes. Upstream may want to drop those two
fixtures; they can only ever assert `0 == 0` from a run that scored nothing.

**Re-measured standing state**, replacing the pre-Zhang-split figures:

| block   | targets | note                                  |
|---------|---------|---------------------------------------|
| core    | 146     | includes Zhang2019BLUP and Wang@EPICv1 |
| fitage  | 28      |                                       |
| packs   | 56      |                                       |
| horvath | 30      | all skipped (DECISIONS 2026-07-25)    |

260 targets + 2 PhysAge + 1 census = **263 blocks**. 228 targets run x 3 expectations, + PhysAge
2 x 6, + census 3 = **699**. Skips are 30 horvath + 2 Wang@450K = **32**.

`core` grew from the recorded 130 to 146. Six of those are accounted for (Zhang2019BLUP x 2 cohorts,
Wang x 2 members x 2 cohorts); the other ten arrived with the uncommitted `R/sysdata.rda`
regeneration and are **not** explained here -- stated as a measurement, not a claim.

---

## 2026-08-02 -- `predict_sex()` reads the karyotype call; it does not reimplement it

Always-on suite **1179 pass / 0 fail / 264 skip / 0 warn**. Parity not run; `R CMD check` not run.
Implies one new export: `export(predict_sex)`. `document()` also re-sorted `export(calc_accel)`
into alphabetical position -- incidental, not a surface change.

`R/predict_sex.R`: `calc_clocks()` on both `DNAmSex_Wang` members, `as.data.frame(long = FALSE)`,
then `apply_karyotype()`. The whole call is `default` + `rules` read off
`group_entry("DNAmSex_Wang")[["routing"]][["karyotype_call"]]` -- **the first consumer of
`mc_groups` in `R/`**, which is why `group_entry()` (`R/accessors.R`) exists at all: group
declarations are read through an accessor like clock ones, never off `mc_groups` in place.
Verified against all four quadrants: `X<0,Y>0 -> Male`, `X>0,Y>0 -> 47,XXY`,
`X<0,Y<0 -> 45,XO`, and `X>0,Y<0 -> Female` **because no rule covers it**, which is also what an
exact 0 on both axes gives. That last pair is the whole reason this is read rather than written: a
symmetric four-quadrant table agrees on the interior and diverges on the boundary, and the boundary
is every sample in `cohort_450K`.

**Two deliberate divergences from the author.**

1. **A sample missing either score gets `NA`, not the default label.** `wateRmelon::estimateSex`
   overwrites with `predicted_sex[which(...)]`, and `which()` drops `NA`, so an unscorable sample
   silently keeps `Female`. We can produce that state -- a matrix carrying the panels but not the
   z-score ref -- and reporting a sex for a sample nothing was measured for is the same error as
   reporting coverage for a sample it is not true of.
2. **The operand -> input binding is checked against the clock id.** The catalog pairs rule keys
   (`chrX`, `chrY`) with `inputs` by declared order only; the actual binding lives in the prose
   `by_chromosome` field, which is not machine-readable. Positional pairing alone means an upstream
   reorder inverts every call with nothing to catch it -- fixture goldens store scores, not calls.
   So `karyotype_inputs()` also asserts `endsWith(tolower(id), tolower(key))` and stops otherwise.
   This is **not** the banned accessor search: the payload is already declared in `inputs`, and the
   string test only cross-checks two declarations against each other.

**The rule engine is unit-tested, against the usual altitude rule, on purpose.** No fixture covers
it in either repo (parity goldens are scores) and `random_betas()` cannot be steered into a
quadrant, so an output-level test can reach the shape and the label set but never the mapping.
`apply_karyotype()` therefore takes `key -> score vector` and is tested directly on synthetic
quadrants including both boundary cases. If a `predicted_sex` fixture ever lands upstream, that
becomes the golden and this drops to a smoke.

---

## 2026-08-02 -- `clock_cpgs()` reports declared moment refs; panels are unchanged

Always-on suite **1173 pass / 0 fail / 264 skip / 0 warn**. Parity not run; `R CMD check` not run.

Settles the question the entry below left open. `clock_cpgs()` is now
`panels_union(...)` unioned with `unlist(resolve_moment_domains(sequence))`.

**Two placements got conflated, and only one of them is dangerous.** Making the ref a **panel role**
inside `clock_panels()` is wrong: `needed_union` is `panels_union(panels)` over both roles, so it
would widen Wang's coverage denominators from 4047 / 284 to 442533 and break "a moment domain is not
a panel". But `clock_cpgs()` is a **leaf** -- `clock_panels_union()` has exactly one caller, and
`mc_spec()` builds `needed_union` from its own `clock_panels()` call -- so adding the ref there
cannot reach coverage at all. The first was briefly taken as an argument against the second, and it
is not one.

**What decided it: `clock_cpgs()` is on a shipped path, not just in front of the simulator.**
`dev/id-streaming-plan.md` sec 8 makes "project with `clock_cpgs()`, block, score, bind" the
supported route for a cohort that does not fit. A user following that advice for `DNAmSex_Wang` would
project away the ref and get an all-`NA` column with no error -- so under-reporting here is a bug on
a documented workflow. The name is also a promise: this function answers "what must I measure", and
a CpG whose absence turns the score into `NA` is measured input by any reading.

**The rule is "is the domain declarable", not a clock.** `resolve_moment_domains()` already carries
`NULL` for the whole-matrix domain, so `Zhang2019`'s `"full"` drops out of `unlist()` on its own --
there is no set to name -- while Wang's ref contributes its 442533. Same distinction the kernel's
`moment_sets` keys on, reused rather than restated: no clock is named, and a future ref-declaring
clock is served without an edit. `Zhang2019EN` is therefore still exactly its 514-CpG panel, and
what a full-panel clock needs beyond that is still carried by `note_full_panel_clocks()` at
`calc_clocks()` time. Both are tested.

**`sim_DNAm()` stays exactly `clock_cpgs()`.** An earlier pass put this union in `sim_DNAm()`
instead; that left the public answer wrong, split one rule across two functions, and would have made
the simulator the only place that knew the real requirement. The blank/NA screen moved up with it,
out of `clock_panels_union()` and onto the whole answer, so it is applied once.

**Costs, accepted.** `clock_cpgs("DNAmSex_Wang_ChrX")` now answers 446580 rather than 4047, and a
Wang sim is that wide (~14 MB at n = 4, transient). `remove =` is diluted for this family -- 5
dropped columns out of 446580 hit a panel probe ~1% of the time, so it stops being a way to force
absent-panel behaviour there. Nothing uses `remove` today. And the smoke tier's 2 expected warnings
are gone, since it now scores Wang for real -- which leaves the `NA`-with-no-reference path covered
**only** by its named test in `test-score-wang.R`. Better place for it, but load-bearing now: do not
delete it as redundant.

---

## 2026-08-02 -- the Wang branch reads its own domain, and says when it cannot

Always-on suite **1167 pass / 0 fail / 264 skip / 2 warn**. The 4 red tests from the two entries
below are green. Parity not run; `R CMD check` not run. The skip count jumping 2 -> 264 is not a
regression: `test-fixtures-parity.R` used to die at **file** level on `score_type()` (one of the 4
failures), so it emitted no targets at all; now it loads and its parity tier skips normally.

`R/score_DNAmSex_Wang.R` + a `(DNAmSex_Wang, wrapper) -> "DNAmSex_Wang"` group hook. The branch is
the declared recipe and nothing else: `sample_scale` against the `DNAmSex_Wang:zscore_ref` domain,
`center_scale` by the declared `center`, `project` onto `rotation`. Verified against an independent
R implementation on both members -- max abs difference 1.4e-13 on scores of magnitude ~40, i.e. 3e-15
relative.

**The operand names are read off the recipe, not hardcoded.** `recipe_step_op(id, op)` finds the step and
the step names its component (`center`, `rotation`); `SystemsAge` hardcodes its component names and
that was the tempting precedent. Reading the declaration costs six lines and means a rename upstream
is a `catalog_bug()` naming the component instead of a silently-NA score vector. Same reason the
branch stops when `center_scale` declares a `scale` it does not apply, and when `center` and
`rotation` do not cover the same CpGs: `center[present]` on a missing name returns `NA`, which would
propagate through the projection and produce an all-`NA` column with no diagnostic.

**A sample with fewer than 2 observed reference CpGs is scored `NA`, noted, and warned about.**
This is the `n < 2` guard in `split_moments()` arriving at a consumer. It has to be said out loud
because **coverage cannot see it**: the ref is not a panel, so a sample can hold 4047/4047 of the
scoring panel, report 100% coverage, and still be unscorable. The mechanism already existed for
exactly this -- `note_scoring_failure()` + a plain `warning()`, the shape `score_normalized()` uses
for a failed BMIQ fit -- so the samples land in `$provenance$scoring_failures` rather than only in a
message the user may have suppressed.

**Consequence, and it is load-bearing for how the smoke tier reads:** `sim_DNAm()` materializes
declared panels, and a moment ref is not one, so the smoke tier hands both Wang members a matrix
their ref meets nowhere and gets `NA` plus that warning back. The tier still does its job (default
configuration, no error), but it is **not** a numeric check on this family, and the 2 warnings in a
green run are expected rather than a smell. `test-score-wang.R` builds panels plus a slice of the
ref to get finite scores. Whether `clock_cpgs()` should report a clock's moment ref as required
input -- which would make `sim_DNAm("DNAmSex_Wang_ChrX")` scorable and fold the ref into the public
"what must I measure" answer -- is left open, not decided by silence: see
`dev/moment-domains-plan.md` sec 12 item 3.

---

## 2026-08-02 -- `DNAmSex_Wang` stays in the callable pool

Decision only; no code. Resolves the contradiction between `dev/moment-domains-plan.md` sec 9
("take Wang out of the pool", left open) and the settled `predict_sex()` shape (thin sugar over
`calc_clocks(DNAm, c("DNAmSex_Wang_ChrX", "DNAmSex_Wang_ChrY"))`, which requires both members to be
callable by name). Both could not hold. **The design wins: `resolve_clocks()` accepts both.**

**Two of the three arguments for removal had expired**, which is what made this lopsided rather than
a genuine toss-up. Removal was offered mainly to buy K = 1 in both front doors -- but the
multi-domain sweep is built and tested, so that simplification no longer exists to be bought. And
the kernel was the stated blocker for the whole family; it no longer is. What was a
cost-of-delay argument in favour of the smaller option is now an argument for nothing.

**The surviving argument against is real but misfiled.** A sex PC in `$scores` is an `n x k` double
a user can hand to `calc_accel()` and get confident nonsense from, since nothing in the type says
only the sign is meaningful. That is a fact about **output typing**, and it generalizes to every
score that is not an age -- `EpiTOC2`'s mitotic index has the same shape. Solving it by removing one
clock from the pool would fix one instance of a general problem while creating specific work:
Wang would never return an `mc_result`, so its coverage record needs a new home, and `predict_sex()`
would have to drive `mc_spec()` / `mc_cohort()` / `score_cohort()` directly instead of the one
public entry point. If the hazard is worth addressing it should be addressed as typing, for all such
clocks, not as pool membership for one.

**Consequence to expect:** `calc_clocks(DNAm, "all")` will span `"full"` and
`DNAmSex_Wang:zscore_ref`, so **K = 2 becomes reachable through the front door** for the first time
-- the case the moment-domain work was built for, and one no test can currently reach.

---

## 2026-08-02 -- moments are keyed domains, not one banked pair

Always-on suite 1154 pass / 4 fail / 2 skip. The 4 are the same parked `DNAmSex_Wang` routing gap
as the entry below (`score_type()` stops on both members, so smoke plus the two census tests that
sweep every catalog clock fail); none are new. Parity not run; `R CMD check` not run.

**`spec[["needs_moments"]]` (logical) becomes `spec[["moment_domains"]]` (key -> cpgs).** The
logical could only ever say "bank the whole-matrix moments", which is why the entry below had to
record Wang banking nothing at all as a deliberate consequence. A request spanning both kinds --
`calc_clocks(DNAm, "all")` the moment Wang is callable -- had no representation, and the failure
mode was silent: Wang would have taken Zhang's whole-matrix moments. Now `mc_cohort()` banks
`key -> list(mean, sd)` and a branch reaches its moments by clock id (`block_domain_moments(block,
id)`, which derives the key), so reaching for an unbanked domain is a `stop()` rather than a wrong
number and **no branch ever spells a domain key**.

**The key is derived, never assigned:** `"full"` (`FULL_MOMENT_KEY`) for a ref-less step, else
`<group_id>:<ref_name>`, minted in `clock_moment_key()` (`R/accessors.R`). `clock_moment_key()` is
split out of `clock_moment_domain()` so a caller wanting only the key -- the `resolve_moment_domains()`
dedup, `clock_needs_full_panel()`, every score branch -- never resolves the ref's CpGs (442533 of
them for Wang) just to throw them away. Minting it in the accessor rather than in
`mc_spec()` is what makes two clocks sharing a ref collapse to one domain by construction -- both
Wang members land on `DNAmSex_Wang:zscore_ref` without anything comparing CpG sets. The `:` is what
keeps a declared domain from colliding with `"full"`, and a group id cannot contain one.

**One pass, whatever K is.** `col_stats(obj, cols, moment_sets)` labels each column with the bitmask
of the sets containing it, accumulates into one Welford atom per distinct mask, and merges atoms per
set with Chan's formula. A column in several domains is read once. Welford rather than the additively
obvious `(n, sum, sum_sq)`: on betas in [0,1] the cancellation is about one digit in sixteen, which
would be unremarkable except that Wang's output is a **sign**, so the error concentrates exactly on
the samples near a quadrant boundary, where a lost digit flips a call instead of nudging a number.
Atoms are indexed by observed mask, not `2^k` -- fine at k = 3, not at k = 8.

**The rejected cheaper option is one sweep per domain**, which needs no kernel change and about
30 lines of R. It re-reads every overlapping column once per domain (~+50% scan on a mixed request).
That was the fallback if the mask proved fiddly; it did not.

**Element validation lives in R (`check_moment_sets()`), not in the kernel.** Not a style call: on
this toolchain `as<IntegerVector>()` on a `NULL` or character element is **not a catchable
condition** -- it terminates the process (exit 127, `tryCatch` bypassed), and the same throw is
catchable under plain `sourceCpp`, so the cause is build configuration (`~/.R/Makevars.win` carries
`-UNDEBUG -g -O0`; the package adds `-fopenmp` with `-static-libgcc`). Validating in R removes the
only path that reaches it. Worth knowing independently: `expect_error()` on any future kernel path
relying on a C++ throw is unusable locally.

**The `n < 2` guard is load-bearing, not defensive.** The kernel reports an unobserved row as
`n = 0, mean = 0, m2 = 0` -- counts disambiguate rather than the kernel emitting NaN -- so R applies
two different thresholds in `split_moments()`: a mean needs `n >= 1`, an sd needs `n >= 2`. `n = 1`
yields `NaN` naturally (`0/0`); `n = 0` yields `sqrt(0 / -1) = -0`, which reads as real zero spread
and divides to `Inf`. Nothing upstream makes this unreachable: the dead-sample gate keys on the
**scoring** panels, and for Wang the scoring panel and the moment ref are **disjoint** (0 shared
probes, 4047 and 284 sex-chromosome probes against a 442533 autosomal ref), so a sample can clear
every existing gate with full sex-chromosome coverage and still have zero ref observations. An
empty domain is likewise a data fact reported as `NA`, not a usage error.

**Coverage is untouched, deliberately.** A moment domain is not a panel: the sets index `DNAm`
directly, so a ref never widens `needed_union` and the declared `n_cpgs` do not move.

**The 1 GiB accumulator guard was written and then removed.** It refused `NS * nr * 20` bytes above
a ceiling, on the theory that `NS` is emergent -- it counts distinct membership *patterns*, not sets
-- so a set system the caller thinks is small could allocate more than it looks like. The bound that
kills the guard is in the same sentence: a signature is a **per-column** fact, so `NS <= ncol`
always, and the accumulators can never exceed `20 * nr * nc` -- **2.5x an `obj` R has already
materialized**. The check could not fire before R itself had failed to allocate, so it was buying a
different error message, not a different outcome. The only shape approaching the ratio is a matrix
narrow enough that nearly every column has a unique pattern (measured: ~430 MB of input at
`nc = 3`), which is unreachable through `scan_missing_cpgs()` -- `col_stats()` is internal, and its
sets come from the catalog, where `K <= 2` gives `NS <= 3`. What survives is the *fact*, as a
three-line comment at the allocation; the ceiling constant, the branch and the four-argument `stop()`
are gone.

**The mask-width check stays, and is not the same kind of thing.** `1u << 8` does not fit a
`uint8_t`, so that bound is correctness. Do not read the removal above as licence to drop it: one
guard was a resource heuristic against an unreachable state, the other is the type's own limit.

**The mask is a `uint8_t`, and K is dynamic.** The width went `uint64_t` -> `uint8_t`: `mask[]` is
one entry per column, so at EPICv2 width that is 7.5 MB -> 0.94 MB, and -- the larger win -- an
8-bit signature has only 256 values, so the signature -> slot map is a flat 256-entry table instead
of an `unordered_map`, and `<unordered_map>` leaves the file. This does **not** undo "index atoms by
observed mask, not by `2^k`": that argument is about the `nr`-length accumulator **blocks**, which
still get compact ids. Only the lookup became a table, because 1 KB is free at any width a byte
holds.

**What K is, since this was gotten wrong once.** K is the number of distinct domains **the requested
sequence** needs: 0 for a clock with no `sample_scale`, **1 for a Zhang arm on its own** (the common
case today), 2 for a run spanning Zhang and Wang. The catalog-wide count of 2 is a *ceiling* over
every clock that exists, not a per-call value. An intermediate version hardcoded `K = 2` and derived
a closed-form slot layout from it; that pinned the kernel to the one value the front door never
produces, so `calc_clocks(DNAm, "Zhang2019EN")` -- the most ordinary call a `sample_scale` clock has
-- failed with "exactly 2 moment sets are required", and the suite went from 4 failures to 18. Only
the *width* is fixed. Do not hardcode K, and do not read a census of the catalog as a contract on
one request.

**The width is enforced where a maintainer sees it.** Dropping 64 -> 8 moves the ceiling closer to a
plausible upstream future, and the refusal would otherwise fire at runtime inside a user's
`calc_clocks()`. So `MAX_MOMENT_SETS` (`R/missingness.R`) mirrors the kernel, one test asserts the
two agree on the boundary, and an always-on census asserts the shipped catalog's distinct domain
count fits. A ninth domain is then a red suite for whoever synced it, not an error for whoever
scored with it -- which is what makes the tighter type safe rather than merely smaller.

**K is still 1 through the front door.** Wang has no scoring branch, so `mc_spec()` only ever mints
`"full"` today; the K > 1 path is covered by unit tests on `col_stats()` and `scan_missing_cpgs()`,
not by `calc_clocks()`. That is the point of doing the seam first -- the Wang branch is now a branch,
not a re-plumbing. Whether Wang stays in the callable pool at all is still open and independent
(`dev/moment-domains-plan.md` sec 9).

---

## 2026-08-02 -- `shared[]` is a declaration, not build scrap; full-panel splits on `ref`

Always-on suite 1074 pass / 4 fail / 2 skip. The 4 are the parked `DNAmSex_Wang` routing gap
(`score_type()` stops on both members, so smoke + the two census tests that sweep every catalog
clock fail); none are new. Parity not run; `R CMD check` not run.

**`shared` leaves `CATALOG_BUILD_ONLY_FIELDS`.** The handoff in `dev/update-DNAmSex.md` asked
upstream to re-declare `zscore_ref` as a `probe_sets[]` entry, on the premise that the trim left
the tensor with no declared pointer. Upstream declined (`dev/reply-DNAmSex.md`) and was right on
the facts, which were checkable here: `build_group_bundles()` runs *before*
`trim_build_only_fields()` (`sync.R:2184` vs `:2187`) and `mc_bundles` / `mc_groups` are assigned
untrimmed (`:2192-93`), so the payload was always in the bundle and the path was always in
`mc_groups[[gid]][["shared_tensors"]]`.

What the trim actually destroyed was narrower and worse-shaped: the **`name` -> `file` binding**.
`shared_tensors` is paths only, so nothing could resolve `recipe[["Xz"]][["ref"]] == "zscore_ref"`
to a file. That is a resolution gap, and the accessor invariant forbids closing it by searching the
bundle. Retaining `shared` closes it exactly, and `SHARED_FIELDS` is already `c("name", "file")`,
so the entry carries the binding and nothing else -- 33 clocks, two strings each. Measured: the
regenerated `sysdata.rda` did not grow.

**Why not accept the `probe_sets[]` copy anyway, since it was additive and cheap.** It would be a
second owner of a fact the recipe operand already owns, and it would be *inert* -- `sample_scale`
resolves through the shared namespace, so nothing would ever compare the two. Drift would be
silent on both sides. Upstream's ownership rule and our "accessors read declarations" rule are the
same rule seen from two ends; the fix belonged wherever the duplication would not be created, and
that was here.

**`shared` is keyed by `name` in `key_catalog_lists()`**, joining `components` / `probe_sets` /
`recipe`. A resolver over an unkeyed list is a scan, and the point of the change was to make the
lookup a lookup.

**`clock_needs_full_panel()` now means "a `sample_scale` step with no `ref`", not "has a
`sample_scale` step".** It was user-visible and wrong: it told a caller that Wang "scores against
every column of `DNAm`" when Wang's moments come from a declared, closed 442533-probe set. The
distinction is not size -- Zhang's moment set is *whatever the caller supplied*, so it has no
closed membership at any time, including sync time. That is a difference in kind, and `ref`
presence is the seam that expresses it.

Note the knock-on, which is deliberate: `spec[["needs_moments"]]` is `length(full_panel) > 0`, so
Wang now banks no moments at all. That is correct while Wang has no scoring branch, and it means
whoever writes that branch cannot get whole-matrix moments by accident -- they have to do the
`col_stats()` work first.

**Still parked, and the block moved.** `col_stats(row_moments = TRUE)` sweeps the subset *plus*
its complement by construction (`col_stats.cpp:180`), so there is no way to obtain moments over a
declared subset today. That blocks the Wang scoring branch, hence the `score_type()` group hook,
hence `predict_sex()`. The upstream contract is no longer the blocker; the kernel is.

**Do not hard-code the karyotype quadrants.** Upstream found our prose wrong while structuring it:
`wateRmelon::estimateSex` assigns `Female` unconditionally and overwrites with three rules, so
`X>0 & Y<0` is never evaluated and a four-quadrant table diverges wherever a score is exactly 0.
`cohort_450K` is that case -- every score 0, all 80 samples labelled Female, which no quadrant
matches. Both suites store scores, not calls, so neither would have caught it. The map now ships
as `mc_groups[["DNAmSex_Wang"]][["routing"]][["karyotype_call"]]` (a `default` plus three rules);
derive from it.

---

## 2026-08-01 -- the multi-batch test is `is_multi_batch()`, and a frame may decline to build

Always-on suite 1051 pass / 0 fail / 260 skip. Parity not run; `R CMD check` not run.

**The one multi-batch test moved out of `drop_single_batch()` into `is_multi_batch()`**, which
`drop_single_batch()` now calls. Nothing about the exit schema changed -- this is a rename plus one
new reader.

**Every exit frame now reads it to skip building the batch column instead of building and then
dropping it.** `shape_scores()`'s long frame is n x k rows and the batch column is one of its four,
so at a single batch -- every record that has not been through `rbind`, i.e. the common case -- a
quarter of the frame was allocated to be deleted three lines later. Measured at 1.2e6 rows: 37.1 MB
/ 11.98 ms -> 28.4 MB / 9.98 ms, and the two forms are `identical()`. Scales linearly, so ~38 MB at
5e6 rows.

**The invariant this looks like it touches is about the output schema, and the output is
unchanged.** What the 2026-07-31 entry protects is the four exit frames agreeing on whether the
join key exists; a frame that never builds a doomed column and one that drops it are the same
frame. The real risk in a conditional build is the *predicate* getting a second home and drifting
-- which is why the test is named and stated once rather than inlined at the build site. That was
the whole cost of the change, and it is ~10 lines.

**`drop_single_batch()` still runs at all four exits, including the one that no longer needs it.**
Assigning `NULL` to an absent column is a no-op, so the call is free, and keeping it means the
output is correct even if a build condition is later changed or removed. Do not "clean up" the
now-redundant call: it is the gate, and the conditional build is an optimization underneath it.

**All four exits decline to build, not just the one with the measurable win.** The saving is
overwhelmingly in `shape_scores()`'s long branch; the coverage frames were converted for
uniformity, and their numbers are small to nil. `samples_coverage()` at 36000 rows: 22.3 MB ->
19.9 MB (-11%), time within noise -- its frame is six columns of mixed type and `rbind` dominates,
so the batch column is a smaller share than in the four-column score frame. `clocks_coverage()` is
(clocks x batches) rows, so a full catalog at one batch saves under a kilobyte. Neither was worth
doing on its own; both were worth doing so that one rule holds at every exit and no reader has to
work out why one frame builds a column it will lose.

**The two coverage frames thread a flag; they cannot add the column afterwards.**
`samples_coverage()` builds parts per batch, `rbind`s them and then NA-filters, so after the filter
there is no per-row batch left to reconstruct. The flag therefore goes down through
`clock_sample_rows()` and `panel_rows()`, and `empty_sample_rows()` takes it too because it seeds
the `rbind` and a mismatched seed is a hard error. `b` still masks the rows in the loop -- only the
label is withheld -- so the masking never depends on the flag. `empty_sample_rows()` takes a
**logical, not a label**: at zero rows there is no value to carry, and reusing the loop's `b` there
would read the last batch, or fail outright on a record with none.

**`clocks_coverage()` keys on the provenance vector, not on `names(per_clock)`.** The count it
iterates and the count that decides the column are different questions, and CLAUDE.md already pins
the answer to the second one. Passing `if (keep) b else NULL` keeps the loop reading `per_clock`'s
names while the decision reads `provenance[[mc_batch_id]]`.

**Inverting to an `add_batch_column()` that only ever adds was considered and rejected** -- it is
the cleaner one-function form, but `samples_coverage()` cannot use it for the reason above, so it
would need a second path for that site and stop being one function.

---

## 2026-07-31 -- `accel_id` names the spec, and one spec per call

Always-on suite 1051 pass / 0 fail / 260 skip. Parity not run; `R CMD check` not run.

**`calc_accel()` long output carries `accel_id` beside `clock_id`**, derived from the call's
`(formula, type)` pair -- `Age_accel`, `Age_Female_accel`, `Female_diff`, and bare `diff` when
there is no rhs. `long = FALSE` is that same frame pivoted: one column per pair, named
`<clock_id>_<accel_id>`. Wide accel column names therefore changed; they were bare clock ids.

**One `(formula, type)` per call. There is no grid argument.** Every per-fit diagnostic is scoped
to one design matrix -- `say_fill_batch()`, the "too few complete samples" warning that names the
dead clocks, `merge_accel_data()`'s conflict error -- so a grid either fires each of them N times
or merges them into something naming neither call. A grid also forces a cross-product-vs-paired
semantics decision, which is an interpreter. Composition covers it instead, and `accel_id` is what
makes the composition safe: `rbind` over long frames from several calls is unambiguous because
`(id, clock_id, accel_id)` stays unique, and wide frames `cbind` because the label is in the names.

**Term labels go in verbatim -- not lowercased, not slugged.** `~ Age + Female` gives
`Age_Female_accel`, because `Age` is the column the caller has in their pheno and the label should
round-trip to it. Lowercasing only makes sense if `calc_clocks()` case-folds pheno names on the way
in, which it does not and which is a much larger change. A non-syntactic term (`I(Age^2)`) lands in
a wide column name as-is; that is the same choice `shape_scores()` already makes for clock ids via
`as.data.frame(optional = TRUE)`.

**The first draft slugged (`tolower` + non-alphanumerics to `_`) and added an `accel_id =`
override to escape the collisions slugging created** -- `Age.x` and `Age_x` both becoming `age_x`.
Verbatim labels do not have that problem, so the override went with it: derived, never assigned,
same as the batch label. The one residual case is a `_` inside a pheno column name colliding with
the `_` join (`~ Age_Female` against `~ Age + Female`), which is a coincidence in the caller's own
data and does not buy an argument. Note the `is_auto_label()` reasoning (2026-07-30) is about
*gating* on a label; nothing gates on `accel_id`, so an override would have been cheap -- it was
dropped for having no remaining job, not because it was unsafe.

**`accel_id` is not `resid_type`.** `type` already names the `type =` argument, and bare
`type = "diff"` with no formula residualizes nothing, so "resid type" is false for exactly the case
it most needs to label. `accel_id` parallels `clock_id`, and both are join keys.

`MC_ACCEL` lives in `R/calc_accel.R`, not `R/constants.R` -- one file uses it, which is that
file's stated rule. It is deliberately **not** reserved against `data =` the way `MC_BATCH` is: it
never enters the pheno or the formula namespace, so a caller's `accel_id` column cannot collide
with it.

---

## 2026-07-31 -- the batch label is multi-batch only at the exits

Always-on suite 1044 pass / 0 fail / 260 skip. Parity not run; `R CMD check` not run.

**`mc_batch_id` now appears in an exit frame only when the record spans more than one batch.**
`as.data.frame()`, `calc_accel()`, `clocks_coverage()` and `samples_coverage()` all drop it
through one helper, `drop_single_batch()` (`R/mc_result.R`). Nothing internal changed: `$provenance`
still always carries the per-sample vector, `calc_accel()` still always offers `mc_batch_id` to
the formula and still always errors when `data =` supplies it.

**This is not a new policy -- `print.mc_result` has done it since it was written.** The printer
emits its `$provenance [N batch(es)]` block only when `N > 1`, under the comment "a single-pass
record has nothing new to say here". The exit frames were the outlier, not the change. Anyone
re-opening this should start from that: the record already declined to show a single-batch label
in the one place a user always looks.

**The first argument against was wrong on the facts and is recorded so it is not re-run.** It ran:
a user who scores two batches separately, exports each with `as.data.frame()`, and stacks the
frames in dplyr a week later gets two single-batch records, so dropping the label removes it from
exactly the workflow it exists for. That does not hold. The label is a bare 16-hex `xxhash64` of
the id column -- `2f6267c9ce1fa180` -- and a user hand-stacking two frames labels them with
something they chose, not with our hash. At one batch the column is a single repeated opaque value:
its information content is zero, not merely low. The join argument fails the same way; with one
batch `clock_id` is already unique, so `merge(clocks_coverage(x), samples_coverage(x))` still
resolves.

**What decided it was comprehension, not noise.** The noise problem was already solved once by
putting the column last (2026-07-31, "one name for the batch label"). The live problem is that
`mc_batch_id` looks alarming, most users never call `rbind`, and **prose docs are deferred**, so
there is currently nowhere for a user to look up what the column is. Shipping an unexplained hash
to the majority case to serve the minority one is the wrong default.

**The surviving cost is accepted knowingly: this is a data-dependent schema.** `df[["mc_batch_id"]]`
returns `NULL` rather than erroring, `dplyr::select(mc_batch_id)` errors, and either way the code
works on the author's record and fails on the user's. That is the price, and the mitigation is that
**every multi-batch-only decision reads one definition** -- `batch_labels()` (`R/mc_result.R`),
`unique(x[["provenance"]][[MC_BATCH]])`. The printer and all four exits go through it. If they
diverged, the two coverage frames could disagree about whether the join key exists, which is
strictly worse than always carrying it.

`per_clock`'s names were the other candidate and are what the printer used before this entry. They
are not the source: they answer "which batches have coverage records", while the column is filled
from the per-sample vector. The two agree except under a 64-bit label collision -- and there,
provenance gives the more honest answer, because two batches sharing a label cannot be told apart
by a column carrying that label anyway.

**Correction while here: the label is 16 hex (64-bit), not 12 hex (48-bit).** `batch_hash()` returns
`digest(algo = "xxhash64")` untruncated. The 2026-07-30 entry below and `CLAUDE.md` both said 12,
and the collision figure was computed from 48 bits (1.8e-11 at 100 batches); the real figure is
~2.7e-16. `CLAUDE.md` is corrected in place. Nothing downstream depended on the width -- the
conclusion was "no gate on labels" and it only gets stronger.

Not chosen: an `as.data.frame(x, batch = FALSE)` argument. Four functions gain an argument for a
cosmetic default, and the caller still reasons about a conditional schema -- just one they picked.

---

## 2026-07-31 -- R/constants.R holds the shared constants only

Housekeeping pass over the SCREAMING_CASE constants. Always-on suite 1036 pass / 0 fail / 260 skip.
Parity not run; `R CMD check` not run.

**A package namespace is already one flat, sealed environment**, so nothing here was ever defined
twice -- `loadNamespace()` locks the bindings and every `R/*.R` sees every other file's top-level
objects. The whole question was where a human looks, which is why the answer is a plain file of
`NAME <- value` and not a list or a locked env: those buy nothing the namespace lock already gives,
and under the `[[`-only rule they would read as `MC[["NORM_SCHEMES"]]` at every site.

**Only cross-file values moved.** `R/constants.R` holds `MC_BATCH` and the four `NORM_*` policy
sets, which are read from four files between them and whose *names* are the entire point --
a bare `"quantile"` in `coverage.R` says nothing, `NORM_SCHEMES_FILL` says why it is there.
Everything used in exactly one file stayed next to its caller (`MC_ASSET_SUFFIX`,
`PACK_SCORE_TYPES`, `WRITE_SIM_EXTS`, `MC_TAGS`, `MIAGE_*`): `ADULT_AGE <- 20` three lines above
its only reader documents better than the same line in a constants file, and moving it makes you
jump files to read one branch.

**`NORM_ROLES` and `STACK_NAMESPACES` stay in `R/accessors.R` on purpose**, and both now carry a
one-line comment saying so. `data-raw/sync.R` does `source("R/accessors.R", local = mc_runtime)`
and reads them straight back out; hoisting them turns `mc_runtime[["NORM_ROLES"]]` into `NULL`,
which is a sync failure with no error at the point of the mistake. This will look like an unfinished
hoist to the next reader -- it is not.

**Constants that were only naming a well-known number were deleted, not moved.** `ADULT_AGE`,
`LOG_AGE_OFFSET`, `WARN_COVERAGE_MARGIN` and `MC_DEFAULT_RELEASE_REPO` each had exactly one use
site within a few lines of their definition, and the name restated what the literal already said.
The Horvath anti-transform is `21 * exp(x) - 1` / `21 * x + 20` to everyone who reads that
function. `MIAGE_LOWER`/`MIAGE_UPPER` collapsed into one `MIAGE_BOUNDS` because the starts grid
derives from the bounds -- that coupling is real and worth a name, two separate scalars were not.

**Four `mc_batch_id` literals survive, and that is not an oversight.** `MC_BATCH` replaced the
lookup keys (`out[[MC_BATCH]]`, `prov(args, MC_BATCH)`) because a typo there silently yields a
wrong column. It cannot replace a `name =` position inside `list()` / `data.frame()` without
`setNames` or `do.call`, which is more machinery than the literal costs. The rule is: keys use the
constant, declarations spell it out.

---

## 2026-07-31 -- one name for the batch label, and finalizers re-finalize

Third pass over the finalizers, driven by walking the rbind -> `calc_accel` workflow. Always-on
suite 1036 pass / 0 fail / 260 skip. Parity not run; `R CMD check` not run.

**The batch label is `mc_batch_id` everywhere, renamed from `batch`.** Two reasons, and the second
is the one that forced it. First, it was about to become a formula variable in `calc_accel`, and
a user's `batch` is their **slides or plates** -- biology, and a covariate they may legitimately
want in the rhs -- while ours only says which samples shared a cohort-mean fill. Shadowing theirs
with ours would be the worst kind of silent wrong answer. Second, and independent of formulas:
`batch` is a common column name, so `clocks_coverage()` / `samples_coverage()` output collided with
a user's own metadata on a join. The `mc_` prefix fixes both, and the rename went everywhere the
user can touch it -- both coverage frames, both finalizer frames, and `$provenance` -- because
"prefixed here, bare there" is exactly the inconsistency being removed. `data =` supplying
`mc_batch_id` is a flat **error**: the name is reserved, which is simpler than a precedence rule
and, given the prefix, cannot happen by accident.

**The label is always minted, so it is never auto-injected into a design.** `construct_mc_result()`
sets one on every run, so a single-batch record's label is a *constant* column -- absorbed by the
intercept, numerically nothing, conceptually noise. And on a multi-batch record injecting a
k-level factor adds k-1 columns and moves **every** residual for **every** clock. Chunking is an
operational decision (memory, streaming); it must not silently change the model. So the label is
offered to the formula and never added to it.

**The prompt to use it is gated on the fill actually happening.** Not "more than one batch" but
"more than one batch **and** some CpG was cohort-mean filled". Verified detectable:
`score_imputed_partial` is 0 on a clean run and >0 once the beta matrix has NAs, per batch. Only
the *partial* fill counts -- `imputed_full` takes the clock's vendored reference, which is the same
constant in every batch. Without that gate the note fires on records where the batches are
numerically irrelevant, which is how a real warning gets trained away.

**Every finalizer re-finalizes; `rbind` still does not.** `calc_accel` was residualizing
per-batch reductions of a bound cross-sample clock (measured: `DNAmPhysAge` accel changes after
`refinalize_clocks()`, `Hannum` does not) and saying nothing. The rbind reason for staying hands-off
does not transfer: `do.call(rbind, ...)` recurses, so re-finalizing there would redo the work at
every intermediate step, but a finalizer is a **leaf** -- it hands back a frame the record cannot be
recovered from, so it must hand back the right numbers. `as.data.frame()` does it too, so the two
finalizers cannot disagree about the same clock on the same record. Wanting the per-batch
reductions is still expressible: finalize each record *before* binding.

**The guard is `say_pending()`'s, not `length(pending)`.** First cut used the latter and printed
"Re-finalized ..." on every single-batch record holding a cross-sample clock -- `calc_clocks()`
retains `pending` even for a single pass, and a single-batch reduction already spans its whole
cohort, so that call was a numerical no-op emitting a message. Matching `say_pending()`'s condition
(pending **and** more than one batch) is also the right symmetry: finalizers re-finalize exactly
where `rbind` would have warned.

---

## 2026-07-31 -- `calc_accel` internals: one QR per row set, and the join gates what a join can break

Second pass over `R/calc_accel.R`. Always-on suite 1020 pass / 0 fail / 260 skip. Parity not run;
`R CMD check` not run.

**`lm()` is gone, and with it the scratch response column.** The first cut fitted per clock with
`lm(.mc_y ~ <rhs>)`, which needed a fixed magic name for the response because the response is a
different clock each iteration, and needed the rhs spliced back into a two-sided formula via
`call("~", quote(.mc_y), formula[[2L]])`. All of that existed only to satisfy `lm`'s interface. But
the rhs is **already one-sided**, so it *is* the design formula: `model.matrix()` consumes it
directly and `qr.resid()` needs no response name at all. Verified identical to `residuals(lm())`,
max abs diff 0, over `~ Age + Female`, `~ Age + I(Age^2)`, a factor covariate, `~ 1`, and a
rank-deficient design. It also makes the degenerate test exact -- `nrow(X) - qr$rank < 1` -- where
before it was a `length(vars) + 1` heuristic (wrong for factors and interactions) backed by a
post-fit `df.residual` check.

**The grouped fit stopped being an optimization.** `dev/calc_accel.md` had deferred "group clocks
by missingness pattern, one QR per pattern" as a performance escape hatch. Under `qr.resid` -- which
takes a **matrix** response -- it is simply the shorter code, and in the common case (no `NA`
scores) there is exactly one group, so the entire residualization is one decomposition over the
whole n x k matrix instead of k separate `lm` calls each rebuilding the same design. Group on the
rows a column *fits over* (`!is.na(resp) & keep`), not on its raw NA pattern: two clocks whose NAs
differ only where `keep` is already `FALSE` share a design and must not get two QRs.

**The cost is accepted, not overlooked:** `qr.resid` returns residuals and nothing else -- no
coefficients, R^2, or p-values. The frame is for plotting, and anyone modelling on it re-adds
`Age`/`Female` to their own rhs. If a future caller wants a slope back, that is `lm` coming back
with it.

**The join gates what a join can break, and nothing else.** Proposed a strict coverage gate --
every record sample must appear in `data`, mirroring `resolve_pheno()`. Wrong framing, and dropped.
`merge_accel_data()` is a **left join**; unmatched left rows are what a left join *is*. The only
thing that makes it ill-defined is a duplicated key on the right, and since `$pheno`'s ids are
already unique by the scoring-time `check_pheno`, "validate 1:1" reduces to the one check that was
already there. So: **multiplicity is a gate, coverage is a report.** The unmatched-id case warns
instead, which also draws a distinction worth having -- an explicit `NA` row in `data` is the
caller saying they know, a *missing* row is the caller possibly not knowing. Same downstream `NA`,
different thing to say.

That report is not cosmetic. Zero id overlap (a `sample1` vs `Sample1` typo) previously produced
two true but misleading warnings -- "Age: 8 samples missing", then "too few complete samples" --
and sent the user looking for phenotype data that was fine.

**Also corrected:** the duplicate-id gate's stated reason. `match()` returns the first match and
cannot multiply rows, so "fan-out" was wrong; the gate exists because the pick would be silent and
arbitrary.

**Type families replaced the numeric/character free-for-all.** `values_agree()` compared
storage-agnostically across *everything*, so `45` and `"45"` agreed. Now `type_family()` returns
`number` (numeric/logical), `string` (character/factor), or the class name -- and a cross-family
pair is a **distinct error naming both classes**, not a disagreement count, because nothing
disagrees in value. Returning the class name for everything else means `Date`, `POSIXct` and list
columns each become their own family and mismatch for free, so no separate atomic guard is needed.
`TRUE` vs `1` still agrees, deliberately: that is `Female`'s storage, not its value.

---

## 2026-07-31 -- building the finalizers: `calc_accel`, no `collapse`, and what the record's pheno cannot answer

Built the surface designed in the entry below (`R/calc_accel.R`,
`tests/testthat/test-clocks-accel.R`). Five things came out differently once it was code, so the
entry below states the design and this one states the build. Always-on suite: 1008 pass, 0 fail,
260 skip. Parity not run; `R CMD check` not run.

**The verb is `calc_accel()`, not `clock_accel()`.** The finalizer family is already
`clocks_coverage()` / `samples_coverage()` -- plural noun, then the quantity. A singular `clock_`
would have been the only one, and the frame is one row per (sample, clock) anyway.

**No `collapse` dependency.** The plan added it for `pivot()` and `join(validate = "1:1")`, and
both evaporated: keeping the whole pipeline matrix-native (`$scores` in, an NA-filled matrix of the
same shape out, `shape_scores()` choosing long or wide once at the end) makes the long/wide exit
five lines of `rep()` and makes the residual rejoin **positional by construction** rather than a
runtime join whose 1:1-ness has to be validated. A dependency that buys nothing base R does not do
in three lines is bloat regardless of how much we like the package; `collapse` stays available for
a future path (`HDW()`) where it would earn its place. This also means `clock_id` is a character
column, not the factor `pivot()` returns.

**`check_pheno()` gained nothing; its NA hint got shorter instead.** The plan wanted a caller-noun
argument so the message could say "score NA" for `calc_clocks()` and "dropped from the fit" for
`calc_accel()`. Built that, then removed it: the two consequences are the same consequence
("Those samples will score NA"), and a format-string parameter to distinguish them is machinery in
place of one general sentence. The **missing-column** message is a different matter and does live
in `calc_accel()`, because `check_pheno()`'s hint ("add it to `pheno`") would send the caller
back to re-score when the fix is `data =`.

**The `$pheno` narrowing stands, and the reason is stronger than the one in CLAUDE.md.** Proposed
widening `resolve_pheno()`'s `keep` to retain every supplied column, on the grounds that discarding
a user's `Age` is a silent discard and that `~ Age + <cell counts>` is the realistic accel formula.
Declined, correctly. `$pheno` does not exist to remember what was fed in -- it exists to remember
**what went into the numbers**. `Age` and `Female` are on the record because the clocks consumed
them, which is exactly what gives `merge_accel_data()`'s conflict check something real to defend:
it can refuse a `data` that would residualize on a different `Age` than the one that produced the
score. Widen it and `$pheno` becomes a general covariate bag, the conflict check starts guarding
columns no score ever saw, and the footgun it exists to prevent comes back. A covariate the scoring
did not use is the caller's to carry, and `data =` is where it goes. (The costs of widening, for
the record, were also real: a new `rbind` gate on pheno columns, and `print.mc_result` cutting the
pheno block instead of printing every column. Neither was the deciding argument.)

**So `calc_accel(res)` failing on a covariate-free record is the design, not a papercut.** The
error names `data =` and that is the whole fix.

**`type = "diff"` with an rhs that spans `Age` is exactly `type = "accel"` on that rhs.** Subtracting
a column of the design from the response changes the coefficients and leaves the residuals
identical, so `diff` + `~ Age` and `accel` + `~ Age` agree to floating point. The
constrained-vs-estimated distinction the design turns on is real but only bites when the rhs does
**not** span Age -- `~ Female`, say, where `diff` fixes the age slope at 1 and `accel` never sees
Age at all. A test asserted the two always differ and failed, correctly. The test now pins both
halves: equal on `~ Age`, different on `~ Female`. Do not "fix" the equality case by special-casing
`diff` -- it is linear algebra, and a `diff` that disagreed with `accel` there would be wrong.

**Also:** `expect_warning(x <- expr)`, never `x <- expect_warning(expr)`. testthat returns the
*condition* when one is caught, so the second form binds a condition object; `acc$SomeClock` is
then `NULL` and `all(is.na(NULL))` is `TRUE`, which is how a degenerate-clock test passed
vacuously here before it was caught.

---

## 2026-07-31 -- the finalizer family: no `augment`, `clock_accel` instead, `data` adds but never changes

**Design decision, nothing built.** Recorded now because the shape was argued out in full and the
reasoning should not be re-derived. Ideas were collected from PR #3 (`dev/pr3-triage.md` sec 4.4,
D1-D3); **no code is being taken from it** -- this is a clean re-implementation, and the surface owes
the PR about two things (the `na.exclude` hazard, and clash detection on a user covariate frame).
This entry settles `dev/pr3-triage.md` sec 5.4 for D1 and D2, which are in. D3 (`codebook`) was
kept out here and that was **reversed on 2026-08-04** -- see the entry of that date; do not read
this paragraph as the current position.

**`mc_result` is the canonical class, so every data.frame-returning function is a *finalizer*** --
a one-way exit. Past it you have no `rbind`, no `refinalize_clocks()`, no coverage, no provenance,
and re-binding what you got is the caller's problem. This is not a new contract: `clocks_coverage()`,
`samples_coverage()` and `as.data.frame(cite_clocks(x))` are already finalizers by this definition.
Naming it is what makes `as.data.frame.mc_result` an ordinary member of an existing family instead
of a novel API decision -- and it means **no finalizer needs a guard**. No `force =`, no round-trip,
no re-attaching a class on the way back.

**No finalizer's output contains another finalizer's output.** Each returns join keys plus its own
payload; the caller joins. This is the rule that kills `augment()`, and the argument is not merely
that `augment(x)` duplicates `as.data.frame(x)` -- it is that the nesting *grows*: add `_resid`, a
second `adjust` set, a `_diff`, and the frame doubles each time while the id and score columns repeat
verbatim. Applied to the derivative verb it also settles what comes back: `clock_accel()` returns id
+ clock + acceleration and **not the covariates it fit on**. `Age` and `Female` are inputs, not
payload, and the caller already holds the frame they came from.

**`augment` was also the wrong name twice over.** `generics::augment` is a real S3 generic
re-exported by broom and every tidymodels package, so a bare `augment <- function(x, ...)` is masked
the moment a user attaches broom afterwards, and the call then dispatches to broom's generic, finds
no `mc_result` method, and errors -- the same collision that made `cite_clocks()` a package-owned
name (2026-07-23, 2026-07-24, 2026-07-25). Beyond the collision, broom's `augment` is singular
because there is one model and one `.fitted`; here the derived quantity is plural and parameterized.
Splitting the join from the derivative lets the derivative be named for what it does.

**`type` sets the response; `formula` sets the right-hand side.** That factorization is the whole
signature: `type = "accel"` with no formula is `resid(score ~ Age)`, the classic; with a formula it
is the residual on that RHS; `type = "diff"` is `score - Age` with no fit at all; `type = "diff"`
plus a formula is `resid((score - Age) ~ ...)`. The last is **not** redundant with the second and
must not be "simplified" away -- constraining the age coefficient to exactly 1 is a different
estimator from estimating it, and the constrained one is a real (and contested, re: regression to
the mean) choice in this literature. `diff` is expressible as a formula in principle
(`resid(lm(score ~ 0 + offset(Age)))` is exactly `score - Age`) and is not offered that way: nobody
reads `~ offset(Age)` as a subtraction, dropping the `0 +` silently yields something else, and a
projection routine will not honour `offset()`. The structural reason is stronger than the ergonomic
one -- **accel is cohort-relative and diff is per-sample**, so diff alone is stable under subsetting
and under `rbind`, and that is the one distinction a reader must see at the call site.

**`type = "diff"` with no formula is identity, not `~ 1`.** `resid(d ~ 1)` is `d - mean(d)`, which
would silently re-introduce cohort dependence into the one arm that does not have it.

**`data` may add a column; it may never change one.** `$pheno` is `unique(c(pheno_id, covariates))`
where the covariates are the ones the run *actually required*, so every non-id column of it is by
construction a value that was fed into scoring. A `data` column that disagrees therefore disagrees
with a scoring input, and residualizing on an `Age` the score was not computed from is incoherent.
Three consequences, no clauses: absent from `$pheno` -> added; present and equal -> silent, `$pheno`
wins; present and different **including NA-vs-value** -> error, pointing at `calc_clocks()`.

The first cut allowed `data` to fill where `$pheno` was `NA`, on the argument that such a sample
either scored `NA` anyway or its clocks never read the covariate. **Rejected: under that rule the
same value has two different roles in one call, and which role depends on which clock you look at**
-- for an Age-needing clock the fill overrides an input that produced `NA`, for an Age-free clock it
supplements a score that never touched it. No rule stated in one sentence can mean both. Scoring is
cheap enough that "re-run `calc_clocks()` with the corrected pheno" is the honest answer, and it is
the same line `rbind` already takes: bind and label, never reconcile, and no `force =`. Nobody is
prevented from anything -- a caller can rebuild the frame in two lines -- we decline to launder it
through the machinery.

Two details this rests on. **The comparator is tolerance-based and storage-agnostic, never
`identical()`**: a pheno re-read from CSV gives integer `50L` where the in-memory one had double
`50`, and `identical()` would hard-error the modal workflow (score with `pheno`, add `Female`, pass
the same frame as `data`) over a difference that is not one. This is the same reasoning as
"never `expect_identical()`, always `expect_equal()`". And **the conflict check is *not* scoped to
the formula's variables** even though everything else is -- a disagreeing `Female` means some scores
were computed from a different `Female` whether or not this call reads it, and scoping would make one
`data` frame acceptable in one call and refused in another, which is exactly the clause-dependent
behaviour the NA-fill rejection removed.

**Two NA drops, different scopes, and neither is `na.omit()`.** The covariate set is exactly
`all.vars(formula)` (plus `Age` for `type = "diff"`) -- `na.omit()` on the frame would drop rows
whose `Female` is `NA` on a `~ Age` call, where `Female` is irrelevant. The pheno drop is **global**
(one eligible sample set, applied before the pivot, so no clock is ever fit against a covariate row
another clock did not see); the score drop is **per-cell**. Order is: drop pheno NA, pivot long, drop
score NA, inner join. That is the same cell set as one `complete.cases()` over the long frame, and it
is written as two steps because the two counts are separately meaningful and separately reported.

**Fit sets are per-clock and are deliberately not unified.** A clock with `NA` scores is fit on
fewer rows, which is honest -- its residual is genuinely relative to the vector that exists. The
alternative (drop any sample `NA` in *any* clock, so residuals are cross-clock comparable) is
rejected because one degenerate clock -- a pack clock `NA` for half the cohort -- would then silently
shrink every other clock's fit set. One bad column must not move Horvath's residuals.

**Consequence to document, not to fix: `accel` is `NA` on a superset of where `score` is.** Worked
example, n = 100, 5 samples with `Age = NA`, GrimAge (needs Age) also 3 coverage-NA, Horvath 2
coverage-NA. Both drop the same 5 globally; GrimAge fits 92, Horvath 93. Horvath's acceleration is
`NA` for 7 samples but its *score* is `NA` for only 2 -- the other 5 have a good score and no
residual because their covariate was missing. A reader seeing a score with no acceleration will
otherwise assume scoring failed.

**Drop-upfront makes `na.action` moot, so use `na.fail` as a free assertion.** The classic
misalignment trap is that `resid()` returns one value per complete case and assigning it into a
full-length column shifts everything; `na.exclude` patches that by padding back out. Having already
dropped the incomplete rows, the residual is length- and order-matched by construction, so nothing
depends on `na.action` semantics -- and setting it to `na.fail` turns any NA that slips through into
a loud error instead of a silent shortening. The rejoin is a left join on `(id, clock_id)` with
multiplicity validated `1:1`; a fan-out is the one remaining way this pipeline could misalign.

**Degenerate clocks yield `NA` and one aggregated warning.** All-NA (zero rows after the drop, `lm`
errors) and `n <= p` (`lm` returns aliased `NA` coefficients rather than erroring, and the residuals
are meaningless) are both caught by a row-count check before the fit. One warning naming the affected
clocks -- 121 separate warnings would be unusable.

**`check_pheno()` is reused rather than re-implemented**, with `extra_columns = all.vars(formula)`.
It already errors on a missing required column, type-checks `Age` (finite numeric) and `Female`
(integerish 0/1) for free, and warns with a per-column count of NAs over the id-joined rows -- which
is the global pheno drop, reported without new code. Two of its strings are clock-specific ("columns
the requested clocks need", "Clocks that need them will score NA") and need a parameterized noun,
safe to change since the test rule forbids asserting wording.

**`long = TRUE` is the uniform default across the family**, including `as.data.frame.mc_result`;
`as.matrix()` is the wide exit for scores. Per-function "natural" defaults were considered and
dropped -- a shared argument whose default flips per function is worse than no shared argument.
Long is also the computation form: one response name makes `reformulate(vars, response = "score")`
uniform across every clock, the per-cell NA drop is a single row filter, and the id survives so the
rejoin is exact.

**Fit N is derivable, so it is not stored.** It is the count of non-`NA` `accel` per clock. Storing
it would breach the no-nesting rule and it cannot be a column in the wide shape anyway.

**`collapse` is accepted as a dependency** -- `pivot()` for the long/wide exit, `join()` for the id
joins, and `HDW()`/`fHDwithin()` as the eventual fast path for linear residualization. Its NA/`fill`
semantics must be verified against a hand-rolled `lm` on a small case before it is used for the fit;
until then `split` + per-clock `lm` is the baseline. That loop rebuilds the model matrix once per
clock, which at this package's scale is about a second and not worth optimizing -- and the escape is
not a rewrite, since `lm` accepts a **matrix response**, so grouping clocks by distinct missingness
pattern gives one QR per pattern (one fit total in the common case) through the same
`residuals()` call.

Still open: nothing blocking. `type = "diff"` needs `Age` spelled exactly that way (the catalog's
only spelling; 25 clocks require it, 16 require `Female`) and gets its error from `check_pheno`
rather than a bespoke check.

## 2026-07-31 -- the paper's fields ride `mc_citations`; the runtime never parses BibTeX

`sync.R` now parses the vendored `clocks.bib` and widens the citation join with the paper's own
fields (`title`, `author`, `year`, `journal`, `volume`, `number`, `pages`, `doi`, `url`), so
`as.data.frame(cite_clocks(x))` is 13 columns instead of 4.

**No information was ever lost -- only its shape.** `bib_entries()` already read every entry
verbatim and `toBibtex()` already round-tripped it complete, fields and all. What did not exist was
a *tabular* view: `$links` is the clock -> paper join (`clock_id`, `pmid`, `role`, `bib_key`), which
answers "which paper" and not "what is that paper". Anyone wanting to sort by year or pull a DOI had
to regex the BibTeX text themselves.

**Parsed at sync, not at read.** The alternative was an `as.data.frame(x, expand = TRUE)` that ran
the extraction on demand. That puts a regex parser in the read path to recover data the build
already had in hand, which is the shape of thing "accessors read declarations; they never search"
exists to prevent. The `.bib` is vendored at sync time, so parsing it there keeps every runtime read
a lookup. Cost is 9 columns x 128 rows denormalized across 41 distinct papers -- immaterial in
`sysdata.rda`.

**The layout is read as declared, never searched.** `clocks.bib` comes from upstream's
deterministic emitter (`scripts/lib_bib.py`, `write_clocks_bib` / `format_entry`), gated by
upstream's own round-trip probe: one field per line, two-space indent, ` = {` assignment, values
collapsed and never wrapped. `read_bib_fields()` parses exactly that and `stop()`s on any deviation
-- unindented, over-indented, wrapped, unclosed, duplicated field, duplicate key, missing required
field. The first cut regex-mined each field out of the entry blob and was replaced: **a tolerant
per-field search cannot tell "field absent" from "layout moved"**, and `volume` / `number` are
legitimately absent on some entries, so drift would have shipped a table of NAs instead of failing.
This is the same rule as "accessors read declarations; they never search", applied at sync.

**Three ways it fails loudly.** A cited `bib_key` absent from the `.bib`; the two independent `pmid`
copies (`clock_citations.csv` and the entry) disagreeing; and any field in `BIB_FIELDS` appearing on
**no** entry at all. That last one is not redundant with the per-entry required-field check:
`BIB_REQUIRED` is only `title`/`author`/`year`/`pmid` because the entry *type* is deliberately not
asserted -- a `@misc` may legitimately carry no `journal` -- so a wholesale rename upstream
(`journal` -> `journaltitle`) passed every per-entry gate and produced an all-NA column. Measured:
it silently NA'd all 43 entries. Presence-somewhere is the strongest claim that holds without
pinning the type. Contrast PR #3's `bibliography()`, which synthesized a placeholder
`@article{key, pmid, url}` for a missing key and so hid exactly the gap worth hearing about.

The strict rewrite is **output-identical** to the regex version on the shipped bibliography
(`identical()` on all 128 rows), so it needed no `sysdata.rda` regeneration.

`pmid` is parsed but not attached -- the join already carries it, and a second copy would be a
column that can disagree with its neighbour. `volume` (6 rows) and `number` (10) are `NA` because
those papers declare neither upstream; every other field is complete across all 128 rows.

---

## 2026-07-30 -- The batch label is derived from the sample ids, and `batch =` is gone

Reverses this same day's "labels are assigned, never derived" and the `CLAUDE.md` line "never hash
anything to make one". `$pheno` is now always materialized, `calc_clocks(batch =)` is removed, and
`construct_mc_result()` derives the label as `batch_hash(pheno[[pheno_id]])`.

**What the assigned label cost.** The sticky-vs-auto rule needed `is_auto_label()`, which
re-derives the caller's *intent* from the label's spelling (`^[0-9]+$`). Nothing in a string can
carry that, so `batch = "2024"` twice silently renumbered to `2024, 2025` while `"T1"` twice threw
-- a wave/year label is the obvious thing to pass and the one that misbehaved. The collision error
also dead-ended: it said "give argument 2 a different name", and doing so hit
`apply_arg_name()`'s refusal to rename a non-auto label, with no way out but re-scoring. Both are
symptoms of storing a label without storing whether a human chose it.

**Deriving it removes the policy rather than fixing it.** A label that is a function of the
record's own ids is stable under re-association by construction: `rbind(rbind(r1, r2), r3)` and
`rbind(r1, r2, r3)` now return `identical()` records, where before the right-hand side renumbered.
So `rbind` mints nothing, renames nothing and renumbers nothing -- `DEFAULT_BATCH`,
`is_auto_label()`, `check_batch_label()`, `next_auto_label()`, `apply_arg_name()`,
`resolve_labels()` and `batch_maps()` all deleted, and the old "sticky name collides loudly" error
with them.

**`rbind` argument names are dropped, not refused.** Refusing them was tried first and reverted the
same session: `split()` names its result by factor level, so
`do.call(rbind, lapply(split(seq_len(n), g), score))` -- the canonical blocking idiom, the workflow
this feature exists for -- died on a wall of `"1", "2", ... "50"`. That is the `is_auto_label()`
mistake again in a new place: a name on `...` no more carries "I meant to label a batch" than a
digit-shaped string does. Names arrive from `split()`, `setNames()` and `Map()` for reasons that
have nothing to do with batching, and those cases outnumber a hand-typed `rbind(early = r1)` by far.
There is no way to tell the two apart -- `do.call` builds a call carrying the names either way --
so a warning would fire mostly on correct code. `unname(list(...))` and nothing else.

**Why the id column and not `$pheno` whole.** Hashing the pheno frame was the original proposal
and is wrong twice over. Before this change `resolve_pheno()` returned `NULL` whenever no pheno was
supplied -- 95 of 120 callable clocks require no covariate -- so *every* such batch hashed
identically and k batches would have collapsed onto one `per_clock` key, silently merging exactly
the per-batch fill regimes the batch axis exists to keep apart. Materializing `$pheno` fixes that,
but hashing the whole frame then folds covariate *values* into batch identity: correcting one
subject's age renames the batch, and `digest` is sensitive to column storage type, so the same CSV
read with integer vs double `Age` hashes differently. The id column has neither problem and answers
the question the label is actually asking -- *which samples were scored together*. In the 95/120
no-covariate case the two are byte-identical anyway.

**Rejected reasons, revisited.** 2026-07-30 rejected id-set hashing on cost, not soundness (it
says outright "hashing the id set is the only sound one"): a `digest` dependency, hex in the two
frames people read, and the `batch_set_id` non-goal. The first is paid (`digest` Suggests ->
Imports; it was already a Suggests). The second is handled by moving `batch` to the **end** of
`clocks_coverage()` and `samples_coverage()` -- it is still the key those frames join on, it just
no longer sits in front of `clock_id` where it reads as noise. The third stands as written and is
the part being reversed: the ban existed to stop anything *joining* on a content-derived id, and
`samples_coverage()` already joins on `batch` regardless of how the label is produced, so deriving
it changes what the label costs to compute and not what it is used for.

**The hashed value is canonicalized, and that is not optional.** `digest(ids)` hashes the id
*sequence and its R representation*, not the id set: reversing three ids, naming them, or handing
in a factor each produce a different label, so re-scoring a block after sorting its rows would
silently relabel it. `batch_hash()` therefore hashes
`paste0(sort(unname(as.character(ids)), method = "radix"), collapse = "\r")`. Every piece is
load-bearing -- `sort` makes it a function of the set, `method = "radix"` keeps that
locale-independent (plain `sort()` collates per-locale, so two machines would disagree),
`unname`/`as.character` drop attributes, and `serialize = FALSE` hashes the bytes so no R
serialization version reaches the label. A batch label appears in published
`clocks_coverage()` output, so it has to survive a change of machine.

Width is 12 hex of `xxhash64` (48 bits). Gate 1 makes the id sets disjoint, so a shared label needs
a real collision: **1.8e-11 at 100 batches, 1.8e-9 at 1000** (birthday bound; an earlier revision
of this entry said 1e-13, which was wrong by ~100x). 200k distinct inputs collided zero times when
measured. Since there is no longer a gate on labels, a collision would be silent -- the trade is
accepted at these numbers, and the fix if it ever mattered is to widen the truncation, not to
re-add the gate. `xxhash64` is non-cryptographic, which is fine here because sample ids are not
adversarial.

**Two unreachable guards were cut rather than kept "just in case."** A `gate_distinct_batches()`
and a partial-`pending` "package bug" stop were both written and both deleted the same session,
39 lines between them. Neither can fire: the first needs the 48-bit collision above, the second
needs two records with equal `$provenance$clocks` and different `spec$cross_sample`, which cannot
happen because `covariates`/`cross_sample` are functions of the clock sequence that gate 2 pins
(measured: `GrimAgeV1` vs `c("GrimAgeV1", "DNAmADM")` give identical output columns *and*
identical pending keys). The rule they now follow is the one that also deleted `gate_same_pheno()`'s
column check and rejected a `covariates_used` gate: **a guard on something an earlier gate makes
impossible is deleted, not kept.** Keeping some and cutting others left no statable principle, and
"which unreachable guards do we keep" is exactly the per-guard re-litigation the "one line decides
every gate" framing exists to stop. If a cross-sample clock ever becomes a routing target, gate 2
stops pinning `cross_sample` and the second guard becomes reachable -- add it back then, with a
test that fires it.

**`$pheno` is always present now.** With none supplied it is the id column alone, which reads
*closer* to the standing invariant ("the aligned pheno narrowed to the id column plus the
covariates the run required") than `NULL` did -- `NULL` carries no id column at all. This also
kills two things: `print.mc_result`'s `is.null(pheno)` branch, and `gate_same_pheno()`'s carried-
column check. That check is now unreachable -- `names(pheno)` is `unique(c(pheno_id, covariates))`
and `covariates` is a pure function of the clock sequence, which gate 2 pins -- so it is deleted
for the same reason the "same id, different covariates" gate was deleted this morning, leaving
`gate_same_pheno_id()`.

**`sim_DNAm(batch =)` is renamed `suffix =`, not removed.** It never labelled anything: it
suffixes the sample ids (`sample1_T1`) so two simulated blocks are disjoint by construction and
clear gate 1. That job is still needed and is now the *only* thing the caller controls about
batching, so it stops sharing a word with the derived label. `$batch` on the returned `mc_sim`
becomes `$suffix`.

---

## 2026-07-30 -- Phase 4 gates: refuse what the caller chose, record what batching forced

Settles the `rbind` design, which then shipped in the same session (`R/bind.R`, `test-bind.R`).
Narrows sec 8's four gates to four different ones and reverses two things sec 8 said.

**The line that decides every gate.** Sec 8 already said "record, never refuse, on differing fill
regimes", but did not say what that generalizes to, so each new question re-litigated it. It is:

> **Record what batching forces. Refuse what the caller chose differently.**

A per-batch fill regime is *forced* -- cohort means are per-run by construction, so a batched user
cannot avoid it and refusing would refuse the whole feature. Every other difference between two
records is a free argument with one obviously right answer across batches, so a difference is a
mistake and the bind says so. An earlier framing, "refuse on ambiguous identity, record on differing
method", was rejected: it puts `normalize=` on the record side, and `normalize=` is exactly the
caller choice that should be refused.

**The gates, and why there are four of them and not sec 8's four.**

1. **Disjoint ids**, checked as `anyDuplicated()` over the concatenated `$provenance$sample_id` --
   **not** over the pheno id column. `$pheno` is `NULL` whenever the run required no covariates,
   which is most runs, so a pheno-side check silently checks nothing exactly where the exposure is
   worst.
2. **Identical score columns**, reordered to the first record or thrown. This **subsumes sec 8's
   gate 3** (comparable coverage denominators): `output_ids` is the compute sequence minus routed
   members and the sequence is a deterministic function of the requested ids, so identical columns
   implies identical sequence implies identical panels. One gate, not two.
3. **Identical `pheno_id`.** Sec 8's gate 4 also required "no id appearing twice with different
   covariates", which is dead as written -- gate 1 forbids an id appearing twice at all.
4. **Identical `$provenance$normalized`.** New. Two records with the same clocks but different
   `normalize=` pass gates 1-3 with identical columns and identical *scoring* panels, while
   `Horvath1` measured 0.114 vs 7.715 absolute against the oracle depending on the setting
   (2026-07-29). The experimenter comparing normalized against raw is **already refused by gate 1**
   -- same samples scored twice, so the ids collide -- and wants the two columns side by side
   anyway, not stacked under a batch label. So the only caller who reaches this gate is one whose
   loop varied the argument across chunks, which is a bug.

**Refusing on overlapping ids is not over-strict.** A warning does not help: the double-count lands
in the *mean and sd*, so a z-score or an age-acceleration residual over the bound record is wrong
for every sample, not just the duplicated ones. After the bind, "one sample scored twice" and "two
samples sharing a name" are indistinguishable, so there is nothing to record instead. The fix is one
line user-side (`rownames(DNAm) <- paste0(rownames(DNAm), "_T1")`) and is the same relabelling the
caller's own downstream analysis needs. No `force =`, and no comparing record *contents* to detect
that two arguments are literally the same record -- the id sets are the whole test.

**Batch labels: no hash, and no stored counter.**

Hashing was considered in three forms and all three are rejected. Hashing `$pheno` cannot work:
`resolve_pheno()` narrows it to the id column plus required covariates, so it is `NULL` or exactly
`data.frame(ID = ids)` for any request with no covariates, and every such batch hashes the same;
where it *does* differ it folds covariate values into batch identity, so correcting one subject's
age silently renames a batch. Hashing the **fill regime** (`partial_fill` + `usable_cols`) is
float-brittle -- 2026-07-24's chunk-invariance measurement already found last-bit drift from
`dgemm` blocking alone -- and collides on the case that matters, two clean matrices over the same
CpG set. Hashing the **id set** is the only sound one (gate 1 makes collision impossible), but it
promotes `digest` from Suggests to Imports for a naming problem, prints hex into the two frames
people read, and is a `batch_set_id` -- the one thing the plan's non-goals name.

A stored counter field is rejected for a different reason: `$provenance$batch` is already a
per-sample vector, so the labels present **are** the counter and the next one is
`max(as.integer(labels)) + 1`. A separate field's only possible behaviour is to drift out of sync
with the vector it summarizes.

So the label is assigned, not derived:

- a user-supplied name (`calc_clocks(batch = "T1")`, or the `rbind()` argument name, which
  `do.call(rbind, named_list)` supplies for free) is **sticky** -- never renumbered, and two records
  claiming the same name throw, because that collision was deliberate;
- otherwise the next free integer, **renumbered on collision**.

`sim_DNAm(batch =)` is the same word doing a different job and is worth not confusing: it
**suffixes the sample ids** (`sample1_T1`) so two simulated batches are disjoint by construction,
and it defaults to `NULL` rather than to an integer, because a default would rename every sample in
every existing call. Threading it also collapsed a latent bug -- `ID` and the matrix rownames were
two independent `paste0("sample", seq_len(n))` expressions that happened to agree.

**Renumbering is safe only because the label is not a key.** `rbind(rbind(r1, r2), rbind(r3, r4))`
brings two sides both carrying batches `{1, 2}`; the right-hand side is renumbered to `{3, 4}`. That
shifts a name and never a row's group, so the partition is preserved -- which is exactly what the
"a batch label is not a `batch_set_id`" non-goal buys. The left-to-right folds
(`Reduce(rbind, ...)`, `do.call(rbind, ...)`, flat `rbind(r1, r2, r3)`) never renumber at all, and
a user who cares about stable labels names them.

**Two things in sec 8 were wrong and are rewritten, not annotated.**

- "Chunk reassembly labels every row one batch" predates the park. It was true when chunk
  reassembly meant the Phase 6 front end, which shares one fill regime across blocks. Post-park
  `rbind` **is** the chunking path, so k records are k batches -- which is precisely why the
  per-block offset has to be recorded.
- Sec 8's gates 3 and 4 are subsumed and half-dead respectively, per above.

**Coverage nests by batch.** `$coverage$per_clock` becomes batch -> clock -> record on **every**
record, single-pass included, and `clocks_coverage()` becomes one row per (clock, batch) with
`batch` leading. Merging was rejected: `score_imputed_partial` counts the panel CpGs in *that run's*
partial cache, and two independently-scored batches almost never share an NA pattern, so a merged
figure would be wrong on nearly every bind rather than in a corner case. The counts stay CpG counts
on the CpG axis; what changed is that the CpG axis is now indexed by batch. This is also what
satisfies sec 8's "per-batch imputation summary" without a hash -- two batches with different fill
regimes are visibly different in `clocks_coverage()` whether or not their labels say so.
`samples_coverage()` carries a `batch` column too. It is redundant against `id` (gate 1 makes the
id determine the batch) and exists anyway, because without it the long frame cannot be joined to
`clocks_coverage()` on the key that frame is now on.

**Not defended against: pack version drift.** Two records scored against different pack payloads
have identical columns and different coefficients, and it is undetectable -- `payload_hash` is
maintainer-side and deliberately never reaches a record (2026-07-24). Reversing that to make one
gate possible is not worth it; a package that moves breaks reproducibility by many other routes.

**Retaining `pending`.** It lands on `$provenance$pending`, keeping the record's four top-level
elements as `CLAUDE.md` states them. It is built **with** `rbind`, not before: retention only pays
off at bind time, so landing it early would add a field to the public record with no consumer.

## 2026-07-30 -- Phase 6 is parked; `rbind` covers the realistic case

Reverses the entry below it, which made the store mandatory and Phase 6 next. **Nothing in sec 5 is
deleted** -- it stays as the record of what was measured -- but the streaming front end is no longer
scheduled, and Phase 4 takes its place as the next thing built.

**The projection number, measured off the committed catalog.** 87 bundled callable clocks: panels
sum to 32,228 and union to **20,430** scoring CpGs, or 38,540 including the 21,368-probe BMIQ norm
panel. Against an 866k array:

| n | full array | projected union | x3 copies |
|---|---|---|---|
| 1,000 | 6.45 GB | 0.15 GB | 0.46 GB |
| 10,000 | 64.5 GB | 1.52 GB | 4.57 GB |
| 50,000 | 322 GB | 7.61 GB | 22.8 GB |

So a request for **every bundled clock is resident past any cohort that exists** on the 4-8 GB box
sec 5 targets. The entry below drew the same conclusion and then built for the exception anyway.
(The PC-scale panel sizes it quotes -- 357,852 for PCBrainAge and so on -- were not re-verified
here; the catalog does not carry `n_cpgs` and the packs were not staged.)

**The premise that fell is gzip.** The store was mandatory because a `.csv.gz` cannot be seeked, so
two passes meant two whole-stream deflates. But anyone with a cohort genuinely too big to hold does
not have a `.csv.gz` -- they have HDF5, Zarr or TileDB, because that is what an append-capable store
is. On a random-access source two passes are free and the store has no remaining justification, and
with it go the ingest contract, the lifecycle/consent surface for an 80 GB user-data object with no
default location, and sec 5.8's blocking benchmark.

**What is actually left, and why it is niche.** One case: a cohort too big to hold *even projected*,
**carrying NA**. The NA is load-bearing -- cohort-mean fill is the only thing sample-blocking can get
wrong, so on a complete matrix a user-side loop is already exact. That conjunction (too big to hold
projected AND missing values AND, in practice, a PC-scale panel) is rare enough not to lead a
roadmap.

**What covers the rest: Phase 4.** A user with a random-access store can already project
(`clock_cpgs()` ships), block, and score; the only missing piece is assembly. Two notes carried into
sec 8:

- Per-batch cohort means are a real difference but second-order as noise (block mean unbiased,
  SE `sd/sqrt(n_block)`). The mechanism worth naming is a probe all-NA within one block but partial
  cohort-wide: that block takes the vendored ref while others take the cohort mean, which is a
  systematic per-block offset, not noise. `rbind` cannot undo it -- the batch label is what makes it
  honest.
- **Cross-sample re-finalization can be exact, which sharpens the "do not re-finalize" rule rather
  than reversing it.** `physage_raws()` is per-sample, so if a record *retains* `pending` instead of
  discarding it after `finalize_cross_sample()`, a bind-time re-finalize reproduces the single-pass
  number rather than approximating it. Still opt-in and still never silent, since it rewrites a
  column the user has already seen, and it keys on `spec$cross_sample` rather than on PhysAge.

Also settled while measuring, and recorded because it outlives this reversal:

- **Bioconductor keeps phenotype out of the HDF5 by design.** `saveHDF5SummarizedExperiment()`
  writes `assays.h5` (assay datasets only, not even dimnames) plus `se.rds` holding colData,
  rowData and dimnames. Dimnames are a `writeHDF5Array()` convention -- `.<name>_dimnames/1` and
  `/2` plus a `DIMENSION_LIST` attribute -- and rhdf5 cannot read that attribute (`VLEN not yet
  implemented`), so the naming convention is the contract and the attribute is decoration. Any
  future HDF5 support adopts that layout rather than inventing one, and pheno stays an R argument.
- **The float32 non-goal is wrong as written.** "No float32 anywhere on the input path" justified by
  `PARITY_REL_TOL = 1e-10` conflates input storage precision with arithmetic precision; that
  tolerance bounds our agreement with an oracle given identical inputs. Refusing float32 input
  declines to score data whose precision was already fixed before we saw it. The rule that does the
  intended work: promote to float64 on read, never compute or store intermediates in float32.

## 2026-07-30 -- the duckdb store is mandatory, not scratch; the resident set is the panel, not the file

Sharpens 2026-07-29 rather than reversing it: duckdb stays the access engine, but the **ingest step
is promoted from a request-scoped `tempdir()` convenience to a required stage**, and the reason is
stated differently.

**The target this is built for.** A small box -- 2-4 cores, 4-8 GB RAM, i.e. a cheap AWS partition
or a laptop -- against a tall `.csv.gz` far larger than RAM, with no ETL the user has to run
themselves. On that box the run is **deflate-bound and serial**: gzip inflate is whole-stream and
single-threaded, so it is the floor no amount of cores moves. That is the whole reason chunking
exists here; if the data fits in memory the question is moot.

**Projection is the first lever, and it was missing from the plan.** Only `needed_cpgs x n_samples`
has to be resident, never the file. Against a 1e6-probe array at n=1e4 (80 GB as doubles) the 101
bundled panels union to 39,025 probes = 3.1 GB. So a 1e6 x 1e4 file scores on an 8 GB box with
*nothing* streaming, and what forces chunking is **requesting a PC-scale panel, not cohort size** --
PCBrainAge's 357,852 probes is 28.6 GB. The sizing table is in the plan (sec 5).

**Why the store is mandatory: gzip is whole-stream, a columnar store is per-page.** A `.csv.gz`
cannot be seeked, so every pass over it re-inflates and re-parses the entire file regardless of how
little that pass needs. Since cohort-mean fill makes two passes non-negotiable (plan sec 2), the
choice is between paying the deflate twice per request or ingesting once. The store is what converts
one whole-stream deflate into per-page reads.

**A per-request projected spill was considered and rejected.** Writing only the resolved panel,
partitioned by sample block, during pass 1 is smaller (3.1-28.6 GB against ~80 GB) and free as a
byproduct of that pass. It loses because the panel depends on the clock request, so it re-deflates
the source **once per request** -- and the realistic session is a sequence of requests over one file
(bundled clocks, then a PC clock added, then a re-run after a pheno fix). One deflate ever beats one
deflate per request, and 80 GB of scratch disk is the cheap side of that trade.

**duckdb, not parquet.** Maintainer's call on two grounds: duckdb handles the wide schema (1e4+
sample columns) better than parquet's per-column-per-row-group footer metadata, and a `PRIMARY KEY`
on the CpG id gives a real index plus an ingest-time duplicate-probe check. **This is a decision, not
a benchmark** -- parquet was not measured at that width, and neither was duckdb (see below).

**Five corrections found while deriving this**, all now in the plan:

1. **The block-width formula was off by ~3x.** Sec 5.3 divided the budget by one copy of the block.
   Peak is three: the block, `build_partial_cache()`'s slice, and `pack_design()`'s `n x |panel|`
   copy. The cache is not a sliver -- at any realistic NA rate nearly every column carries at least
   one NA, so it is a full second copy.
2. **Projecting first would silently redefine `sample_scale`.** `scan_missing_cpgs()`'s row moments
   span *every* column of the matrix, and the two Zhang2019 arms z-score each sample over the whole
   array; substituting the panel union was measured at 1.8e1 absolute / 82% relative. So pass 1's
   probe filter is on **iff `spec$needs_moments` is FALSE**.
3. **duckdb's `avg()` does not treat `+/-Inf` as missing; `col_stats()` does.** A SQL route to the
   per-sample moments therefore disagrees with the kernel on a file carrying an `Inf`, and needs the
   value gate ordered ahead of it. The kernel gives them in the traversal pass 1 already performs.
4. **`pheno` is not the id source** -- it is `NULL` for any request with no covariates. Sample ids
   are the header, which is a real simplification: under the tall orientation id uniqueness and the
   pheno subset are settled by a header-only read, not by the cross-block accumulation sec 5.4
   assumed.
5. **The expression-depth wall is about SQL arithmetic, not about width.** Returning 1e4 columns is
   fine (measured to 20,000); writing `S1+...+S10000` is not. Pass 1 asks for rows and sums them in
   `col_stats()`, so the limit never applies. Sec 5.3 said this; it reads as a workaround and is
   actually a refusal to use SQL for arithmetic at all.

**Recorded, not scheduled: probe-axis accumulation.** A pack-scored clock is a batched weighted sum,
and a weighted sum is additive over the contracted axis -- so accumulating `n x k` partial scores
over probe chunks needs ~10 MB resident at n=1e4 instead of a 28.6 GB projected panel, and would
remove sample-major blocking for exactly the clocks that need it most. The blocker is the branch
contract, not the arithmetic: "a branch returns only its score" would become "a branch returns a
partial score plus a merge rule", and every branch in the closed set would need an additivity
classification. It sits beside sample-blocking rather than replacing it, so it stays an open.

**Nothing here is measured at 1e4 samples.** Every figure above is arithmetic off the 2026-07-29
baseline (500 samples / 20k columns). Ingest carries ~10 ms per column of fixed overhead and 20,000
columns took 96.5s at only 5M cells, so a 1e4-column ingest sits in an overhead-dominated regime
nothing has been measured in. If it is bad, the fix is partitioning the store by sample block, which
changes the schema -- so that benchmark comes before the build (plan sec 5.8).
