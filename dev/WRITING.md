# Writing rules: roxygen prose and cli messages

**This is the single source for how user-facing text is written in this package.** It is
self-contained on purpose: an agent picking up one batch of topics should need this file and
nothing else.

**`CLAUDE.md` points here rather than restating.** It keeps the cli-versus-`stop()` audience line,
because that is a design invariant. Everything about how the text itself is written -- R1 to R8,
the cli mechanics, the roxygen template -- lives here and only here (2026-08-03).

This file is tracked.

---

## 0. Hard prohibitions

Binding on every agent, no exceptions, no "just once to check":

- **Never run `R CMD check` or `devtools::check()`**, by any entry point, including `Rscript`,
  `pkgbuild`, or a background shell. It hangs in this environment.
- **Never run the parity tier.** No `test_parity()`, no `MC_PARITY=1`. It is minutes long and is
  the maintainer's call.
- **Never hand-edit `NAMESPACE` or `man/*.Rd`.** They are generated. Own the tags.
- **Never add, remove, or retarget a `@seealso` link.** The groups are closed. See section 6.
- **Do not run `devtools::document()` or `devtools::test()` when working as one of several
  parallel agents.** Both write shared files and the runs race. The coordinating session runs
  them once after every batch is in.

---

## 1. Which channel: cli or `stop()`

**The line is audience, not transport.**

- A message about **input the user chose** is `cli`, whatever function raises it. This includes a
  message raised deep inside an internal helper, if the value it complains about arrived through
  an exported argument.
- A message about a **package defect** -- a "cannot happen" state, a catalog or sync gap, a
  missing dispatch branch -- is a plain `stop(..., call. = FALSE)`.

A defect message is **hard-coded and greppable**: a fixed prefix first, values appended after it,
so a bug report can be located from pasted text with no stack trace. Not necessarily a constant
string, but the fixed part leads.

Today: assets lifecycle, discovery printers, public S3 refusals, and the whole `calc_clocks()`
front door are cli. Accessors, score branches, pack dispatch, catalog bugs and citation internals
are `stop()`.

---

## 2. The English: R1 to R8

These bind **every word a user can see**: cli message text, roxygen prose, **`vignettes/*.Rmd`
and `README`**. They do **not** bind code comments, `stop()` text aimed at developers,
`data-raw/`, or `dev/` docs. In those, ASCII `--` is still required and these rules do not apply.

The prose files are the newest members of that set and the easiest to write carelessly, because
markdown invites a chattier register than a `@details` block. R2 and R3 are the two that slip:
no first person, no contractions, no `--`, no `;`. See section 10 for what else they must satisfy.

- **R1. ASD-STE100 Simplified Technical English.** One instruction per sentence. About 20 words
  for an instruction, 25 for a description. The simple word over the elaborate one. One word for
  one meaning. Articles present. No noun cluster longer than three words. No ambiguous `-ing`
  form. `coverage` and `superseded` are accepted exceptions, because both are public API names.
- **R2. No first person and no "please". No contractions.** Prefer the data or the object as the
  subject over a bare passive. A function name is allowed as a grammatical subject
  (`calc_clocks()` narrows `pheno` ...). Second person is allowed where the user really is the
  actor.
- **R3. No `--` and no `;` in user-facing text.** Period, comma, colon, and a single spaced
  hyphen. This is an accessibility requirement, not a style preference. A `;` is how a long
  sentence gets smuggled past a period: prefer a period and a short new sentence.
- **R4. Describe the problem and give an actionable next step.** Never state that scoring
  continues or that nothing was stopped. Name the function to call **next**, never the one
  currently running.
- **R5. No vector sized by user input reaches cli.** `MC_MSG_CAP` is 10, in `R/utils.R` with
  `capped_bullets()` / `capped_vals()` / `capped()`. There is no "... and N more" tail: the true
  total is already in the lead line. Where the cost is upstream of the render, cap the input
  count, not the rendered output.
- **R6. Every R language object in prose carries markup.** cli `{.arg}` / `{.fn}` / `{.val}` /
  `{.code}` / `{.path}` / `{.field}`; roxygen backticks and `[fn()]` links. **Only R language
  objects.** Domain and platform terms (`EPICv2`, `MSA`, "M-value", "CpG") stay plain prose,
  because markup that resolves to a cross-package link breaks when that package is absent.
- **R7. Keep a marked span short.** Mark the identifier, not the surrounding phrase and not a
  whole call with long arguments. roxygen2 renders a long backticked span badly and it breaks the
  PDF manual.
- **R8. No internal vocabulary.** Name an object the way the reader's own code names it. A word
  that is all over `CLAUDE.md` and `DECISIONS.md` but nowhere in the API is our word, not the
  reader's: **"record", "closed set", "the block", "routed member", "cohort-mean fill"**.
  The test is mechanical: if the word is not a function name, an argument name, a component name
  (`$scores`, `$pheno`, `$coverage`, `$provenance`), a column name in a returned frame, or a word
  already in a message the user sees, it is jargon. Say "the returned value" or "an `mc_result`
  object", never "the record".

  **"panel" passes the test and is allowed.** An earlier draft of this rule banned it. That was
  wrong on the rule's own terms: `samples_coverage()` returns a column literally named `panel`,
  and the count columns are `score_*` and `norm_*` against it. A word the reader sees in their own
  output is theirs, not ours.

### What no rule covers

Wordiness is not rule-shaped. "Packs already loaded by `load_mc_assets()` are used as they are"
breaks nothing above and is still bad; "Packs pre-loaded by `load_mc_assets()` are used directly"
is the fix, and it is a judgement. **The instrument for this is `dump_roxygen()`** (section 7).
Read the rendered page, not the tags. Both weak sentences in the first draft of `calc_clocks()`
were obvious on one reading of the output and invisible while writing the source.

---

## 3. cli mechanics

- **`sprintf` output must never become cli input.** cli parses every message and bullet as a
  template, so a `{` arriving inside a built string is read as syntax. Build the line with
  `cli::format_inline()` and hand data in as interpolated values, which are never re-parsed.
  `sprintf` is still right for a plain `stop()`, for `cli_verbatim()`, and for an `askYesNo()`
  prompt.
- **`bullets()` escapes braces**, because rendered text goes back through cli as a template.
  Every bullet path goes through that one door.
- **Cap first, then format** (`capped_bullets(x, fmt)`). At most ten elements are ever formatted,
  which is what makes per-element markup affordable.
- **`say_*` emits to the user. `note_*` records into the block's collector.** Never use `note_`
  for something that prints.
- **Bind every `{?}` plural marker with an explicit `cli::qty()`** unless the quantity is the
  interpolation immediately before it. cli binds a marker to the last interpolated value earlier
  in the same string, and the quantity does not carry across elements of a `c()` message vector.
  Get this wrong and the handler throws in place of the real diagnostic, or worse, silently
  pluralizes against the wrong value. Safe form:
  `"Add {cli::qty(need)}{?it/them} to {.arg pheno}."`
- **cli reflows whitespace.** A pre-aligned block collapses onto one line. Use `cli_verbatim()`
  where alignment matters; inside `cli_abort()` / `cli_inform()` bullets carry no alignment at all
  and emit one self-contained bullet per row.
- **Never hand `askYesNo()` a multi-line prompt.** Print context with cli first, then ask one
  short line. `mc_ask_yes_no()` in `R/mc_data.R` is the one door.
- **Tests assert *that* a message errors, never its wording.**

---

## 4. Roxygen: the template

### Tag order, every topic

```
Title Case Noun Phrase
(blank)
One sentence. What it does.
(blank)
@inheritParams mc-params   (only if it takes a shared param)
@param      one per remaining formal, in signature order
@details    only if needed
@returns    always
@examples   or @examplesIf, always
@export
```

- **Titles are noun phrases.** "The Default Scatterplot Function", not "Plot a Default
  Scatterplot". So "Epigenetic Clock Scores", not "Score Epigenetic Clocks". Title case, no
  trailing period.
- The sentence under the title is the description. It is a real sentence with a verb and a period.
  **No `@description` tag.**
- **`@returns`, not `@return`.** Both work; `@returns` is what this package uses.
- **No `@author`, no `@keywords`** on an exported topic.
- Markdown is on (`Roxygen: list(markdown = TRUE)`), so backticks and `[fn()]` links work.

### `@param` form

```
@param name <type>. <What it is, one sentence.> Default is <value>.
```

| kind | type fragment |
|---|---|
| matrix | `A numeric matrix.` |
| character vector | `A character vector.` |
| single string | `A string.` |
| threshold | `A number between 0 and 1.` |
| flag | `A boolean.` |
| count | `A single whole number.` |
| numeric vector | `A numeric vector.` |
| data frame | `A data.frame.` |
| list | `A list.` or `A named list.` |
| catch-all S3 argument | `Any object.` |
| a method that only throws | `Nothing.` (`@returns` only) |
| record | ``An `mc_result` object.`` |
| citation | ``An `mc_citation` object.`` |
| simulation | ``An `mc_sim` object.`` |
| formula | `A one-sided formula.` |
| named logical | `A named logical vector.` |
| enum | `One of "accel" or "diff".` |

**This table is `DOC_TYPES` in `R/dev-utils.R` and the linter checks against it.** If a genuinely
new kind of argument appears, add the fragment there in the same commit.

- **State the default only when there is one.** `Default is NULL.`, `Default is TRUE.`,
  `Default is "ID".` An argument with no default gets **no** default sentence, and never
  "Required." The linter checks both directions against the actual formals.
- **One sentence of description.** If it needs two, the second belongs in `@details`.
- **Union types:** a param taking more than two forms names its **primary** type, and a shared
  `@section` enumerates the rest. `ext_data` is the only one: type is `A string.`, and the three
  forms live in the donor's `The assets directory` section, pulled in with
  `@inheritSection mc-params The assets directory` by every topic that takes it. A `@section` on
  the donor works exactly like `@param` there, and is the right shape when the text is longer than
  a sentence and must not drift between topics.
- **`...` gets one of three fixed sentences:** `Passed to [calc_clocks()].` when it is passed
  through; `Not used.` for generic consistency; `Two or more mc_result objects.` for
  `rbind.mc_result`.
- **No `list(**item**: ...)` blocks.** If named values need spelling out, use a plain bulleted
  list with no bolded lead-in, each item a fragment. Prefer one more sentence to any list.

### `@returns`

Names a type from the same table, but is **not** held to the opening-fragment form. It describes
a value, so it reads as a sentence: ``An `mc_result` object. It holds the scores, the narrowed
`pheno`, the coverage counts, and the provenance of the run.`` The linter only checks that a type
appears somewhere in the text.

---

## 5. The shared-parameter donor

`R/mc-params.R` holds the text for every param used by more than one topic. A topic that takes
one writes **`@inheritParams mc-params`** and does not repeat the text. Currently in the donor,
nine: `DNAm`, `x`, `clocks`, `pheno`, `normalize`, `ext_data`, `ask`, `all_columns`, `groups`.
Re-read `R/mc-params.R` rather than trusting this list.

- **It must keep `@keywords internal`.** Measured: a `@keywords internal` donor generates an
  `.Rd` with no `\usage` and draws no check NOTE. A `@noRd` donor **breaks** -- roxygen resolves
  inheritance against topics, `@noRd` produces none, and the recipients end up with undocumented
  arguments.
- **A param joins the donor only if its default is the same, or absent, at every call site --
  unless a topic overrides it.** roxygen fills only the params a topic has not documented itself,
  so a topic with a different default writes its own `@param` and wins. `groups` works this way:
  the donor carries it with `Default is "all"`, and `load_mc_assets()`, which has no default,
  overrides it locally. That is worth doing when the shared text is long enough to be worth
  sharing, as the four asset-group names are; it is not worth doing for a one-line param.
- **The donor's `x` is an `mc_result`.** A topic whose `x` is something else --
  `print.mc_sim`, `cite_clocks.character`, the `mc_citation` methods -- **must not** write
  `@inheritParams mc-params`. Inheritance matches on name alone and yields confidently wrong text
  rather than an error. This is the one live footgun in the scheme.
- **Argument naming is settled: `x` stays `x`.** Renaming to `mc_result` was possible for only 5
  of 12 topics, because `print`, `as.matrix`, `as.data.frame`, `rbind` and `toBibtex` are locked
  to the base generic's argument names by the "S3 generic/method consistency" check. Uniform won.

---

## 6. `@seealso`: the groups are closed, do not add a link

The groups were decided once, with the whole surface in view, on 2026-08-03. **Do not add,
remove, or retarget a `@seealso` link on your own.** Not "obviously related", not "the same verb
pair". If you notice a pairing worth linking, **say so in your report** and leave the tag alone.

The reason is the invariant the set carries:

**A topic's `@seealso` is the union of the groups it belongs to, minus itself.** Symmetry then
holds by construction, including where a topic is in two groups. The closed set is:

| group | members |
|---|---|
| discovery | `list_clocks`, `list_clock_tags`, `clock_cpgs`, `list_mc_assets` |
| coverage | `clocks_coverage`, `samples_coverage` |
| record verbs | `as.data.frame.mc_result`, `as.matrix.mc_result`, `rbind.mc_result`, `refinalize_clocks` |
| analysis | `calc_accel`, `score_associations` |
| assets | the six `*_mc_assets*` verbs, including `list_mc_assets` |

`list_mc_assets` is in two groups, so it carries the union at eight links. That is the widest
list in the manual and it is intended.

**Untagged on purpose:** `calc_clocks` is the entry point and points nowhere, `cite_clocks`
carries its methods on one topic instead, and `sim_DNAm`, `predict_sex` and the three `print`
methods are left out by decision (DECISIONS 2026-08-03).

**Form:** one sentence for a single link, a plain bulleted list for two or more, each item
`[fn()]` plus a short clause saying what the reader gets there. The tag sits after `@returns` and
before `@examples`.

**A dangling link is an `R CMD check` WARNING**, and a one-sided link is invisible to every tool.
`lint_seealso()` (section 7) catches both. Run it after any change to a doc block.

---

## 7. Tooling

Both live in `R/dev-utils.R`, which is `.Rbuildignore`d, so `roxygen2` is a dev-time dependency
only and never reaches `DESCRIPTION`.

- **`lint_roxygen(path = ".")`** returns a data frame of violations, zero rows when clean. It
  checks: every `@param` opens with a `DOC_TYPES` fragment; a formal with a default states it; a
  formal without one does not; and `@returns` names a type. It reads the parsed call rather than
  the evaluated object, so it needs no `load_all()`.
- **`lint_seealso(path = ".")`** returns a data frame of cross-reference violations, zero rows
  when clean. It catches the two failures nothing else sees: a link whose topic does not exist
  (an `R CMD check` WARNING) and a link that is one-way (invisible to every tool). It reads
  `man/`, not `R/`, because `\link` targets only exist after `document()`.
- **`dump_roxygen(path = ".")`** renders every `man/*.Rd` through `tools::Rd2txt()` into one text
  file and returns the path. This is the read-through instrument for everything the linters
  cannot check.

```r
devtools::document()
lint_roxygen(".")            # must be empty
lint_seealso(".")            # must be empty
file.edit(dump_roxygen(".")) # then read it as a document
```

---

## 8. Worked exemplar

`calc_clocks()` in `R/calc_clocks.R` is written and is the cadence to copy. Read it before
starting. It shows the noun-phrase title, the one-sentence description, `@inheritParams` plus
three own params, three `@details` paragraphs, a typed `@returns`, and a runnable two-block
example.

### Examples

- **No `\dontrun{}`, ever.** Use `@examplesIf`, which renders to a `\dontshow{if (cond) ...}`
  wrapper: the example stays real, visible in the manual, and simply does not execute where the
  condition is false.
- **The six assets verbs get `@examplesIf interactive()`.** They prompt, download, delete, or
  print a per-machine path. `interactive()` is `FALSE` on CRAN.
- **Everything else runs unconditionally**, so it must work with no network and no packs. Build
  inputs with `sim_DNAm()` and **bundled clocks only** (`Horvath1`, `Hannum` are safe). Never
  `SystemsAge`, `PCClocks`, `PCBrainAge`, or the `Zhang2019` BLUP arm.
- **`sim_DNAm()` returns an `mc_sim` list, not a matrix.** The matrix is `sim[["DNAm"]]`, and it
  is **samples in rows, CpGs in columns**.
- No `set.seed()`. Nothing checks the printed numbers.
- About five lines for most topics. One `sim_DNAm()`, one call, one look.
- Examples do not share state across topics, so each block builds what it needs.

---

## 9. Auditing the manual: known-good exceptions

For an independent reader going through `dump_roxygen()` output or `man/*.Rd`. Everything in this
section is **intended and verified**. Reporting one of these as a defect is a false positive.

### The three shape rules, and where the package stands

Measured 2026-08-03 over all 28 topics. Re-measure rather than trust this line.

| rule | state |
|---|---|
| Every exported topic carries `@examples` or `@examplesIf` | 26 of 26 |
| No internal topic carries an example | 0 of 2 carry one |
| Every exported topic carries `@returns` | 26 of 26 |

The two internal topics are `mc-params` (the shared `@param` donor, section 5) and
`methylCIPHERv2-package` (roxygen-generated). Neither has `@returns` or `@examples`, and neither
needs them: **both are `@keywords internal` and generate no `\usage` block, so `R CMD check` does
not ask for either.** An exported topic missing either tag **is** a real finding.

### Known skips

0. **`vignettes/assets.Rmd` shows `download_mc_assets()` and `clear_mc_assets()` in code blocks
   that are `eval = FALSE`.** Every block that downloads, deletes, or prints a per-machine path
   is shown with its output written out by hand. This is the `@examplesIf interactive()` policy
   in the other syntax, and it is required: **the package builds with no network**, and a
   vignette is built by `R CMD check`. Only one block in that file evaluates, a `list_clocks()`
   call that reads the shipped catalog. Do not "fix" a block to run.
1. **`@examplesIf interactive()` looks like an unguarded example, and is not.** It renders as
   `\dontshow{if (interactive()) withAutoprint(\{ # examplesIf}`, and `Rd2txt` prints the body
   with the guard stripped, so the rendered manual reads as though the call runs every time. It
   does not. The six asset verbs use it -- `get_mc_assets_dir`, `set_mc_assets_dir`,
   `list_mc_assets`, `download_mc_assets`, `load_mc_assets`, `clear_mc_assets` -- because they
   prompt, download hundreds of megabytes, delete files, or print a per-machine path.
   `interactive()` is `FALSE` under `R CMD check` and on CRAN, which is the whole point. Do not
   "fix" one of these to run, and do not report it as an untested example.
2. **`cite_clocks` shows four `\usage` lines on one topic.** The generic and its three methods are
   merged with `@rdname` (DECISIONS 2026-08-03). One `@param x` covers both accepted input types
   by design.
3. **`@seealso` is a closed set.** See section 6. An apparently missing cross-reference is a
   decision, not an omission.
4. **`calc_clocks` carries no `@seealso`** although nearly every topic relates to it. It is the
   entry point, and the inherited `x` param text already links it.
5. **`--` and `;` are required in this file and banned in package messages.** The R3 ban is scoped
   to text a user can see. `dev/` docs, code comments and developer `stop()` text are outside it.

### What an auditor should actually check

Run both linters first (section 7); they cover the mechanical rules. Then read for the things no
linter sees: a `@param` sentence that does not match the formal it names, a `@details` paragraph
that describes behaviour the code no longer has, a title that is not a noun phrase, an example
that would need the network or an asset, and any breach of R1 to R8 in text a user can see. **R8
is the one with the most live breaches**, and no linter catches it, so give it the closest read.

---

## 10. The prose files: vignettes and README

`vignettes/*.Rmd` and `README.Rmd` are public-facing text. R1 to R8 bind them in full. Everything
below is in addition.

### The two files are not built the same way, and that sets most of the rules

This is the distinction to get right before writing a chunk. An earlier draft of this section
applied one rule set to both and would have forbidden the README that exists (2026-08-04).

- **`vignettes/*.Rmd` is built by `R CMD check` and on CRAN.** It runs on a machine with no
  network, no assets, and no staged cohort. Every constraint below that mentions offline
  execution comes from this and only this.
- **`README.Rmd` is `.Rbuildignore`d and is rendered by the maintainer alone.** Nothing on CRAN
  executes it, ever. It may therefore evaluate a real scoring run against a staged asset. The
  rendered `README.md` **is not** `.Rbuildignore`d and does ship in the tarball, so what lands in
  it still has to be clean.

| | `vignettes/*.Rmd` | `README.Rmd` |
|---|---|---|
| built by `R CMD check` and CRAN | yes | **no** |
| may reach the network | no | no |
| may read a staged asset | no | **yes** |
| clocks an evaluated chunk may score | bundled only | any |
| ships in the tarball | yes | no, but `README.md` does |

**Rendering `README.Rmd` wants the assets staged**, though not because a chunk scores an external
clock. The scoring example is `Horvath1`, which is bundled. It is the `list_mc_assets()` pair that
needs them, because the point of that pair is `downloaded` moving from `FALSE` to `TRUE`. With an
empty assets directory the render still succeeds and the table simply reads `FALSE` twice, so the
failure is a wrong narrative rather than a broken build. Not worth designing around: the package
tracks exactly one `README.md`, so whoever regenerates it is by definition the maintainer, and a
maintainer has the assets.

### Chunks

- **In a vignette, a chunk that needs the network, writes to disk, deletes, prompts, or prints a
  per-machine path is `eval = FALSE`**, with its output pasted below the call. This is the same
  policy as `@examplesIf interactive()` and it is not negotiable.
- **In a vignette, a chunk that does evaluate must run offline with no asset.** Read the shipped
  catalog, or build inputs with `sim_DNAm()` and bundled clocks only. `Horvath1` and `Hannum` are
  safe, the four external groups are not.
- **In `README.Rmd`, prefer a chunk that evaluates.** Real output cannot rot. Reserve
  `eval = FALSE` for the calls that genuinely must not run during a render, which today is the
  download itself.
- **Pasted output is a claim about behaviour and rots exactly like a `@details` sentence.** Run
  the call, then paste what it actually printed. Never retype it and never adjust it. The
  vignette's hand-written asset sizes were wrong by about 9x for an unknown period and fed a bad
  premise into a design decision (DECISIONS 2026-08-04).
- **Pin anything random.** `set.seed()` at the top of the first chunk that simulates, or the
  rendered file churns on every build. This is the opposite of the rule for `@examples`
  (section 8), where nothing checks the numbers and no seed is wanted.
- **Say what the reader does, not what the package does internally.** These files are where the
  pull toward internal vocabulary is strongest, because the mechanism is interesting. R8 still
  applies: it is an **asset**, never a pack, and there is no "closed set".

### Line breaks: do not hard-wrap the prose

**Do not line-break prose outside a code block. One paragraph is one line in the source, and the
editor soft-wraps it.** A hard-wrapped paragraph produces a diff on every reflow, so a one-word
edit rewrites five lines and the real change disappears into the noise.

**Inside a code block, wrap at 80 columns**, because that text is read as code and is not
reflowed by anything.

This is scoped to `README.Rmd` and `vignettes/*.Rmd`. **`dev/` docs, including this file, stay
hard-wrapped** and are unaffected.

### ASCII

The rule binds what an author types, not what a run prints.

- **Never hand-write a non-ASCII character**, in prose or in a chunk you author.
- **Non-ASCII that cli generated in captured output is fine and stays.** Do not set
  `options(cli.unicode = FALSE)` to launder it. Suppressing it makes the rendered file disagree
  with what the reader's own console shows, which is the one thing pasted output exists to
  demonstrate.
- **Pandoc manufactures non-ASCII from ASCII input if it is allowed to.** Its `smart` extension
  turns a straight apostrophe into a curly one, so hand-written prose arrives non-ASCII through no
  fault of the author. `README.Rmd` disables it in the YAML:

  ```
  output:
    github_document:
      md_extensions: -smart
  ```

### Rendering

- Vignette: `knitr::knit("vignettes/x.Rmd", output = tempfile())` after `devtools::load_all()`.
  It catches a chunk that should have been `eval = FALSE` by failing.
- README: `devtools::build_readme(".")`, which installs the package into a temporary library and
  compiles `src/`. **`Rscript --vanilla` cannot find pandoc**, so set `RSTUDIO_PANDOC` first. On a
  machine with RStudio installed, `<RStudio>/resources/app/bin/quarto/bin/tools` holds it.

Current files: `vignettes/assets.Rmd` and `README.Rmd`.
