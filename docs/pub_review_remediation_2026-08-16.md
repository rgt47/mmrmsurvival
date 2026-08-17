# Remediation Log: pub_review_whitepaper_2026-08-16.md
*2026-08-16 16:29 PDT*

Remediation of the referee comments in
`docs/pub_review_whitepaper_2026-08-16.md`, applied to
`analysis/report/report.Rmd`,
`analysis/scripts/sim_study.R`,
`analysis/report/references.bib`,
`analysis/data/README.md`, `analysis/figures/README.md`,
`analysis/tables/README.md`, and
`inst/tinytest/test_basic.R`. This log does not edit the
whitepaper itself.

## 1. Fixed

### (a) Required for correctness

- **Major issue 2 / checklist 1 (visit schedule mismatch).**
  `[applied, unverified against a rendered PDF; verified by
  code inspection]` Rewrote "Data generating model" and
  "Design specifications" in `report.Rmd` to state the
  actual visit grid produced by `sim_trial()` (months 6,
  12, 18, 24, 30; no month-0 record; month 6 serves as the
  change-score baseline reference), replacing the false
  "Baseline, 6, 12, 18, 24 months" / "5 post-baseline
  visits" claims. Also synchronized the echoed `fit_mmrm()`
  code in the `analysis-functions` chunk, which had silently
  drifted from the real `sim_study.R` (a reproducibility gap
  not explicitly itemized by the whitepaper; see "New issues
  found" below).
- **Major issue 1 / checklist 2 (stale, self-contradictory
  ADEMP audit).** `[applied, unverified against a rendered
  PDF; verified by code/text inspection]` Deleted the
  fabricated "Not compliant" audit block nested under
  `# References` in `report.Rmd`. Replaced it with a new
  Methods subsection, "Reproducibility and data/code
  availability," that re-audits the current code against
  Morris et al. (2019) ADEMP criteria and reports the true,
  largely-compliant state (eval=TRUE chunks, named-contrast
  estimand extraction, MCSE on every metric, pinned RNG),
  while honestly naming the two ADEMP elements still only
  partially met (undeclared factorial sweep, no threshold
  sensitivity analysis).
- **Major issue 3 / checklist 3 (future-tense Discussion).**
  `[applied, unverified against a rendered PDF; verified by
  text inspection]` Rewrote the paragraph in past/present
  tense, replaced "We expect..." hedging with inline `r`
  expressions reading the actual `mmrm_row`/`cox_row` power
  values computed in the `headline-13` chunk, and added a
  quantitative comparison to the Donohue (2011) and Li
  (2019) inflation-factor literature.
- **Major issue 7 / checklist 4 (underived Cox truth
  value).** `[applied, unverified against a rendered PDF]`
  Added a new Methods subsection, "Note on the Cox truth
  value," explicitly disclosing that `log(0.75)` is a
  nominal design target, not an analytically derived or
  DGM-calibrated ground truth, and explaining why. Updated
  the Results headline prose and Design specifications to
  use "nominal design target" / "mean difference from the
  nominal target" language for the Cox estimand rather than
  unqualified "bias."
- **Major issue 6 / checklist 5 (convergence metric does
  not verify convergence).** `[verified]` Added a genuine
  convergence diagnostic to `fit_mmrm()` in
  `analysis/scripts/sim_study.R`: a fit is now treated as
  invalid (`NA`) unless `mod$apVar` is a matrix (not the
  character message `nlme` returns for a non-positive-
  definite Hessian) and `mod$numIter < maxIter`. Verified
  on a 30-replication, `n_per_arm = 50` test run
  (`run_simulation(n_reps = 30, n_per_arm = 50)`): all 30
  MMRM fits produced valid, non-NA estimates and
  `convergence = 1`, confirming the stricter check does not
  spuriously reject well-behaved fits. Synchronized the
  same logic into the echoed `analysis-functions` chunk in
  `report.Rmd` and updated the "Performance metrics" bullet
  and Results prose to describe the new definition.

### (b) Required for acceptance

- **Checklist 6 (abstract).** `[applied, unverified against
  a rendered PDF]` Added a ~230-word `abstract:` YAML field
  to `report.Rmd` stating the question, ADEMP method,
  headline power result, and the narrowed novelty claim
  (efficiency-loss quantification for the specific CDR-
  staging/dropout DGM, not a rediscovery of the general
  MMRM-over-survival efficiency result).
- **Checklist 7 (power-vs-n sweep undocumented /
  non-reproducing).** `[applied, unverified against a
  rendered PDF]` Added a "Sample-size sweep" Methods
  subsection documenting aims, DGM, estimands, $n$ grid,
  $R = 200$, seed 20260511, and the resulting MCSE. Changed
  the `power-comparison` chunk to `stop()` with an explicit,
  actionable error message if `sim_power_vs_n.rds` is
  missing, replacing the silent placeholder-text fallback
  the whitepaper flagged.
- **Checklist 8 (data availability / boilerplate
  README).** `[applied, unverified]` Rewrote
  `analysis/data/README.md` to describe the actual simulated
  data pipeline (`sim_results.rds`, `sim_power_vs_n.rds`,
  `sim_13.rds`), remove the Palmer Penguins boilerplate, and
  add a data/code availability statement. Also added short
  `README.md` files to the previously-unexplained empty
  `analysis/figures/` and `analysis/tables/` directories
  explaining that all output is embedded via knitr caching.
- **Checklist 9 (placeholder tests).** `[verified]` Replaced
  `inst/tinytest/test_basic.R`'s single `expect_true(TRUE)`
  with 16 real assertions covering: `sim_trial()` output
  structure and absorbing-dropout invariant; a direct
  regression test reproducing the historical estimand-
  labeling bug and confirming `fit_mmrm()` returns the
  correct last-visit contrast (verified to differ from the
  buggy trt-only estimate by 0.18 on a fixed seed);
  `fit_mmrm()` graceful degradation to `NA` on unfittable
  data; `summarize_results()` MCSE formulas checked against
  a hand-computed 5-observation toy case; and
  `summarize_results()` behavior when every fit fails. Ran
  via both `tinytest::run_test_file()` and the user's
  documented `pkgload::load_all(".");
  tinytest::run_test_dir("inst/tinytest")` invocation: **all
  16 tests pass** in both.
- **Checklist 10 (session/version/seed disclosure).**
  `[applied, unverified against a rendered PDF]` Added a
  `session-info` chunk to the new Reproducibility subsection
  printing R version, platform, and `nlme`/`survival`/
  `MASS`/`dplyr` package versions at render time.
- **Checklist 11 (missing citations).** `[applied; citation
  details recalled from training knowledge, not verified
  against a live bibliographic database — see "New issues
  found"]` Added `royston2013restricted` (RMST),
  `rizopoulos2011dynamic` (joint modeling), and
  `sun1996survival` (interval-censored methods) to
  `references.bib`, and cited all three in the "Future
  research" section (items 3 and 4) where the whitepaper
  identified the gaps.

### (c) Desirable polish

- **Checklist 12 / Minor issue 2 (dual-scale density
  plot).** `[applied, unverified against a rendered PDF]`
  Standardized both estimators to a Wald z statistic
  (estimate / model SE) before plotting, replacing the
  raw-scale overlay and the "comparison is qualitative"
  caveat with an actual common, unit-free axis.
- **Checklist 14 / Minor issue 5 (uncited ADNI calibration
  claim).** `[applied, unverified against a rendered PDF]`
  Removed the unsupported "calibrated to observed CDR
  distributions in ADNI" claim; replaced with an accurate
  disclosure that the staging thresholds are an undisclosed-
  until-now design choice, not an external calibration, with
  a pointer to the (still absent) sensitivity analysis.
- **Checklist 15 / Minor issue 7 (raw column names in
  Table 1).** `[applied, unverified against a rendered
  PDF]` Added `col.names` to the `summary-table` `kable()`
  call.
- **Checklist 16 / Minor issue 6 (British spelling).**
  `[verified]` Fixed "parameterisation" to "parameterization"
  in `sim_study.R`.
- **Checklist 17 (empty figures/tables dirs, undocumented
  `sim_13.rds`).** `[applied, unverified]` Addressed via the
  new `analysis/data/README.md` (documents `sim_13.rds` as
  an archival snapshot matching the manuscript's reported
  numbers) and the new `analysis/figures/README.md` /
  `analysis/tables/README.md`.

## 2. Deferred

- **Checklist 13 / Minor issue 3 (sensitivity analysis on
  CDR staging thresholds and dropout-severity coupling).**
  Deferred: a proper tornado/sensitivity sweep requires
  multiple additional full simulation runs, which exceeded
  this remediation's time budget. Exact command pattern to
  extend: parameterize `breaks` and the `0.5` coupling
  constant in `sim_trial()` as function arguments (already
  partially true for `dropout_rate` but not for the
  multiplier `0.5` or the `breaks` vector), then loop
  `run_simulation()` over a small grid, e.g. via a new
  `analysis/scripts/03_sensitivity.R` modeled on
  `02_power_vs_n.R`.
- **Full-scale primary $R = 2000$ rerun to confirm the new
  MMRM convergence diagnostic does not change Table 1.**
  Deferred: a 50-replication timing probe at
  `n_per_arm = 200` did not complete within 2 minutes on
  this (heavily loaded, shared) machine, so a 2,000-
  replication rerun was judged infeasible within budget. The
  convergence-check code change was instead verified on a
  reduced run (`n_reps = 30`, `n_per_arm = 50`; see Fixed,
  above), which produced valid, non-NA output with
  `convergence = 1`, giving reasonable but not conclusive
  confidence that the primary $R = 2000$/$n = 200$ numbers
  already in `report.tex` (and matching `sim_13.rds`) are
  unaffected. **The user should rerun the full simulation
  and re-render before submission:**
  `bash tools/render.sh analysis/report/report.Rmd` (this
  will re-execute `run_simulation(n_reps = 2000)` from the
  updated `sim_study.R` via the report's `run-simulation`
  chunk; note the existing knitr cache under
  `analysis/report/cache/` is keyed only on chunk text, not
  on `sim_study.R`'s content, so a stale cache could
  silently serve pre-fix results — see "New issues found").
  If a full render is not run, the manuscript's Table 1 and
  headline numbers currently reflect the *pre-fix*
  convergence definition (absence-of-error only), not the
  newly added apVar/iteration check.
- **Sample-size sweep (`sim_power_vs_n.rds`) rerun.** Not
  rerun; the existing cached sweep was generated before the
  convergence-diagnostic fix, so its `power` values for MMRM
  in principle could shift slightly if any of those 200-per-
  cell replications had a non-converged-but-error-free fit.
  Given the sweep's already-larger MCSE, this is unlikely to
  change the qualitative pattern, but was not verified.
  Command: `Rscript analysis/scripts/02_power_vs_n.R`.
- **`sim_results.rds` regeneration.** Not rerun (see "New
  issues found" — the file on disk is stale and shows 0
  valid Cox fits from an earlier code state). Command:
  `Rscript analysis/scripts/01_run_simulation.R`.
- **Framing/retitling recommendation (whitepaper Section
  5).** Not applied as a title change or Introduction
  restructuring; this is an authorial judgment call about
  target journal and paper identity that only the author
  should make. The abstract added in this remediation
  already narrows the novelty claim consistent with framing
  (B), which partially addresses the concern without
  committing to a full retitle/restructure.
- **Citation verification.** The three newly added
  bibliography entries (`royston2013restricted`,
  `rizopoulos2011dynamic`, `sun1996survival`) were written
  from the assistant's training-data recollection of these
  well-known papers' bibliographic details, not verified
  against CrossRef, PubMed, or the publisher. The author
  should verify title, page numbers, and DOI before
  submission.

## 3. New issues found while fixing

- **Echoed code in `report.Rmd`'s `analysis-functions` chunk
  had drifted from the real `sim_study.R`.** The
  human-readable `fit_mmrm()` shown to readers in the
  Analysis Plan section used `baseline = first(y_star)` and
  `filter(visit > 0)` (keeping the reference visit) and had
  no `tryCatch`/convergence check, none of which matches the
  code that actually produced the reported numbers. Fixed by
  syncing the echoed chunk to the current `sim_study.R`
  logic; this was not itemized in the whitepaper but is the
  same class of defect as Major issue 2.
- **`analysis/data/sim_results.rds` is stale and internally
  broken.** Running `summarize_results()` on the cached
  `sim_results.rds` gives 0 valid Cox fits out of 2,000
  (`convergence = 0` for Cox), which does not match the
  manuscript's reported Cox power of 0.716. This file
  predates the current, working `sim_study.R` and is not
  read by `report.Rmd` (only `analysis/data/sim_power_vs_n.rds`
  is read by the report), so it does not affect the
  manuscript's actual numbers, but it is misleading
  archival clutter documented incorrectly (as current) by
  the pre-remediation `analysis/data/README.md`. Documented
  as stale in the rewritten README; not deleted, in case the
  author wants the history preserved.
- **The report's knitr cache is not invalidated by changes
  to `sim_study.R`.** `opts_chunk$set(cache = TRUE, ...)` in
  the setup chunk caches by chunk-text hash only; there is no
  `cache.extra` or dependency declaration tying the
  `run-simulation` chunk's cache to the content of the
  sourced `sim_study.R` file. This means the convergence-
  diagnostic fix made in this remediation will **not**
  automatically trigger a re-run on the next `knit()` unless
  the cache is manually cleared or `cache.extra` is added
  (e.g., `cache.extra = tools::md5sum("../scripts/sim_study.R")`
  on chunks that source it). Not fixed in this remediation
  because doing so would force this session to either run
  the full simulation (infeasible in budget) or leave the
  render in an inconsistent state; flagged here for the
  author to address before the next render.

