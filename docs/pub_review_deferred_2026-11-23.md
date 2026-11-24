# Deferred-item pass: 13-mmrm-vs-survival-ad
*2026-11-23 12:49 PST*

Source log: `docs/pub_review_remediation_2026-08-16.md` (the only
`pub_review_remediation_*` file in the workspace), answering
`docs/pub_review_whitepaper_2026-08-16.md`. The later
`analysis/report/referee-remediation-log-2026-11-19.md` (a separate
referee cycle that redesigned the simulation) was used as evidence of
what changed, not as inventory. Its working-tree changes, and the
full render of 2026-11-19 (`report.pdf`, staged as
`share/report-2026-11-19-2056-d102506-wip.pdf`), are uncommitted.
Nothing was staged or committed in this pass.

## 1. Inventory

| # | Source log item | Class |
|---|---|---|
| 1 | Checklist 13 / Minor 3: sensitivity on CDR staging thresholds and dropout-severity coupling | DO-NOW (script) + LONG-RUN (full grid) |
| 2 | Full-scale R = 2000 primary rerun under the new MMRM convergence check | RESOLVED-SINCE |
| 3 | Sample-size sweep (`sim_power_vs_n.rds`) rerun | NO-LONGER-APPLIES |
| 4 | `sim_results.rds` regeneration | LONG-RUN (optional; not read by the report) |
| 5 | Framing / retitling (whitepaper section 5) | AUTHOR-DECISION |
| 6 | Verification of three bib entries added in August | RESOLVED-SINCE |
| 7 | New issue in August log: knitr cache not keyed on `sim_study.R` | RESOLVED-SINCE |

Item 7 was listed under "New issues found", not "Deferred", but it
was explicitly left unfixed, so it is included.

Evidence per item:

- **Item 2.** `report.Rmd` now computes every number at render time.
  The rendered `analysis/report/report.tex` (2026-11-19 22:13, later
  than `sim_study.R` and `report.Rmd`, both 2026-11-19 00:29) reports
  Table `tab:summary-table` with 2000 valid fits and convergence
  1.000 for all three analyses, under the `apVar` and `maxIter`
  convergence check in `fit_mmrm()`. [inspected: file timestamps and
  rendered tex; I did not rerun it]
- **Item 3.** `report.Rmd` no longer reads `sim_power_vs_n.rds`
  (no `readRDS` call remains); the sweep is computed in the
  `run-sweep` chunk, whose cache dates from the 2026-11-19 render.
  `analysis/data/README.md` marks the `.rds` stale. [inspected]
- **Item 4.** No file in `analysis/report/` reads
  `sim_results.rds`; the README marks it stale. Regeneration is
  optional housekeeping; see Long runs.
- **Item 6.** The 2026-11-19 log states that
  `royston2013restricted`, `rizopoulos2011dynamic`, and
  `sun1996survival` were re-verified against Crossref. [inspected
  log claim; not re-verified in this pass]
- **Item 7.** `report.Rmd` lines 72-89 set
  `cache.extra = unname(tools::md5sum(sim_path))` globally and
  `autodep = TRUE`. [inspected]

## 2. Done now

- **Item 1, script half.** `sim_trial()` already exposes
  `dropout_coupling` and `breaks` as arguments (added in the
  November redesign), so no change to `sim_study.R` was needed.
  This matters: `sim_study.R`'s md5 is in every chunk's cache key,
  so any edit to it, even a comment, would force a full rerun of
  about 4,500 replications at the next render. I added
  `analysis/scripts/03_sensitivity.R`, a one-at-a-time grid with
  five cells: the primary scenario (coupling 0.5, CDR 0.5/1.0 cut
  at 0.75) as an internal reference, coupling 0 and 1, and cut 0.65
  and 0.85. Each cell recomputes its own large-sample targets with
  `large_sample_targets()` under the same breaks, because the
  stage-MMRM and Cox targets depend on them. Defaults are one core,
  R = 500 per cell, seed 20261123; output goes to
  `analysis/data/sim_sensitivity.rds`. The grid values are
  illustrative and need author confirmation (see Author decisions).
  `[verified]` structurally: a reduced run
  (`Rscript analysis/scripts/03_sensitivity.R 4 2e4 1 <scratch>`)
  completed all five cells in 72 s single-threaded with
  `n_valid = 4` for every method and cell. The output went to the
  scratchpad, not to `analysis/data/`. With R = 4 the numbers carry
  no information and are not reported.
- **Results TODO (new item, not acted on).** `report.Rmd` line 727
  asks the author to confirm the abstract's power ordering against
  the rendered tables. The rendered values confirm it: power 0.998
  (latent MMRM), 0.945 (stage MMRM), 0.845 (Cox) in the primary
  scenario, MCSE at most 0.008, and the same ordering in the
  smaller-effect (0.666, 0.466, 0.304) and MCAR (1.000, 0.934,
  0.838) rows. [inspected: `report.tex`] I left the TODO in place
  so `report.Rmd` and `report.pdf` stay consistent; deleting it
  needs a re-render, which should hit the cache but would put
  hours of cached results at risk if the cache missed.
- **Tests.** `pkgload::load_all('.');
  tinytest::run_test_dir('inst/tinytest')`: 32 results, all OK
  (9.3 s). One assertion calls `run_simulation(n_reps = 2,
  n_cores = 2L)`, so that test briefly forks 2 workers; the rest
  ran single-threaded.

## 3. Long runs

Timing (single core, this machine, 2026-11-23): `run_simulation()`
took 32.5 s for 10 replications at 200 per arm (about 3.25 s per
replication), and `large_sample_targets(2e5)` took 29.5 s.

- **Sensitivity grid (item 1).**
  `Rscript analysis/scripts/03_sensitivity.R 500 2e5 1`. That is
  2,500 replications plus 5 target runs, about 2.3 h on one core.
  On an idle machine, `... 500 2e5 8` should take about 20 to 30
  min (extrapolated, not timed). Afterwards, add a short
  sensitivity table to Results, read inline from
  `sim_sensitivity.rds`. Then revise the sentences that say no
  sensitivity analysis was done: `report.Rmd` about lines 420-427
  (thresholds "a design choice ... no sensitivity analysis"),
  525, 661-663 (ADEMP "partly met"), and 1063-1066 (Discussion).
  Then render with `bash tools/render.sh analysis/report/report.Rmd`.
  If `sim_study.R` is left untouched, the existing chunks should
  restore from cache; back up `analysis/report/cache/` first.
- **`sim_results.rds` regeneration (item 4).**
  `Rscript analysis/scripts/01_run_simulation.R`. The script
  hardcodes `n_cores <- 8L`. It runs 2,000 replications plus a
  200,000-per-arm target run, about 1.8 h of single-core work, so
  roughly 15 to 30 min on 8 free cores. Nothing in the manuscript
  needs updating afterwards. The alternative is to leave the file
  marked stale, as it is now.

## 4. Author decisions

- **Framing and title (item 5).** The whitepaper recommended
  framing (B): a simulation study, retitled to drop "Historical
  Review", with the history sections trimmed. The abstract now
  follows framing (B), narrowing the contribution to separating
  the cost of a coarse instrument from the cost of the
  time-to-event summary. The title is still "MMRM vs. Survival
  Analysis for AD Clinical Trials: A Historical Review and
  Simulation Study". The options are (A) keep the dual
  review-plus-simulation identity, (B) retitle and trim per the
  whitepaper, or (C) split it into a short perspective piece and a
  methods paper. Recommendation: (B). The abstract and the
  three-analysis design already argue for it, and the title is now
  the main thing out of step with them.
- **Sensitivity grid values (item 1).** The coupling values
  {0, 0.5, 1} and CDR 0.5/1.0 cut points {0.65, 0.75, 0.85} in
  `03_sensitivity.R` are my illustrative choices. Narrowing the
  0.5 band also changes the screened population, and the cut
  points 0.25, 1.5, and 2.5 and the residual SD are not varied.
  Recommendation: confirm or widen the grid before the long run.
  One-at-a-time variation is enough to answer Minor issue 3.

## 5. New issues found

- **The full-scale results exist only in the working tree.** The
  2026-11-19 render (`report.pdf`, `report.tex`, staged share PDF)
  and the code that produced it are uncommitted, and
  `analysis/report/cache/` is gitignored. If the working tree or
  cache is lost, regenerating means about 4,500 replications (the
  last render ran for about 77 min on 8 cores). Recommend that the
  author commit the November changes.
- **Driver scripts hardcode `n_cores = 8L`.**
  `01_run_simulation.R` and `02_power_vs_n.R` take no core-count
  argument, so they cannot be run politely on a shared machine.
  Results do not depend on core count (per-replication L'Ecuyer
  streams, covered by a test), so a CLI argument would be safe.
  Not changed: these are reproducibility drivers, and the change
  is cosmetic.
- **Overlap with open November items, for cross-reference.** The
  2026-11-19 log also defers the residual-error and trial-duration
  sensitivity (M3 iii), relative sample size for 80% power (M3
  iv), and a Schoenfeld proportional-hazards check for the Cox
  model (M6). Not acted on here.
