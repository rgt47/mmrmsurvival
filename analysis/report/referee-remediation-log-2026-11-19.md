# Referee remediation log: report.Rmd
*2026-11-19 00:22 PST*

Remediation of `referee-report-2026-11-18.md` for
`analysis/report/report.Rmd`. Nothing was staged, committed,
pushed, or rendered. The voice and humanizer edits of
2026-11-18 were kept wherever the text survived; several
Methods, Results, and Discussion passages were rewritten
because the design changed.

The central change is to the simulation itself. The report
now computes every reported number during the render, so no
number in the manuscript was typed by hand. The full render
(about 4,500 replications) was not run in this session; see
Deferred.

## Fixed

- **M1 (MMRM analyzed the latent variable).**
  `sim_study.R`, `report.Rmd`. `fit_mmrm()` gained an
  `outcome` argument; `run_simulation()` now fits a second
  MMRM on the observed global CDR stage. Methods state that
  the latent-score MMRM is an approximate upper bound and that
  CDR-SB is not simulated. The false sentences ('primary
  continuous outcome is the global CDR score', 'global CDR
  ... as a continuous outcome', 'not a continuous surrogate',
  the CDR-SB density caption) were removed or corrected.
  [verified: reduced runs and tests]
- **M2 (population did not match description).**
  `sim_trial()` now simulates a month-0 screening and
  randomization visit and enrolls only candidates at observed
  CDR 0.5. Results report event fractions, month-6 conversion
  fractions, and completion by arm, all computed inline. The
  Discussion no longer says conversion is modest or censoring
  heavy; it lists the actual sources of Cox power loss.
  [verified on reduced runs; see below]
- **M3 (one scenario, no null).** Added null
  ($\beta_3 = 0$), smaller-effect ($\beta_3 = -0.0375$), and
  MCAR-dropout scenarios at R = 500 each, with a scenario
  table and inline type I error rates. The power ratio was
  dropped; the text now says the latent MMRM is at the
  ceiling. [applied; chunk code executed in a reduced dry
  run, see Checks]
- **M4 (estimand and baseline).** Change is now from month
  0 (pre-randomization). The latent truth is
  $\beta_3 \times 2.5 = -0.1875$, confirmed by a large
  simulated trial (-0.1871 at 200,000 per arm). The MMRM now
  includes a baseline-by-visit term, because the single pooled
  baseline slope misspecified the mean (see New issues). An
  ICH E9(R1) estimand paragraph was added. [verified]
- **M5 (dropout informative for Cox).** New 'Dropout'
  subsection says the dropout is MAR only for the latent
  MMRM and violates independent censoring for Cox (cites
  Kalbfleisch and Prentice 2002). MCAR scenario added.
  [applied]
- **M6 (Cox truth).** `large_sample_targets()` computes the
  large-sample limit of the Cox log hazard ratio (and the
  observed-stage MMRM target) from one 200,000-per-arm trial
  without dropout. The nominal log(0.75) is gone. In a
  scratch run the Cox limit was -0.334 and the stage target
  -0.194. [verified]
- **M7 (review errors).** (a) Study 201 now described as a
  Bayesian ADCOMS analysis at 12 months (Swanson 2021);
  Clarity AD (van Dyck 2023) is the CDR-SB MMRM example;
  Dhadda 2022 described as a consistency analysis.
  (b) 'nearly every' replaced by 'many'; TRAILBLAZER-ALZ 2
  added. (c) PMRM attributed to Raket 2022; Jönsson 2024 as
  application. The unverified sentence 'Its efficiency is
  comparable to that of MMRM' was removed. (d) New paragraph
  on Sano 1997, Petersen 2005, and Thomas et al. 2000.
  (e) Placebo-dropout claim rewritten as conditional, with
  the cholinesterase-inhibitor counterexample. (f) 'OS is
  binary' corrected to death. [applied]
- **M8 (reproducibility).** Global
  `cache.extra = tools::md5sum(sim_study.R)`; echoed code now
  comes from `sim_study.R` via `knitr::read_chunk()` (labels
  confirmed with `knitr:::knit_code`); the stale sweep `.rds`
  is no longer read, the sweep is computed in the report.
  Each simulation chunk sets its own seed, and each
  replication has its own L'Ecuyer stream, so results do not
  depend on core count. `autodep = TRUE` added so edits to a
  simulation chunk invalidate the chunks that use it.
  `run_simulation()` stops if any parallel replication is
  lost (a killed `mclapply` child returns NULL).
  The 'per-replication RNG states retained' claim and all
  changelog language were removed. `01_run_simulation.R`,
  `02_power_vs_n.R`, and `analysis/data/README.md` updated to
  match. [verified: purl, read_chunk, and a tinytest that two
  calls with one seed match across core counts]
- **Section 6 references.** Added and checked against
  Crossref (DOI, volume, issue, first page, authors):
  `swanson2021lecanemab`, `vandyck2023lecanemab`,
  `sims2023donanemab`, `raket2022pmrm`,
  `sano1997selegiline`, `petersen2005vitamine`,
  `thomas2000longitudinal`, `altman2006cost`,
  `kalbfleisch2002statistical`. Also re-verified
  `jonsson2024progression`, `dhadda2022lecanemab`,
  `li2019relative`, `donohue2011relative`, and the three
  August entries (`royston2013restricted`,
  `rizopoulos2011dynamic`, `sun1996survival`); all match
  Crossref. Li et al. now described as a simulation from an
  ADNI-fitted model, not observational data. [verified]
- **m1** changelog language removed. **m2** 'not an accident
  of history' replaced by a hedged view that precedent
  matters at least as much as biology. **m3** Leber global
  measure described as a clinician's impression. **m4**
  realized completion reported inline; base-rate figure
  (about 74%) derived from the hazard. **m5** Wald z
  reference and absence of KR/Satterthwaite df stated. **m6**
  sweep caption no longer asserts dominance; Wilson
  intervals via `wilson_ci()`. **m7** density caption fixed.
  **m8** Future research item 1 renamed. **m11** DPM
  sensitivity sentence removed. **m12** opener already gone.
  [applied]
- **Tests.** `inst/tinytest/test_basic.R` rewritten for the
  new design: 32 assertions, all passing (see Checks).

## Deferred

- **Full render.** Not run (instructions forbid render, and
  the runs exceed 10 minutes). Command:
  `bash tools/render.sh analysis/report/report.Rmd`. Under the
  load seen in this session (load average 65 to 93 on 8
  cores), 144 replications took 4.9 minutes on 8 cores, so the
  render may take several hours on a loaded machine. After it,
  check the abstract's qualitative results sentence and the
  Discussion's 'conversion was common and censoring light'
  against the rendered numbers.
- **M3 (iii) and (iv).** Residual error, staging thresholds,
  dropout coupling, and trial duration are not varied;
  relative sample size for 80% power is not reported. Both
  disclosed in the text. Needs a scenario grid; budget.
- **M6 PH check.** No Schoenfeld-residual assessment.
- **Section 7 suggestions.** Simulated CDR-SB, ordinal or
  interval-censored models, RMST, PMRM comparator, and
  calibration to ADNI not done. CDR-SB would need invented
  rating-error values; left for the author.
- **Unverifiable items.** The 2023 Raket/Jönsson conference
  abstract (not cited; full citation unknown). Ito and
  Hutmacher's '1-point change in about 2 years' figure (m9):
  the number was removed, the citation kept. The exact Li et
  al. power ratio (m9): kept as 'roughly doubled'.
  `leber1990guidelines`, `fda2018earlyad`, and
  `fda2018oncology` still lack identifiers or URLs (m10); not
  verified, not invented.
- **Facts for author confirmation** (from knowledge, not
  searched): Sano 1997 primary endpoint (time to composite
  milestones); Petersen 2005 primary endpoint (time to AD
  diagnosis); EMERGE/ENGAGE and Clarity AD primary analyses
  (CDR-SB, MMRM); TRAILBLAZER-ALZ 2 primary model (natural
  cubic spline on iADRS) and its reported risk of progression
  to the next CDR-G stage; the cholinesterase-inhibitor
  adverse-event dropout pattern; the closing judgment that
  the field moved away from milestone endpoints in part for
  efficiency reasons (hedged as 'perhaps' and 'likely').
- **Stale data files.** `sim_results.rds`,
  `sim_power_vs_n.rds`, and `sim_13.rds` were not deleted or
  regenerated; README marks them stale.

## Illustrative values chosen in this remediation (confirm)

- Smaller effect $\beta_3 = -0.0375$ (12.5% slowing),
  labeled illustrative in Methods.
- MCAR dropout at a constant 15% per year.
- Entry rule: observed global CDR 0.5 on one month-0
  assessment; months 0 and 6 always observed.
- Replication counts: R = 500 per secondary scenario,
  R = 200 per sweep cell; large-sample targets at 200,000
  per arm; seeds 20260309, 20260312, 20260313 (new) plus
  20260310 and 20260511 (existing); `n_cores = 8`.

## New issues found while fixing

- **Pooled baseline slope (suspected).** Under the new
  month-0 design, a 48-replication scratch run of
  `change ~ trt * visit_f + baseline` gave latent-MMRM bias
  -0.011 (about 1.8 MCSE) and coverage 0.875. In principle
  the regression of change on baseline differs by visit
  under random slopes, so a pooled term misspecifies the
  mean, which can bias the contrast under outcome-dependent
  dropout. Changed to `baseline * visit_f`; an 80-replication
  run then gave bias 0.001 (MCSE 0.005). The two runs used
  different random streams, so this is suggestive, not
  shown. The same concern applies in principle to the old
  design (referee M4).
- **Cache and seeds.** knitr caches `.Random.seed` with a
  chunk, so the old single-seed setup was probably not
  broken. Per-chunk seeds were added as a robustness choice,
  so that each simulation can be reproduced by the driver
  scripts in isolation.
- **The abstract misnamed a source.** It cited 'Donohue and
  Sabbagh (2011)'; the bib entry is Donohue, Gamst, Thomas
  et al. Corrected to 'Donohue et al. (2011)'.

## Reduced-run evidence (scratch, not in the manuscript)

Final code, 80 replications per scenario, 200 per arm,
8 cores (6.0 minutes). Targets used: -0.1875 (latent),
-0.1939 (stage), -0.3340 (Cox), from a 200,000-per-arm
scratch run. These numbers are not in the manuscript; the
render recomputes everything at full R.

| Scenario | Method | Bias (MCSE) | Power | Coverage |
|---|---|---|---|---|
| Primary | MMRM latent | 0.001 (0.005) | 1.000 | 0.912 |
| Primary | MMRM stage | 0.004 (0.007) | 0.887 | 0.925 |
| Primary | Cox | 0.023 (0.013) | 0.750 | 0.963 |
| Null | MMRM latent | 0.005 (0.004) | 0.025 | 0.975 |
| Null | MMRM stage | 0.001 (0.006) | 0.088 | 0.912 |
| Null | Cox | 0.008 (0.011) | 0.025 | 0.975 |

Primary-scenario descriptors (mean over replications):
events 83% placebo / 74% active; month-6 conversion 31% /
25%; completion 70% / 71%. Smaller-effect and MCAR paths
were exercised only by the reduced dry run of all report
chunks (see Checks); an earlier run of those scenarios used
superseded code and is not reported.
The stage MMRM's null rejection rate (7/80) and its
coverage are worth watching in the full run; at R = 80 they
are within Monte Carlo error of nominal but not tightly.

## Checks

- tinytest (`pkgload::load_all("."); run_test_dir(...)`):
  32 results, all OK.
- `knitr::purl()` on `report.Rmd`: succeeds (after making the
  `read_chunk()` path literal, since purl evaluates that call
  alone).
- Every `@key` in `report.Rmd` exists in `references.bib`
  (not a symlink; 35 keys cited).
- Em dashes in `report.Rmd`: 0 before, 0 after.
- YAML front matter parses with `yaml::yaml.load()`.
- Reduced dry run of every report chunk: the purled report,
  with R cut to 24 (primary), 16 (each scenario), 8 (sweep
  cells at n = 100 and 200), and targets at 20,000 per arm,
  ran without error, and all 52 inline expressions
  evaluated. This is a structural check of the chunk code,
  not evidence about the full-scale numbers.
- Sentence-length scan of body prose: remaining hits over
  25 words are scanner artifacts (bullet lists joined into
  one sentence, sentences that begin with a citation key).
- `git status`: only the seven intended files changed, plus
  this log.

## Abstract voice pass

Old: 231 words, 7 sentences (longest 53). New: 226
words, 16 sentences (longest 23 by a whitespace count; the
splitter breaks at 'et al.', so the count is approximate).
No structured headings were present or required; no word
limit is named, so length was held at the old count. YAML
block scalar and two-space indentation kept;
`yaml::yaml.load()` parses the header.

Factual changes:

1. 'Donohue and Sabbagh (2011)' corrected to 'Donohue et
   al. (2011)' (wrong author in the old abstract).
2. 'in observational ADNI data' removed (Li et al. 2019 is a
   simulation; M7/section 6).
3. 'staged categorical instrument (CDR global stage, not a
   continuous surrogate)' removed: the old MMRM analyzed the
   latent score (M1). Now three analyses are named, and the
   latent MMRM is called an approximate upper bound.
4. 'realistic, severity-dependent missing-at-random dropout'
   replaced by 'dropout rises with latent severity' and a
   statement that it is informative for the Cox model (M5).
5. Entry at CDR 0.5 stated (M2; now true by screening).
6. Secondary scenarios named (M3).
7. 'Across 2,000 simulated trials ... MMRM substantially
   outperformed' replaced by 'In the primary scenario, power
   was highest for the latent-score MMRM and lowest for the
   Cox model.' This was checked only on an 80-replication
   run; a visible TODO at the top of Results asks the author
   to confirm it after rendering.
8. The novelty claim was narrowed to the decomposition of
   the loss (coarse instrument versus time-to-event
   summary); 'a combination prior work has not simulated in
   this precise form' removed.
9. Opening 'converged on divergent primary analyses of
   arguably comparable disease processes' simplified to
   'analyze progressive disease in different ways'.

No sentence exceeds 25 words.
