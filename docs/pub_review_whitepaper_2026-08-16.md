# Referee Review: MMRM vs. Survival Analysis for AD Clinical Trials

*Review date: 2026-08-16 16:03 PDT*

Manuscript under review: `analysis/report/report.Rmd` (rendered
`analysis/report/report.pdf`, staged copy
`analysis/report/share/report-2026-08-15-1638-4b5c7a8.pdf`),
"MMRM vs. Survival Analysis for AD Clinical Trials: A Historical
Review and Simulation Study," R. G. Thomas. Single manuscript in
this repository; no companion reports were found under `analysis/`,
`docs/`, `paper/`, `manuscript/`, or `vignettes/`.

## 1. Summary of the work under review

The manuscript combines a narrative history of analytic conventions
in Alzheimer's disease (AD) and oncology trials with a Monte Carlo
simulation study. The historical section (Introduction, roughly
half the document) traces the regulatory origin of the FDA
dual-outcome requirement for AD trials (Leber 1990), the empirical
case against last-observation-carried-forward (LOCF) imputation
(Mallinckrodt and coauthors, the 2010 NRC panel report), the
resulting dominance of the mixed model for repeated measures (MMRM)
in AD trials, and the contrasting persistence of time-to-event
survival analysis in oncology. It then poses the question motivating
the simulation: given that the Clinical Dementia Rating (CDR) scale
supports both a change-score and a time-to-event operationalization,
what is the relative statistical efficiency of MMRM versus a Cox
model applied to the same underlying disease process? The Methods
section specifies an ADEMP-structured (Morris, White, and Crowther
2019) simulation: a two-arm trial, $n=200$ per arm, a bivariate
random-intercept/random-slope latent CDR trajectory discretized into
observed CDR stages, MAR dropout, an MMRM on change from baseline in
latent CDR, and a Cox model on time to first observed CDR $\geq 1.0$.
The Results section reports point estimates for bias, empirical SE,
model SE, power, coverage, and convergence for both methods across
2,000 replications (verified: `analysis/report/report.tex` lines
744-757 show MMRM power 0.994 vs. Cox power 0.716, both with
reported Monte Carlo SEs), plus two figures: power versus per-arm
sample size (from a separately cached 200-replicate-per-cell sweep)
and a density plot overlaying the two estimators' sampling
distributions. The Discussion interprets the findings as consistent
with the historical literature favoring continuous-outcome
efficiency, and lists caveats (asymmetric information use, visit-grid
censoring of the survival endpoint, threshold-dependence, and
distinct clinical interpretations of the two estimands). A Future
Research section lists seven extensions. The manuscript closes with
a self-audit against the Morris et al. (2019) ADEMP reporting
checklist, embedded oddly under the "References" heading, which
declares the simulation "Not compliant" and lists gaps that,
on inspection of the current code and rendered output, appear to
have already been remediated.

## 2. Major issues

1. **Stale, internally contradictory self-audit left in the
   manuscript.** *Location:* `report.Rmd` lines 751-770 (section
   "Morris et al. (2019) ADEMP Compliance," nested inside the
   `# References` heading), referencing
   `docs/morris-audit-2026-04-17.md`. *Problem:* the embedded audit
   (dated 2026-04-17) declares the simulation "Not compliant" and
   lists as unresolved: (a) report chunks are `eval=FALSE` so
   "results are never produced by the rendered report"; (b) the
   estimand extractor at `sim_study.R:92` falls back to
   `grep('trt')[1]`, risking silent mislabeling; (c) "No Monte Carlo
   SE on any metric." *Inspected*, all three are false as of the
   current `report.Rmd` and `sim_study.R`: the simulation chunks
   (`run-simulation`, `summary-table`, both figure chunks) are
   `eval=TRUE` and produce the numbers in Table 1 (*verified* via
   `report.tex`); `sim_study.R` lines 91-119 use an explicit named
   linear combination (`L["trt"]`, `L["trt:visit_f<last>"]`) with no
   grep fallback, and the accompanying code comment even documents
   the historical bug and its fix; and the summary table carries
   `mcse_bias`, `mcse_emp_se`, `mcse_power`, `mcse_coverage` columns
   populated with real values. *Why it matters:* a referee reading
   the manuscript's own closing self-assessment will conclude the
   simulation is not trustworthy, directly undercutting the
   Results section three pages earlier — this is a correctness and
   presentation defect that actively damages the paper's own claims.
   Leaving a superseded audit unedited in a submitted manuscript is
   the kind of internal inconsistency an referee would flag.
   *Remediation:* delete the audit block from the manuscript body
   (it is process documentation, not a manuscript section), or
   replace it with a brief, current, and honest paragraph in
   Methods/Supplement stating what ADEMP elements are met and citing
   an up-to-date audit file; do not nest it under "References."

2. **Visit schedule and "baseline" as described in Methods do not
   match the data-generating code.** *Location:* `report.Rmd` Methods,
   "Data generating model" (lines 305-313, "5 post-baseline visits");
   "Design specifications" (lines 372-380, "Visit schedule: Baseline,
   6, 12, 18, 24 months" and "MMRM endpoint: Change from baseline...
   at 24 months"); code in `analysis/scripts/sim_study.R` lines 6-22
   (`times <- seq(0.5, 2.5, by = 0.5)`, i.e., 6, 12, 18, 24, and 30
   months, with no $t=0$ record generated at all) and lines 64-73
   (`baseline = y_star[visit == min(visit)]`, i.e., the 6-month value
   is treated as "baseline," and the record at `visit == 1` — the
   6-month value itself — is then dropped by `filter(visit > 1)`
   before fitting). *Inspected, verified by code reading.* The text
   asserts a baseline (month 0) assessment and a final visit at
   month 24; the code has no month-0 record, treats the month-6
   value as baseline, and its last visit is month 30, not month 24.
   The internal arithmetic is self-consistent (the "24-month" change
   score is actually the change between the 6-month and 30-month
   observations, a 24-month span, which is why `true_mmrm = -0.075 *
   2` is numerically correct relative to the code's own definitions),
   but the prose describing the design to the reader is wrong on its
   face and cannot be reconciled with the code without working
   through the estimand algebra. Likewise, "5 post-baseline visits"
   over "24 months every 6 months" is arithmetically inconsistent on
   its own terms (that schedule has 4 post-baseline visits, not 5).
   *Why it matters:* reproducibility and correctness — a reader or
   independent reimplementer using the stated visit schedule would
   not replicate the reported estimand or truth value. This is
   exactly the kind of DGM-documentation gap ADEMP review is meant to
   catch. *Remediation:* either (a) add a true $t=0$ baseline visit
   to `sim_trial()` and align the last visit to 24 months, updating
   `true_mmrm` accordingly, or (b) correct the prose to state the
   actual schedule (6, 12, 18, 24, 30 months, first visit serving as
   the baseline reference for the change score) and drop the
   "5 post-baseline visits over 24 months" phrasing.

3. **Discussion is written in the future/hypothetical tense and does
   not engage with the actual computed results.** *Location:*
   `report.Rmd` lines 672-684 ("We expect the simulation results to
   confirm the theoretical predictions from the literature: MMRM
   should yield substantially greater power..."). *Inspected.* The
   Results section immediately above (lines 545-565) already reports
   the realized numbers (power 0.994 vs. 0.716, etc.), so this
   paragraph is a leftover from an earlier draft state in which the
   simulation had not yet been run (consistent with the April audit
   showing `eval=FALSE` chunks at that time) and was never updated
   after the simulation was activated. *Why it matters:* a referee
   would read this as evidence the manuscript was not proofread after
   a substantive revision; it also wastes an opportunity to interpret
   the actual magnitude of the effect (power gap of 0.278, roughly
   40% relative reduction in power for Cox) against the cited
   analytic-inflation-factor literature (Donohue and Sabbagh 2011; Li
   et al. 2019), which would strengthen the paper's contribution.
   *Remediation:* rewrite in past/present tense, quote the actual
   power/bias/coverage numbers, and connect them quantitatively to
   the inflation-factor predictions from the cited prior work (e.g.,
   does the observed power ratio match Donohue's derived inflation
   factor under comparable assumptions? If not, why not?).

4. **No abstract.** *Location:* entire document; the YAML header
   (lines 1-29) has no `abstract:` field and the body proceeds
   directly from the title to "# Introduction." *Inspected.*
   *Why it matters:* every venue in scope (JASA, Biometrics, JCGS,
   Statistics in Medicine — whose CSL is used here) requires a
   structured or unstructured abstract; its absence alone would
   generate a desk-reject or an editorial request before review even
   begins. *Remediation:* add a 150-250 word abstract stating the
   question, method (ADEMP simulation), headline quantitative result,
   and the paper's actual claim (historical synthesis + relative
   efficiency demonstration).

5. **The one figure with a substantive quantitative claim (power vs.
   sample size) is generated by a simulation that is neither described
   in the Methods section nor executed within the report's own build.**
   *Location:* `report.Rmd` "Figures" (lines 588-621, chunk
   `power-comparison`) reads a cached `analysis/data/sim_power_vs_n.rds`
   produced by `analysis/scripts/02_power_vs_n.R`; the chunk falls
   back to placeholder text ("*Awaiting `sim_power_vs_n.rds`...*") if
   that file is absent, i.e., the figure is not guaranteed to
   regenerate from a clean `knit()`. *Inspected.* `02_power_vs_n.R`
   uses its own seed (`set.seed(20260511)`, distinct from the
   manuscript's primary seed 20260310), sweeps $n \in
   \{100,150,200,300,400\}$ at only 200 replications per cell
   (Monte Carlo SE for power near 0.5 is $\sqrt{0.25/200}\approx
   0.035$, three times the primary analysis's precision), and is
   never mentioned in the ADEMP Methods block (which documents only
   the $n=200$, $R=2000$ primary design). The figure caption asserts
   "MMRM uniformly dominates Cox across the explored range" — a
   claim resting on a simulation whose design, replication count, and
   precision are undisclosed in the manuscript text. *Why it matters:*
   reproducibility and completeness — a referee cannot audit the
   design behind a headline figure from the manuscript alone, and the
   report is not self-contained (a `bookdown::render()` from a clean
   checkout without first manually running `02_power_vs_n.R` silently
   downgrades a results figure to a placeholder sentence, which would
   likely go unnoticed at submission). *Remediation:* add an ADEMP-style
   description of the sample-size sweep (aims, $n$ grid, $R=200$,
   seed, and the resulting MCSE) to Methods, increase replications to
   match or justify the discrepancy, and either integrate the sweep
   into the main render (with caching) or make its absence fail loudly
   rather than silently degrading to placeholder text.

6. **"Convergence rate" is defined as absence of a thrown R error, not
   as verified optimizer convergence.** *Location:*
   `analysis/scripts/sim_study.R` lines 75-89 (`fit_mmrm`, `tryCatch`
   around `gls()`, returns `NA` only `error = function(e) NULL`) and
   the "Performance metrics" definition in `report.Rmd` lines 520-521
   ("Convergence rate: Proportion of replications producing valid
   estimates"). *Inspected.* `gls()` from `nlme` does not always throw
   an R error on failure to reach an optimum; it can return a fitted
   object after exhausting `maxIter` or reaching a boundary solution
   without signaling non-convergence via a caught exception. The code
   never inspects `mod$apVar`, iteration counts, or optimizer
   convergence codes. The manuscript reports 100% convergence for both
   methods across 2,000 replications (Table 1), which is plausible for
   this well-identified low-dimensional design but is not actually
   demonstrated by the reported metric — a referee versed in `gls`
   behavior would ask whether "100% convergence" reflects 2,000
   genuinely converged fits or simply 2,000 fits that did not error.
   *Why it matters:* under Morris et al. (2019), convergence/failure
   diagnostics are a required performance measure precisely because
   silent non-convergence biases downstream operating characteristics.
   *Remediation:* check and record `gls` convergence diagnostics
   explicitly (e.g., compare fitted log-likelihood against a
   restarted fit, or flag boundary/near-singular correlation
   estimates), and report the two failure modes (hard error vs.
   converged-but-suspect) separately.

7. **The hazard-ratio truth value (`log(0.75)`) used to compute Cox
   bias and coverage is asserted, not derived or empirically
   calibrated within the simulation.** *Location:* `report.Rmd` line
   381 ("Treatment effect: 25% slowing of progression (hazard ratio
   $\approx$ 0.75 for the survival endpoint)"); used as `true_cox =
   log(0.75)` in the `headline-13` chunk (line 539) and in
   `summarize_results()` calls throughout. *Inspected.* There is no
   analytic derivation in the manuscript connecting the linear
   slope-reduction parameter $\beta_3=-0.075$ to a hazard ratio of
   exactly 0.75 for the derived, staged, MAR-censored conversion-time
   outcome; the true hazard ratio implied by this DGM (nonlinear
   staging thresholds plus severity-dependent dropout) is not in
   general equal to a simple function of the mean-model slope ratio,
   and is not a well-defined constant under a non-proportional-hazards
   DGM unless verified. The small reported Cox bias (-0.002) is
   consistent with 0.75 being approximately correct for this
   parameterization, but "approximately correct because bias came out
   small" is circular as a derivation. *Why it matters:* Bias and
   coverage for the Cox estimand are only interpretable relative to a
   correctly specified true value; an underived nominal target
   understates or overstates apparent bias depending on how far it is
   from the DGM's actual truth. *Remediation:* either derive the
   implied marginal hazard ratio analytically (e.g., via a large,
   separate high-replication run with no MAR dropout and continuous
   CDR resolution, or via numerical integration of the DGM), or
   explicitly relabel the reported "bias" as "difference from a
   nominal design target" rather than bias against ground truth, and
   discuss the approximation's adequacy.

8. **Empty/boilerplate data documentation and empty `figures/` and
   `tables/` directories inconsistent with the zzcollab scaffold.**
   *Location:* `analysis/data/README.md` (lines 1-142) is the
   unmodified zzcollab template describing the Palmer Penguins
   dataset — species, bill length, flipper length — with placeholder
   fields such as `[Analyst name, email]` and `[YYYY-MM-DD]`, and
   `analysis/figures/`, `analysis/tables/` are both empty (all output
   is embedded via knitr caching under `analysis/report/cache/` and
   `report_files/` instead). *Inspected.* *Why it matters:* this is a
   data/code availability and provenance gap required for submission:
   nothing in the repository documents that all data are simulated
   (no real patient data), the DGM parameter provenance ("calibrated
   to observed CDR distributions in ADNI," `report.Rmd` line 333, is
   asserted with no citation, table, or code showing the calibration),
   or what `sim_13.rds` (148 KB, in `analysis/data/`, never referenced
   by any script found) actually contains. *Remediation:* replace the
   data README with an accurate description of the simulated data
   generation pipeline and a data/code availability statement for the
   manuscript itself (required by most statistics journals); either
   document or remove `sim_13.rds`; either populate or remove the
   unused `figures/`/`tables/` directories.

## 3. Minor issues

1. **No test coverage of the simulation machinery.**
   `inst/tinytest/test_basic.R` contains only `expect_true(TRUE)`
   (*inspected*). None of `sim_trial()`, `fit_mmrm()`, `fit_cox()`,
   or `summarize_results()` in `sim_study.R` — the code that produces
   every number in the manuscript — has a unit test (e.g., verifying
   the L-contrast picks the correct coefficient name, verifying
   `summarize_results()` MCSE formulas against a hand-computed toy
   case, verifying `fit_mmrm` degrades gracefully with all-dropout
   data). Given that the manuscript's own audit history documents a
   real, previously shipped estimand-labeling bug (Major issue 2 in
   the April audit, now fixed in code), the absence of a regression
   test for exactly that bug is a missed opportunity and a
   reproducibility risk for future edits.

2. **Density plot overlays estimators on different scales without a
   dual axis or standardization.** `report.Rmd` lines 623-647: MMRM
   estimates (CDR-SB change units) and Cox log-hazard-ratio estimates
   are plotted on one shared x-axis, with the caption acknowledging
   "the comparison is qualitative." A referee is likely to ask why
   this figure is included at all if it cannot support a quantitative
   claim; consider standardizing both to a common scale (e.g.,
   estimate divided by its own SE) or dropping the figure in favor of
   a table of standardized effect sizes.

3. **No sensitivity analysis on the CDR staging thresholds or
   dropout-severity coupling coefficient.** The breakpoints
   `c(-Inf, 0.25, 0.75, 1.5, 2.5, Inf)` and the dropout hazard
   multiplier `(1 + 0.5*(y_star - 0.5))` (`sim_study.R` lines 38-47)
   are stated as calibrated but presented as point values with no
   tornado/sensitivity sweep, despite the Discussion explicitly
   flagging "the relative efficiency depends on the threshold chosen"
   (line 693) as a caveat that is never investigated empirically.

4. **Seed/version/hardware disclosure incomplete.** The primary
   simulation pins `RNGkind("L'Ecuyer-CMRG")` and `set.seed(20260310)`
   (good practice, *verified*), but the manuscript nowhere states the
   R version, package versions (`nlme`, `survival`, `MASS`), or
   wall-clock/hardware used, despite `renv/` being present in the
   repository (an `renv.lock` reference or session-info appendix would
   close this gap cheaply).

5. **"ADNI" is invoked as calibration evidence with no citation.**
   `report.Rmd` line 333 ("threshold parameters calibrated to observed
   CDR distributions in ADNI") and line 247 area citing Donohue and
   Li's ADNI-based work for the *motivation*, but no ADNI data-use
   citation, acknowledgment statement, or table of the actual observed
   CDR-stage proportions used for calibration is provided.

6. **Citation-style/British-spelling artifacts in code comments carried
   into repository documentation.** `sim_study.R` line 94 uses
   "parameterisation" (British spelling) in a comment; not
   manuscript-facing but inconsistent with the US-English standard the
   author otherwise follows in the prose.

7. **Table 1 column labels are raw R variable names.** `report.Rmd`
   lines 570-584: the `kable()` call for the summary table does not
   supply `col.names`, so the rendered table (`report.tex` line 752)
   shows raw internal names (`n_valid`, `mcse_bias`, `mean_se`,
   `mcse_power`, `mcse_coverage`) rather than publication-quality
   labels. Cosmetic but relevant to "Presentation" review criteria.

## 4. What remains to be done

**(a) Required for correctness**

1. Reconcile the visit-schedule/baseline description in Methods with
   the actual `sim_trial()` code (Major issue 2), or fix the code to
   match the stated design.
2. Remove or correct the stale ADEMP self-audit embedded in the
   manuscript (Major issue 1); if retained, it must reflect the
   current, largely-remediated state of the code.
3. Rewrite the Discussion to reference the actual computed results in
   present/past tense rather than "we expect" hedging (Major issue 3).
4. Derive or explicitly caveat the `log(0.75)` Cox truth value used
   for bias/coverage (Major issue 7).
5. Verify and report genuine optimizer convergence diagnostics for the
   MMRM fits, not merely absence of a thrown error (Major issue 6).

**(b) Required for acceptance**

6. Add an abstract (Major issue 4).
7. Describe the power-vs-sample-size sweep in Methods (design, $n$
   grid, replications, seed, resulting MCSE) and make its
   reproducibility non-optional in the report build (Major issue 5).
8. Add a data/code availability statement and correct the
   `analysis/data/README.md` (Major issue 8).
9. Add minimal unit tests for `sim_trial`, `fit_mmrm`, `fit_cox`, and
   `summarize_results`, at least covering the historical estimand-label
   bug as a regression test (Minor issue 1).
10. Add session/version/seed disclosure appendix (Minor issue 4).
11. Survey and add missing literature the surveyed material implies but
    does not cite: restricted mean survival time / RMST as a competing
    continuous-scale survival summary (Royston and Parmar), a joint
    longitudinal-survival modeling reference to support the "Future
    research" joint-modeling item (e.g., Rizopoulos), and an
    interval-censored survival methods reference to support the
    "alternative survival formulations" item (e.g., Sun, or the
    `icenReg` methodology paper). These are directly relevant to claims
    made in Future Research but currently uncited.

**(c) Desirable polish**

12. Standardize or drop the dual-scale density plot (Minor issue 2).
13. Add sensitivity analysis on CDR staging thresholds and the
    dropout-severity coupling parameter (Minor issue 3).
14. Add an ADNI calibration citation/table (Minor issue 5).
15. Clean up table column labels via `col.names` (Minor issue 7).
16. Fix stray British spelling in code comments (Minor issue 6).
17. Populate or remove the empty `analysis/figures/` and
    `analysis/tables/` directories; document or remove
    `analysis/data/sim_13.rds`.

## 5. Recommended framing

**Two components, one paper, unclear identity.** The manuscript
currently reads as two loosely joined pieces: (i) a historical/
regulatory synthesis of why AD trials use MMRM and oncology trials
use survival analysis, and (ii) a small simulation study quantifying
the power cost of dichotomizing a continuous AD endpoint into a
survival endpoint. Neither piece, taken alone, would clear a
statistics journal bar as currently written; the framing decision
determines which piece is the paper and which becomes background.

**Plausible framings for what this material could become:**

- **(A) Historical/regulatory review paper** (à la a *Statistics in
  Medicine* or *Clinical Trials* perspective piece). This framing
  would foreground the Introduction's synthesis and treat the
  simulation as one illustrative example among several. Problem: the
  literature already covers this ground thoroughly and is cited
  extensively by the manuscript itself — Donohue et al. (2012),
  Mallinckrodt et al. (2001, 2003), the NRC (2010) panel report, and
  ICH E9(R1) (2019) already establish the MMRM-vs-LOCF history and
  the change-score-vs-event-time efficiency argument in the cited
  Donohue (2011) and Li (2019) papers. A pure historical review adds
  synthesis value but not a new empirical result, and the "what
  practitioners currently use" question (MMRM/PMRM as primary,
  survival essentially unused in AD) is already settled and
  uncontested in the field — there is no live debate this review
  would resolve.

- **(B) Simulation study establishing/refining the relative-efficiency
  result under a specific, realistic AD DGM** (à en a *Statistics in
  Medicine*, *Pharmaceutical Statistics*, or *JBS* methods paper).
  This is the framing the ADEMP Methods section, the Results table,
  and the two figures are already built for. The genuine gap here is
  narrower than "MMRM beats Cox" (already shown analytically and
  empirically by Donohue and Sabbagh 2011 and Li et al. 2019 in real
  ADNI data) — it is instead the specific question of efficiency loss
  when a *staged categorical instrument* (CDR global stage, not a
  continuous surrogate) is discretized into a *conversion-time*
  endpoint under *outcome-dependent MAR dropout*, which the cited
  prior work does not simulate under this precise combination. This
  is a real, defensible, narrower contribution if the DGM
  documentation and estimand issues (Major issues 2, 6, 7) are fixed.

- **(C) Software/tools paper** around a reusable simulation and power
  utility for AD trial design (foreshadowed by Future Research item
  5). Not viable now: there is no `R/` package code (the `R/`
  directory is empty; *verified* via directory listing), no exported
  functions, no documentation, and no vignette — only analysis scripts.
  This framing would require substantial new development and is better
  treated as a distinct, future deliverable rather than a redirection
  of the current manuscript.

**Recommendation: framing (B), narrowed and re-titled.** The
literature survey the authors have already assembled shows that (i)
the historical "why MMRM in AD, why survival in oncology" narrative
is settled and well-covered, and (ii) the general claim that
continuous repeated-measures analysis is more efficient than
time-to-event analysis of the same disease process is already
established by Donohue and Sabbagh (2011) and Li et al. (2019) using
real ADNI data. The paper's actual value-add, once the correctness
issues above are fixed, is the specific simulation evidence — under a
transparent, fully reported ADEMP design — of how large that
efficiency gap is for the CDR global-stage MCI-to-dementia conversion
endpoint specifically, incorporating realistic severity-dependent MAR
dropout and a categorical (not continuous surrogate) staging
instrument, which is a materially different DGM from the ADNI
observational analyses in the cited prior work.

**Implications of framing (B):**

- *Title:* de-emphasize "Historical Review" (which oversells the
  Introduction's contribution) in favor of something like "Efficiency
  Loss from Dichotomizing a Longitudinal AD Severity Measure into a
  Time-to-Conversion Endpoint: A Simulation Study," with the
  historical material repositioned as motivating background.
- *Abstract/Introduction:* state directly, with citation, that the
  MMRM-over-LOCF and continuous-over-event-time efficiency questions
  are each already resolved in the literature, and that this paper's
  contribution is the ADEMP-compliant quantification for the specific
  CDR-staging DGM with dropout — not a rediscovery of either point.
  Trim the "fall of LOCF" and "regulatory origins" subsections (lines
  95-208) substantially; they are accurate but are restating settled,
  well-cited history rather than motivating a new result, and at
  current length crowd out space that should go to DGM justification
  and sensitivity analysis.
- *Comparators:* the current MMRM-vs-Cox comparison is fair but thin
  for a methods paper. Add at least one alternative continuous-outcome
  competitor already flagged as relevant by the Future Research
  section — e.g., the PMRM of Jönsson et al. (2024), already cited —
  or a restricted-mean-survival-time alternative on the survival side,
  so the comparison is not a two-method horse race but a positioned
  comparison against the current state of the art in each camp.
- *Target journal:* given the narrowed methodological framing,
  *Statistics in Medicine*, *Pharmaceutical Statistics*, or
  *Statistical Methods in Medical Research* are better fits than JASA
  or Biometrics — the contribution is an applied simulation
  demonstration relevant to trial design practice, not a new
  general statistical method or asymptotic theory.
- *Material to move to supplement:* the historical/regulatory
  narrative (Introduction sections "Regulatory origins," "The fall of
  LOCF," "MMRM as the primary analysis in AD," "Why survival analysis
  persists in oncology") should be compressed to a half-page motivating
  paragraph in the main text with the full narrative moved to
  supplementary material; the ADEMP self-audit (once corrected) belongs
  in supplementary material or a reproducibility appendix, never inline
  under "References"; the power-vs-$n$ sweep, once properly documented
  in Methods, can stay in the main text since it is the paper's most
  practically useful deliverable (a sample-size planning figure).

## 6. Assessment

**Major revision**, bordering on reject-and-resubmit if the framing
change in Section 5 is adopted. Justification: the simulation
machinery is fundamentally sound and largely well-engineered — proper
ADEMP scaffolding, Monte Carlo SEs, a pinned RNG stream, a documented
history of at least one caught and fixed estimand-labeling bug, and a
clean render with real (not placeholder) numbers flowing into the
manuscript (*verified* via `report.tex`). That is real, creditable
work. But the manuscript as it stands cannot go to a referee in its
present state: it contains an internally contradictory, stale
self-audit that tells the reader not to trust the very results the
paper reports; a Methods description of the trial visit schedule that
does not match the code that generated the data; a Discussion section
apparently un-updated since before the simulation was run; no
abstract; and a headline figure whose generating simulation is
undocumented and not guaranteed to reproduce from a clean build. None
of these are fatal to the underlying science, but each is the kind of
finding that stops a referee from being able to certify correctness,
and collectively they indicate the manuscript was assembled across
multiple disjoint editing passes without a final consistency pass.
Once the correctness items in Section 4(a) are fixed and the framing
in Section 5 is adopted to sharpen the contribution against the
already-strong cited literature, this is a publishable, useful applied
simulation paper for a pharmaceutical/medical statistics journal.

## 7. Revision history

- 2026-08-16: Initial review. No prior `pub_review_whitepaper_*.md`
  file existed in `docs/`; this is the first referee-style pass on
  `analysis/report/report.Rmd`. Identified 8 major issues (stale
  self-contradictory ADEMP audit embedded in the manuscript; visit
  schedule/baseline description mismatched with code; Discussion
  written in unrevised future tense despite computed results being
  present; missing abstract; undocumented and non-reproducing
  power-vs-$n$ sweep figure; convergence metric that does not verify
  actual optimizer convergence; undeived Cox truth value; boilerplate/
  irrelevant data documentation) and 7 minor issues (no unit tests for
  simulation code, dual-scale density plot, no threshold sensitivity
  analysis, incomplete version/seed disclosure, uncited ADNI
  calibration claim, stray British spelling, unlabeled table columns).
  Recommended narrowing the paper to a simulation-methods framing
  (Section 5, option B) rather than the current historical-review-plus-
  simulation hybrid, on the grounds that the historical narrative and
  the general MMRM-over-survival efficiency claim are both already
  well established in the cited literature (Donohue and Sabbagh 2011;
  Li et al. 2019; Mallinckrodt et al. 2001, 2003; NRC 2010), while the
  specific CDR-staging-with-dropout DGM quantification is the paper's
  genuine, defensible contribution.
