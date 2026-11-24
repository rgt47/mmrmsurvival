# Referee report: 'MMRM vs. Survival Analysis for AD Clinical Trials: A Historical Review and Simulation Study'
*2026-11-18*

Manuscript: `analysis/report/report.Rmd` (with the sourced
`analysis/scripts/sim_study.R`, `02_power_vs_n.R`, and
`references.bib`). Line numbers refer to `report.Rmd` at `HEAD`.
Evidence labels: **verified** (re-derived or executed by me),
**inspected** (read and judged), **suspected** (plausible, not
confirmed).

## 1. Summary of the manuscript

The manuscript has two parts. The first is a narrative history
of why Alzheimer's disease (AD) trials settled on MMRM analyses
of change from baseline while oncology trials settled on
time-to-event analyses. It covers the 1990 Leber guidelines,
the move away from LOCF, the 2010 NRC report, ICH E9(R1), and
recent amyloid antibody trials. The second part is a Monte
Carlo study, reported under the ADEMP framework of Morris,
White and Crowther (2019). A latent continuous severity
process with random intercepts and slopes is discretized into
global CDR stages, and severity-dependent dropout is applied.
The study then compares an MMRM on change in the latent score
with a Cox model for time to first CDR >= 1. In one scenario
(200 per arm, 25% slowing, R = 2000) the MMRM has power 0.994
and the Cox model 0.716. A secondary sweep at R = 200 per cell
shows the same ordering for 100 to 400 per arm. The abstract
concedes that the general efficiency result is known (Donohue
et al. 2011; Li et al. 2019). It claims as the contribution a
transparent quantification for the specific combination of a
staged CDR instrument, a conversion endpoint, and
outcome-dependent dropout.

## 2. Overall assessment and recommendation

**Recommendation: major revision** for a specialty venue
(e.g. *Contemporary Clinical Trials*, *Alzheimer's & Dementia:
TRCI*, *Pharmaceutical Statistics*); **reject** in its current
form for *Statistics in Medicine* or a higher-tier journal.

The paper is clearly organized and the ADEMP reporting is more
careful than most applied simulation papers. However, the
central comparison is not the one the abstract describes. The
MMRM is given the latent, nearly noise-free continuous
variable, which no trial ever observes, while the Cox model is
given the coarsened stage. The power gap therefore mixes two
different things: continuous-versus-event analysis, and
latent-versus-observed measurement. Several statements about
the simulated population (all participants enter at CDR 0.5;
conversion is modest; censoring is heavy) are contradicted by
the data-generating mechanism as coded. There is a single
scenario, no null scenario, and the MMRM arm sits at a power
ceiling, so the reported 'power ratio' carries little
information. The historical review contains at least two
citation errors that a referee in the AD field will catch at
once (lecanemab Study 201; the origin of PMRM), and it omits
the AD trials that did use time-to-event primary endpoints.
None of these problems is fatal. Each is fixable with an
observed-scale MMRM arm, a modest scenario grid, a null
scenario, and corrections to the review.

## 3. Significance and novelty

The general result is well established. Donohue et al. (2011)
derived an inflation factor for time-to-threshold versus rate
of change and checked it by simulation from ADNI-based models.
Li et al. (2019) simulated continuous outcomes from a Bayesian
joint mixed model fit to ADNI, derived progression events from
those trajectories, and found power roughly doubled for the
continuous analysis, including under informative missingness.
Note that Li et al. is a simulation study, not an analysis of
'real, observational ADNI data' as stated at lines 914-916 and
in the abstract. That makes it a closer competitor than the
manuscript allows. In addition, a 2023 *Alzheimer's &
Dementia* conference abstract from the Raket/Jönsson group
simulated amyloid-positive MCI trials and compared PMRM, MMRM,
Cox, and AFT models for time to CDR-global progression across
200 to 700 per arm and 18 to 36 months. It reported Cox and
AFT power rarely above 40% (located by search; full text not
read). This is very close to the present design and is not
cited.

The remaining novelty is therefore narrow: an ADEMP-reported
simulation with a threshold-staged CDR and outcome-dependent
dropout. Even that claim is weakened by Major comment 1,
because the continuous arm does not analyze a staged
instrument at all. In my view the paper's most useful
contribution could be the historical review, if it is made
accurate and complete, combined with a simulation that
isolates the sources of the efficiency loss (coarsening,
dichotomization at a threshold, interval censoring, informative
censoring). As it stands the simulation reads as a
confirmation of known results under a single, favorable
configuration.

Contemporariness: the manuscript does not engage with the
current alternatives that a 2026 referee would expect.
These include PMRM (Raket 2022) as a primary-analysis
competitor, time-to-event secondary endpoints in recent
pivotal trials (e.g. time to CDR-global worsening in
TRAILBLAZER-ALZ 2 and Clarity AD), estimand-based handling of
intercurrent events, and joint or multi-state models. Several
appear only in 'Future research'.

## 4. Major comments

**M1. The MMRM analyzes the latent variable, not the CDR.
(Methods, lines 344-357, 519-533, 547-556; abstract; Discussion
lines 899-909.) Inspected.** `fit_mmrm()` computes
`change = y_star - baseline`, where `y_star` is the
unobservable latent severity with residual SD 0.1. The Cox
model uses `cdr`, the five-level staging of the same latent
variable. The Methods text (line 354) says 'the primary
continuous outcome is the global CDR score observed at these
five visits', and line 302 says the global CDR is used 'as a
continuous outcome'. Both are inaccurate. The abstract's
emphasis on 'a staged categorical instrument (CDR global
stage, not a continuous surrogate)' describes the Cox arm
only. The density-plot caption (line 842) and the Discussion
(line 937) speak of CDR-SB, which is never simulated. As a
result, the power gap confounds (a) continuous versus
time-to-event analysis with (b) latent versus coarsened
measurement. The latter is an advantage no real MMRM
analysis has.
*Remedy:* add an MMRM arm on an observable outcome. Options
are the observed global CDR stage, or a simulated CDR-SB
(e.g. a finer discretization of `y_star` onto 0-18 in 0.5
steps with added rating error). Report the latent-scale MMRM
only as an upper bound. Correct the Methods, caption and
Discussion text accordingly.

**M2. The simulated population does not match its
description. (Lines 343-344, 389-391, 918-923.) Verified.**
The text says 'All participants enter with a latent global
CDR of 0.5'. In the code, 0.5 is the mean; the random
intercept has SD 0.2. Simulating the DGM with 200,000
participants and no dropout, I find that 6.6% are at CDR 0
and 31.9% are already at CDR >= 1 at the first visit (month
6). With dropout on (5,000 per arm), 35% of placebo and 29%
of active participants have their event at the first visit.
The overall event fractions are 81% (placebo) and 72%
(active). This contradicts the Discussion's statements that
'conversion rates are modest', that 'many participants ...
will not reach the CDR >= 1.0 threshold', and that 'heavy
right-censoring' limits the Cox model. In fact the Cox model
here suffers from the opposite problem: roughly a third of
the sample contributes an event at the first assessment, a
time when the treatment effect is only -0.0375 on the latent
scale. These events are driven almost entirely by baseline
heterogeneity. Note also that there is no month-0 assessment,
so CDR 0.5 at entry is never checked; a real MCI trial would
exclude such participants at screening.
*Remedy:* either impose the entry criterion (simulate a
screening visit and condition on CDR = 0.5), or describe the
population honestly. Report the event fraction and the
distribution of event times by arm. Rewrite the Discussion
explanation of the Cox power deficit to match.

**M3. One scenario, no null, and a ceiling on MMRM power.
(Results, Table 1; Discussion lines 899-913.) Inspected.**
All conclusions rest on one parameter configuration. MMRM
power is 0.994 (MCSE 0.002), so the 'power ratio' of about
1.39 is bounded by the ceiling and cannot be compared with
the 'two-fold' figure of Li et al. There is no scenario with
beta_3 = 0, so the type I error of the MMRM (with a Wald z
reference and a post-randomization baseline covariate) and of
the Cox model (with five tied event times and informative
censoring; see M5) is not shown. ADEMP asks for the
data-generating mechanisms to be varied in a factorial way.
The manuscript acknowledges this in lines 646-653 but still
draws general conclusions.
*Remedy:* add (i) a null scenario; (ii) a smaller effect
size so that MMRM power is near 0.8; (iii) at least two
values each of the residual/rating error, the CDR threshold
(e.g. 0.5 to 1 versus 0.5 to 2), and trial duration; (iv) the
relative sample size needed for 80% power, which is a more
interpretable efficiency measure than a power ratio.

**M4. The MMRM estimand and its 'true' value. (Lines
326-329, 434-435, 547-556, 736-740.) Partly verified.** The
reference 'baseline' is the month-6 value, which treatment
already affects (expected difference -0.0375). It is then
used both to form the change score and as a covariate. With
complete data and a visit-specific baseline adjustment, the
treatment contrast at month 30 is -0.135, not -0.15 (verified
by a 200,000-participant regression). In a 60-replication
run of the actual `fit_mmrm()` (seed 424242), the mean
estimate was -0.156 (MCSE 0.004), consistent with the stated
-0.15. So the reported near-zero bias holds empirically, but
for a model-specific reason: the single pooled baseline slope
in `change ~ trt * visit_f + baseline` under GLS weighting.
The stated truth therefore depends on the exact covariate
specification, which the paper does not explain. In addition,
adjusting for a post-randomization variable is not a
recommended estimand practice (ICH E9(R1)).
*Remedy:* simulate a true month-0 assessment, define change
from randomization, and state the estimand in ICH E9(R1)
terms (population, variable, intercurrent event strategy for
dropout, summary measure). The true value then becomes
beta_3 x 2.5 = -0.1875 on the latent scale, with no
dependence on the covariate structure.

**M5. Dropout is informative for the Cox model. (Lines
507-509, 919-923.) Inspected.** Dropout after visit j depends
on the latent `y_star` at visit j. For the MMRM, which
analyzes `y_star`, this is MAR and ignorable. For the Cox
model, which sees only the stage, censoring depends on
information within a stage that the model does not condition
on. The independent-censoring assumption is therefore
violated, and the direction of the resulting bias is not
obvious. The manuscript describes the dropout as MAR without
qualification, which is true only from the MMRM's point of
view. This also means the comparison under dropout is not
symmetric.
*Remedy:* state this explicitly, and add a scenario with
dropout that depends only on observed stage (or MCAR
dropout) to separate the two effects.

**M6. Cox 'truth' and the use of a nominal target.
(Lines 442-469, 745-766.) Inspected.** The authors are
candid that log(0.75) is a nominal target and that the DGM
is not proportional hazards. That candor is welcome.
However, Cox bias and coverage computed against a nominal
value still appear in Table 1 alongside the MMRM bias, and a
reader will compare them. The required computation is cheap:
simulate one very large trial (e.g. 10^6 per arm) without
dropout and fit the same Cox model, to obtain the
large-sample limit of the Cox estimator under this DGM.
*Remedy:* replace the nominal target with that limit, or
remove Cox bias and coverage from Table 1 and report only
power and empirical SE for the Cox model. Also consider an
assessment of proportional hazards (e.g. Schoenfeld
residuals averaged across replications).

**M7. Errors and omissions in the historical review.
Inspected (with web searches noted in section 8).**
(a) Lines 236-241: lecanemab Study 201 did not use 'MMRM with
change from baseline in CDR-SB as the primary endpoint'. Its
primary endpoint was change in ADCOMS at 12 months under a
Bayesian analysis; Dhadda et al. (2022) is a paper about
sensitivity analyses of that trial. The phase 3 Clarity AD
trial (van Dyck et al. 2023) is the correct example for
CDR-SB MMRM.
(b) Line 234: 'virtually every pivotal AD trial uses MMRM'.
TRAILBLAZER-ALZ 2 (Sims et al. 2023) used a natural cubic
spline model for iADRS as primary, so the claim needs
qualification.
(c) Lines 285-289: PMRM was introduced by Raket (2022,
*Statistics in Medicine* 41:5537-5557). Jönsson et al. (2024)
applied and compared it.
(d) The review omits the AD trials that used time-to-event
primary endpoints. These include the ADCS vitamin E/selegiline
trial (Sano et al. 1997, time to a composite milestone) and
the ADCS MCI trial of vitamin E and donepezil (Petersen et
al. 2005, time to AD diagnosis). These are central to the
manuscript's own question, and the author's 2000 *Statistics
in Medicine* paper on the Sano trial data is directly
relevant.
(e) Lines 177-186: the claim that placebo participants
'tend to drop out at higher rates' in AD is uncited, and in
symptomatic trials (e.g. cholinesterase inhibitors) the
active arm often had more dropout from adverse events.
(f) Line 248: 'OS is binary' should refer to death, not to
overall survival, which is a time-to-event variable.
*Remedy:* correct (a)-(c) and (f), cite (d) and discuss
why the field moved away from milestone endpoints, and cite
or qualify (e).

**M8. Reproducibility claims do not match the artifacts.
(Lines 615-656, 695-703, 725-731.) Inspected.** Line 632
states every chunk runs 'with eval=TRUE on each render, so
the reported results are generated by the rendered document
itself'. However, `cache = TRUE` is set globally, and the
cache key does not include `sim_study.R` (as the 2026-08-16
remediation log notes). The rendered `report.tex` and
`report.pdf` (dated 2026-08-15) predate the convergence-check
change to `fit_mmrm()` made on 2026-08-16. The echoed
`fit_cox()` in the Analysis plan chunk differs from
`sim_study.R` (no `tryCatch`, unnamed coefficient handling).
The sweep file `sim_power_vs_n.rds` (dated May 2026) also
predates the August code changes. `sim_results.rds` is
documented as stale. The paper also describes 'per-replication
RNG states retained', but with `L'Ecuyer-CMRG` and a serial
loop this is redundant with the master seed; it is harmless.
*Remedy:* add `cache.extra = tools::md5sum(...)` for the
sourced script, re-run the primary and sweep simulations from
the current code, and make the echoed functions identical to
the sourced ones (or echo them via `knitr::read_chunk`).
Remove changelog language from the manuscript (see m1).

## 5. Minor comments

m1. Lines 626-655, 697-700, 767-770: the text contains
revision-history language ('as of this revision', 'supersedes
the 2026-04-17 audit', 'now stops the render', 'now also
requires'). This belongs in a response letter, not in a
journal manuscript.

m2. Lines 116-119: the claim that the divergence 'is not an
accident of history' but reflects 'fundamental differences in
disease biology' is asserted without support; the review that
follows mostly describes regulatory precedent.

m3. Line 137-140: Leber's guidelines required a clinician's
global impression (CIBIC), which is not a continuous measure;
line 877-878 describes them as 'continuous cognitive and
functional co-primary endpoints'.

m4. Lines 409-410 and 507-508: the dropout is described as
'15%/year', but the hazard is scaled by
`1 + 0.5 (y* - 0.5)`, so the realized rate varies by severity
and arm. About 29-30% of participants do not complete in both
arms (verified, 5,000 per arm). Report the realized dropout.

m5. Line 333 and line 586: the MMRM p-value uses a normal
reference; with 400 participants this is adequate, but state
it, and say whether Kenward-Roger or Satterthwaite df were
considered.

m6. Line 802 (caption): 'uniformly dominates' and 'gap widest
at n = 100' match the cached sweep (gaps 0.50, 0.43, 0.31,
0.11, 0.05; verified from `sim_power_vs_n.rds`), but at
R = 200 per cell three MMRM cells sit at 0.995-1.000. The
error bars at power = 1 have zero width, which understates
uncertainty; use a Wilson interval.

m7. Line 842 (caption): 'continuous CDR-SB rate' should be
'change in latent CDR at month 30'.

m8. Future research item 1 is headed 'Additional covariance
structures' but discusses MNAR mechanisms.

m9. Lines 283 and 293 report figures from Li et al. and Ito
and Hutmacher with 'approximately'; please give the exact
quantities and the setting. The Ito and Hutmacher 'median
time to a 1-point change of about 2 years' was not verified
by me.

m10. `references.bib`: the ICH author is 'International
Council for Harmonisation' (the official spelling; keep it).
`leber1990guidelines` lacks a document identifier or URL.
The three entries added in August were flagged by the prior
remediation as unverified against a bibliographic database.
`fda2018oncology` and `fda2018earlyad` lack docket numbers or
URLs.

m11. Line 238: 'all employed MMRM ... Sensitivity analyses
from DPM and other longitudinal models were used' attributes
DPM sensitivity analyses to the aducanumab trials as well as
Study 201; please check.

m12. Line 444, 'This point warrants emphasis', and similar
openers can be removed.

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Lines 236-239 | Study 201 primary endpoint misstated | Swanson CJ et al. A randomized, double-blind, phase 2b proof-of-concept clinical trial in early Alzheimer's disease with lecanemab, an anti-Abeta protofibril antibody. *Alzheimer's Research & Therapy* 2021;13:80. van Dyck CH et al. Lecanemab in early Alzheimer's disease. *N Engl J Med* 2023;388:9-21. |
| Line 234 | 'virtually every' needs qualification | Sims JR et al. Donanemab in early symptomatic Alzheimer disease: the TRAILBLAZER-ALZ 2 randomized clinical trial. *JAMA* 2023;330:512-527. |
| Lines 285-289 | PMRM attributed to Jönsson 2024 | Raket LL. Progression models for repeated measures: estimating novel treatment effects in progressive diseases. *Statistics in Medicine* 2022;41(28):5537-5557. doi:10.1002/sim.9581 |
| Intro, history | AD trials with time-to-event primaries omitted | Sano M et al. A controlled trial of selegiline, alpha-tocopherol, or both as treatment for Alzheimer's disease. *N Engl J Med* 1997;336:1216-1222. Petersen RC et al. Vitamin E and donepezil for the treatment of mild cognitive impairment. *N Engl J Med* 2005;352:2379-2388. |
| Section 1.7 / novelty | Closest prior simulation not cited | 2023 *Alzheimer's & Dementia* abstract (Raket/Jönsson group) comparing PMRM, MMRM, Cox and AFT for time to CDR-global progression; full citation to be confirmed by the author. |
| Lines 914-916, abstract | Li 2019 described as observational ADNI data | Li et al. 2019 (already cited): describe correctly as a simulation from an ADNI-fitted model. |
| Discussion, dichotomization | General cost of dichotomizing a continuous variable | Altman DG, Royston P. The cost of dichotomising continuous variables. *BMJ* 2006;332:1080. |
| Lines 177-179 | Differential dropout by arm in AD, uncited | A cited source or a qualification is needed; no specific source suggested. |
| Discussion, informative censoring | Cox under covariate-dependent censoring | A standard reference (e.g. Kalbfleisch and Prentice, *The Statistical Analysis of Failure Time Data*, 2nd ed., Wiley 2002) |

Citations in the table other than Raket (2022) and Li et al.
(2019) are from my own knowledge and were not checked against
a bibliographic database during this review.

## 7. Suggestions for strengthening the paper

1. Decompose the efficiency loss. Run, on the same simulated
   trials: (a) MMRM on latent `y_star`; (b) MMRM on a
   simulated CDR-SB; (c) MMRM or ordinal mixed model on the
   global stage; (d) Cox on the threshold; (e) an
   interval-censored model; (f) RMST. This would show how
   much of the gap comes from measurement coarsening and how
   much from the time-to-event summary. That is a genuinely
   new and useful quantity.
2. Add PMRM (`pmrm` package) as a comparator, since it is the
   current competitor that expresses effects as time.
3. Calibrate the DGM to a published placebo-arm trajectory
   (e.g. ADNI or ADCS MCI placebo data), or report a
   sensitivity analysis over the staging thresholds.
4. Report relative sample size for 80% power across
   scenarios, rather than power at one sample size.
5. Present the sweep at R = 1000 per cell; the cost is
   modest (about 1 s per replication on my run).

## 8. Scope of this review

Verified (executed): the DGM entry-state distribution and the
first-visit and overall event fractions (scratch scripts
`check_dgm.R`, `check_cens.R`); the complete-data
baseline-adjusted contrast (-0.135); a 60-replication run of
the actual `run_simulation()` at 200 per arm (mean MMRM
estimate -0.156, Cox power 0.77, about 55 s); the cached
sweep values in `sim_power_vs_n.rds`.

Inspected: all of `report.Rmd`, `sim_study.R`,
`01_run_simulation.R`, `02_power_vs_n.R`, `references.bib`,
`analysis/data/README.md`, the 2026-08-16 remediation log,
and the Results section of `report.tex`.

Not checked: the full R = 2000 primary run was not re-run;
Table 1 values are taken from the 2026-08-15 `report.tex`. I
did not render the document. I did not read the full text of
Leber (1990), Donohue et al. (2011), Ito and Hutmacher
(2014), or the remaining cited papers; claims attributed to
them were judged only for plausibility, except as noted.

Literature searches run (web search): (1) lecanemab Study 201
primary endpoint and Dhadda 2022; (2) Raket PMRM 2022;
(3) Li, Iddi, Aisen, Thompson, Donohue 2019 findings and
design; (4) simulations comparing MMRM with time to CDR-global
progression (2020-2023). Sano 1997, Petersen 2005, Sims 2023,
van Dyck 2023 and Altman and Royston 2006 are cited from my own
knowledge and were not searched.
