# Morris et al. (2019) ADEMP Audit: 13-mmrm-vs-survival-ad
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/scripts/sim_study.R`
- `analysis/report/report.Rmd`

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | comparison described narratively; no ADEMP header |
| DGMs documented | Met | longitudinal + survival DGM in `sim_study.R` |
| Factors varied factorially | Partial | scenario grid not labelled |
| Estimand defined with true value | Partial | treatment effect parameterised but extractor may mislabel |
| Methods justified | Met | MMRM vs survival AD compared |
| Performance measures justified | Partial | power and bias planned but not rendered |
| n_sim stated | Partial | script sets a default; rendered report doesn't run sim |
| n_sim justified via MCSE | Not met | no derivation |
| MCSE reported per metric | Not met | not computed |
| Seed set once | Partial | seed set in script but chunks `eval=FALSE` so unverified in report |
| RNG states stored | Not met | not stored |
| Paired comparisons | Met | methods applied to same dataset per rep (when sim runs) |
| Reproducibility | **Not met** | most report chunks `eval=FALSE`; placeholder text at report.Rmd:493 |

## Overall verdict

**Not compliant.**

## Gaps

- Report chunks `sim-function`, `analysis-functions`, `run-simulation`,
  `summary-table`, `power-comparison`, `density-plot` are all
  `eval=FALSE` (`report.Rmd` around line 493).
- **Estimand extractor bug risk** at `sim_study.R:92`:
  `grep("^trt$", rownames(fe))` with a fallback
  `grep("trt", rownames(fe))[1]` at line 94 will silently pick
  `trt:visit_f2` (a specific interaction term) when the model specifies
  `trt * visit_f` and no main-effect row is named exactly `trt`. This
  can mislabel the reported treatment effect.
- No Monte Carlo SE on any metric.
- Placeholder text at `report.Rmd:493` indicates rendered report
  describes a simulation never run.
- `RNGkind()` not pinned.

## Remediation plan

1. Fix the estimand extractor in `sim_study.R:92-94`: replace the
   `grep("trt")[1]` fallback with an explicit contrast definition.
   Recommended pattern: construct the treatment-effect estimate via
   `emmeans::emmeans(fit, ~ trt | visit_f)` at the primary endpoint
   visit, or pull the coefficient by exact name without fallback.
2. Activate the six `eval=FALSE` chunks in `report.Rmd`. Cache outputs
   to `analysis/data/derived_data/`.
3. Remove the placeholder narrative at `report.Rmd:493`.
4. Add `mcse_*` columns (bias, empirical SE, coverage, power) per
   Morris Table 6.
5. Add an n_sim justification block deriving the chosen n_rep from a
   target MCSE.
6. Pin `RNGkind("L'Ecuyer-CMRG")` and store `.Random.seed` per rep.
7. Add ADEMP Methods section to `report.Rmd`.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/13-mmrm-vs-survival-ad/mmrmsurvival/docs/morris-audit-2026-04-17.md*
