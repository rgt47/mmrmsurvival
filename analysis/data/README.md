# Data Directory

This directory holds simulated data products for the
MMRM-vs-Cox AD trial simulation study
(`analysis/report/report.Rmd`). **No real patient data are
used anywhere in this project.**

`raw_data/` and `derived_data/` are unused, empty
placeholders inherited from the zzcollab project template.

## The report does not read these files

Since the 2026-11-19 revision, `report.Rmd` computes every
reported number, table, and figure during the render: the
estimand values (`targets` chunk), the primary simulation
(`run-simulation`), the secondary scenarios
(`run-scenarios`), and the sample-size sweep (`run-sweep`).
None of the `.rds` files below is read by the report.

That revision also changed the data-generating mechanism
(a month-0 screening visit with a CDR 0.5 entry criterion,
change from month 0, a baseline-by-visit term in the MMRM,
an MMRM on the observed global CDR stage, and per-replication
RNG streams). All three files below predate it and do not
describe the current design.

| File | Producing script | Status |
|------|------------------|--------|
| `sim_results.rds` | `analysis/scripts/01_run_simulation.R` | Stale (April 2026 code state; 0 valid Cox fits). Regenerate with `Rscript analysis/scripts/01_run_simulation.R`. |
| `sim_power_vs_n.rds` | `analysis/scripts/02_power_vs_n.R` | Stale (May 2026; two methods only). Regenerate with `Rscript analysis/scripts/02_power_vs_n.R`. |
| `sim_13.rds` | no current script | Archival snapshot of the pre-revision run (MMRM power 0.994, Cox power 0.716 under the old design). Not reproducible from current code. |

Both scripts must be run from the package root. They use
the same seeds as the report and so reproduce its primary
results and sweep.

## Data/code availability statement

All simulation parameters, data-generating code, and
analysis code for this manuscript are contained in this
`mmrmsurvival` research compendium
(`analysis/scripts/sim_study.R`,
`analysis/scripts/01_run_simulation.R`,
`analysis/scripts/02_power_vs_n.R`,
`analysis/report/report.Rmd`), with package versions pinned
via `renv.lock`. No human-subjects data were used. The
results are reproduced by rendering the report with
`bash tools/render.sh analysis/report/report.Rmd`.
