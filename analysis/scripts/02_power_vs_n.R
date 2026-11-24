## Power-vs-N sweep for project 13, outside the report.
## The report computes the same sweep in its `run-sweep` chunk with
## the same seed and core count, so this script is a convenience for
## inspecting the sweep without rendering. Run from the package root.

suppressPackageStartupMessages({
  library(MASS)
  library(nlme)
  library(survival)
  library(tidyverse)
})

source(file.path(here::here(), 'analysis', 'scripts', 'sim_study.R'))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260511)

power_rows <- run_power_sweep(
  sample_sizes = c(100, 150, 200, 300, 400),
  n_reps = 200,
  n_cores = 8L
)

out_path <- file.path(
  here::here(), 'analysis', 'data', 'sim_power_vs_n.rds'
)
saveRDS(power_rows, out_path)
message(sprintf('wrote %s', out_path))
print(power_rows)
