## Power-vs-N sweep for project 13.
## Run once; caches result to analysis/data/sim_power_vs_n.rds.
## Modest scope: n_reps_sweep = 200 replicates per sample size.

suppressPackageStartupMessages({
  library(MASS)
  library(nlme)
  library(survival)
  library(tidyverse)
})

source('analysis/scripts/sim_study.R')

set.seed(20260511)

sample_sizes <- c(100, 150, 200, 300, 400)
n_reps_sweep <- 200

power_rows <- purrr::map_dfr(sample_sizes, function(n) {
  message(sprintf('--- n_per_arm = %d ---', n))
  sim_n <- run_simulation(
    n_reps = n_reps_sweep,
    n_per_arm = n
  )
  smy <- summarize_results(
    sim_n,
    true_mmrm = -0.075 * 2,
    true_cox  = log(0.75)
  )
  tibble::tibble(
    n      = n,
    method = smy$method,
    power  = smy$power,
    mcse_power = smy$mcse_power
  )
})

out_path <- 'analysis/data/sim_power_vs_n.rds'
saveRDS(power_rows, out_path)
message(sprintf('wrote %s', out_path))
print(power_rows)
