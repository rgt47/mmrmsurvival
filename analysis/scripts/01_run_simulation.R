#!/usr/bin/env Rscript

# Primary MMRM-vs-Cox simulation, outside the report. The report's
# `targets` and `run-simulation` chunks use the same seeds and core
# count. Writes `analysis/data/sim_results.rds` and prints the Morris
# Table 6 summary. Run from the package root.

suppressPackageStartupMessages({
  library(MASS)
  library(nlme)
  library(survival)
  library(tidyverse)
})

source(file.path(
  here::here(), "analysis", "scripts", "sim_study.R"
))

RNGkind("L'Ecuyer-CMRG")
n_cores <- 8L

set.seed(20260309)
targets <- large_sample_targets(n_per_arm = 2e5)

set.seed(20260310)
message("Running simulation with 2000 replications...")
sim_results <- run_simulation(n_reps = 2000, n_cores = n_cores)

out_path <- file.path(
  here::here(), "analysis", "data", "sim_results.rds"
)
dir.create(
  dirname(out_path), showWarnings = FALSE, recursive = TRUE
)
saveRDS(list(targets = targets, sim_results = sim_results), out_path)
message(sprintf("Results saved to %s", out_path))

print(targets)
print(summarize_results(sim_results, truth = targets))
