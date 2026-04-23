#!/usr/bin/env Rscript

# Driver script to run the MMRM-vs-Cox simulation. Produces
# `analysis/data/sim_results.rds` and prints Morris Table 6 summary.

suppressPackageStartupMessages({
  library(MASS)
  library(nlme)
  library(survival)
  library(tidyverse)
})

source(file.path(
  here::here(), "analysis", "scripts", "sim_study.R"
))

# Morris et al. (2019) §4.1: RNGkind pinned and seed set once.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260310)

message("Running simulation with 2000 replications...")
sim_results <- run_simulation(n_reps = 2000)

out_path <- file.path(
  here::here(), "analysis", "data", "sim_results.rds"
)
dir.create(
  dirname(out_path), showWarnings = FALSE, recursive = TRUE
)
saveRDS(sim_results, out_path)
message(sprintf("Results saved to %s", out_path))

metrics <- summarize_results(sim_results)
print(metrics)
