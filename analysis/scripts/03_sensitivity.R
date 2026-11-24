#!/usr/bin/env Rscript

# One-at-a-time sensitivity analysis for project 13 (whitepaper
# 2026-08-16, Minor issue 3 / checklist 13): varies the dropout-
# severity coupling and the latent cut point between global CDR 0.5
# and CDR 1.0, holding everything else at the primary scenario. Not
# yet read by report.Rmd. Run from the package root:
#
#   Rscript analysis/scripts/03_sensitivity.R [n_reps] [n_target] \
#     [n_cores] [out_path]
#
# Defaults: n_reps = 500, n_target = 2e5, n_cores = 1, and
# out_path = analysis/data/sim_sensitivity.rds. The grid values below
# are illustrative and should be confirmed by the author.

suppressPackageStartupMessages({
  library(MASS)
  library(nlme)
  library(survival)
  library(tidyverse)
})

source(file.path(here::here(), 'analysis', 'scripts', 'sim_study.R'))

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 500L
n_target <- if (length(args) >= 2) as.numeric(args[2]) else 2e5
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 1L
out_path <- if (length(args) >= 4) {
  args[4]
} else {
  file.path(here::here(), 'analysis', 'data', 'sim_sensitivity.rds')
}

base_breaks <- c(-Inf, 0.25, 0.75, 1.5, 2.5, Inf)
with_cut <- function(cut) replace(base_breaks, 3, cut)

# The primary scenario (coupling 0.5, cut 0.75) is included so the
# grid has an internal reference computed with the same R.
grid <- tibble::tribble(
  ~factor,    ~level, ~dropout_coupling, ~breaks,
  'primary',  'base', 0.5,               list(base_breaks),
  'coupling', '0',    0,                 list(base_breaks),
  'coupling', '1',    1,                 list(base_breaks),
  'cut',      '0.65', 0.5,               list(with_cut(0.65)),
  'cut',      '0.85', 0.5,               list(with_cut(0.85))
)

RNGkind("L'Ecuyer-CMRG")
set.seed(20261123)

t0 <- Sys.time()
rows <- purrr::pmap(grid, function(factor, level, dropout_coupling,
                                   breaks) {
  brk <- unlist(breaks)
  stopifnot(length(brk) == length(base_breaks))
  targets <- large_sample_targets(
    n_per_arm = n_target, breaks = brk
  )
  res <- run_simulation(
    n_reps = n_reps, n_cores = n_cores,
    dropout_coupling = dropout_coupling, breaks = brk
  )
  smy <- summarize_results(res, truth = targets)
  message(sprintf('%s = %s done (%.1f min elapsed)', factor, level,
                  as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
  list(
    summary = tibble::tibble(factor = factor, level = level, smy),
    targets = targets,
    results = res
  )
})

out <- list(
  grid = grid,
  summary = purrr::map_dfr(rows, 'summary'),
  targets = purrr::map(rows, 'targets'),
  results = purrr::map(rows, 'results'),
  n_reps = n_reps,
  n_target = n_target,
  seed = 20261123
)
stopifnot(all(out$summary$n_valid > 0))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, out_path)
message(sprintf('wrote %s', out_path))
print(out$summary[, c('factor', 'level', 'method', 'bias',
                      'power', 'coverage', 'n_valid')])
