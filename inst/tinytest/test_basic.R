library(tinytest)

# Tests for the simulation code that produces the manuscript's numbers
# (analysis/scripts/sim_study.R). That file lives under
# analysis/scripts/, not R/, so it is sourced directly.

find_sim_study <- function() {
  d <- getwd()
  for (i in 1:6) {
    candidate <- file.path(d, "analysis", "scripts", "sim_study.R")
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NULL
}

sim_study_path <- find_sim_study()

if (is.null(sim_study_path)) {
  exit_file(
    "analysis/scripts/sim_study.R not found; skipping simulation tests"
  )
}

suppressPackageStartupMessages({
  have_deps <- requireNamespace("MASS", quietly = TRUE) &&
    requireNamespace("nlme", quietly = TRUE) &&
    requireNamespace("survival", quietly = TRUE) &&
    requireNamespace("tidyverse", quietly = TRUE)
})

if (!have_deps) {
  exit_file("simulation dependencies not available; skipping")
}

suppressPackageStartupMessages(source(sim_study_path))

## ---- sim_trial(): structure, screening, dropout ------------------------

set.seed(1)
dat <- sim_trial(n_per_arm = 40)

expect_true(
  all(c("id", "trt", "visit", "time", "y_star", "cdr") %in% names(dat)),
  info = "sim_trial() output has the columns downstream code expects"
)
expect_equal(
  as.vector(table(dat$trt[dat$visit == 0])), c(40L, 40L),
  info = "40 participants per arm are enrolled"
)
expect_true(
  all(dat$cdr %in% c(0, 0.5, 1, 2, 3)),
  info = "observed CDR is confined to the staged category set"
)
expect_true(
  all(dat$cdr[dat$visit == 0] == 0.5),
  info = "screening enrolls only participants at CDR 0.5 at month 0"
)
expect_equal(
  sort(unique(dat$time)), c(0, 0.5, 1, 1.5, 2, 2.5),
  info = "visit grid is months 0, 6, 12, 18, 24, 30"
)
expect_true(
  all(table(dat$id[dat$visit %in% 0:1]) == 2),
  info = "months 0 and 6 are always observed"
)
visits_by_id <- split(dat$visit, dat$id)
expect_true(
  all(vapply(visits_by_id, function(v) {
    identical(sort(v), seq.int(0L, length(v) - 1L))
  }, logical(1))),
  info = "dropout is absorbing: observed visits are 0..k for each id"
)

set.seed(2)
dat_unscreened <- sim_trial(n_per_arm = 500, screen = FALSE)
expect_true(
  any(dat_unscreened$cdr[dat_unscreened$visit == 0] != 0.5),
  info = "screen = FALSE does not condition on the month-0 stage"
)

set.seed(3)
dat_complete <- sim_trial(n_per_arm = 30, dropout_rate = 0)
expect_equal(
  nrow(dat_complete), 60L * 6L,
  info = "dropout_rate = 0 gives complete data"
)

set.seed(4)
dat_mcar <- sim_trial(n_per_arm = 30, dropout = "mcar")
expect_true(
  nrow(dat_mcar) < 60L * 6L && nrow(dat_mcar) >= 60L * 2L,
  info = "MCAR dropout removes some, but not the first two, visits"
)

## ---- fit_mmrm(): estimand is the last-visit contrast -------------------
## An earlier implementation returned the `trt` main effect alone, which
## under trt * visit_f is the effect at the reference visit. This test
## refits the model by hand and checks the returned value is the
## last-visit linear combination and differs from the trt-only value.

set.seed(42)
dat_est <- sim_trial(n_per_arm = 60)
res_est <- fit_mmrm(dat_est)

dat_manual <- dat_est |>
  group_by(id) |>
  mutate(
    baseline = y_star[visit == 0],
    change = y_star - baseline
  ) |>
  ungroup() |>
  filter(visit > 0) |>
  mutate(visit_f = factor(visit))
mod_manual <- gls(
  change ~ trt * visit_f + baseline * visit_f,
  data = dat_manual,
  correlation = corSymm(form = ~1 | id),
  weights = varIdent(form = ~1 | visit_f),
  control = glsControl(maxIter = 200, opt = "optim")
)
b <- coef(mod_manual)
expected_est <- unname(b["trt"] + b["trt:visit_f5"])

expect_equal(
  unname(res_est["est"]), expected_est,
  tolerance = 1e-6,
  info = "fit_mmrm() returns beta_trt + beta_trt:visit_f5"
)
expect_false(
  isTRUE(all.equal(unname(res_est["est"]), unname(b["trt"]),
                   tolerance = 0.01)),
  info = "fit_mmrm() estimate differs from the trt-only coefficient"
)

res_cdr <- fit_mmrm(dat_est, "cdr")
expect_true(
  all(is.finite(res_cdr)),
  info = "MMRM on the observed stage fits with a constant baseline"
)

res_degraded <- fit_mmrm(dat_est[dat_est$visit == 0, ])
expect_true(
  all(is.na(res_degraded)),
  info = "fit_mmrm() returns NA when no post-baseline data exist"
)

## ---- make_surv_data() and trial_descriptors() on hand-built data ------

toy_dat <- data.frame(
  id = rep(1:3, times = c(4, 3, 6)),
  trt = rep(c(0, 1, 1), times = c(4, 3, 6)),
  visit = c(0:3, 0:2, 0:5),
  time = c(0:3, 0:2, 0:5) * 0.5,
  y_star = 0,
  cdr = c(0.5, 0.5, 1, 1,  0.5, 0.5, 0.5,  0.5, 1, 1, 2, 2, 2)
)
sd_toy <- make_surv_data(toy_dat)
expect_equal(
  sd_toy$event_time, c(1, 1, 0.5),
  info = "event time is the first visit at CDR >= 1, else last visit"
)
expect_equal(
  sd_toy$event, c(1, 0, 1),
  info = "event indicator marks conversion to CDR >= 1"
)
desc_toy <- trial_descriptors(toy_dat)
expect_equal(
  unname(desc_toy[c("ev_pbo", "ev_act", "ev1_act", "comp_act")]),
  c(1, 0.5, 0.5, 0.5),
  info = "descriptors: event, month-6 event, and completion by arm"
)

## ---- summarize_results(): Morris Table 6 formulas ----------------------

toy <- data.frame(
  rep = 1:5,
  mmrm_est = c(-0.10, -0.14, -0.16, -0.20, -0.12),
  mmrm_se  = rep(0.05, 5),
  mmrm_p   = c(0.04, 0.01, 0.005, 0.001, 0.06),
  cox_est  = c(-0.20, -0.30, -0.25, -0.28, -0.35),
  cox_se   = rep(0.15, 5),
  cox_p    = c(0.20, 0.03, 0.09, 0.04, 0.01)
)
smy <- summarize_results(toy, truth = c(mmrm = -0.15, cox = log(0.75)))
mmrm_row <- smy[smy$method == method_labels[["mmrm"]], ]
cox_row <- smy[smy$method == method_labels[["cox"]], ]

expect_equal(
  mmrm_row$bias, mean(toy$mmrm_est) - (-0.15),
  tolerance = 1e-8, info = "bias = mean(estimate) - true value"
)
expect_equal(
  mmrm_row$emp_se, sd(toy$mmrm_est),
  tolerance = 1e-8, info = "empirical SE = sd of estimates"
)
expect_equal(
  mmrm_row$mcse_bias, sd(toy$mmrm_est) / sqrt(5),
  tolerance = 1e-8, info = "MCSE(bias) = empirical SE / sqrt(n_sim)"
)
expect_equal(
  mmrm_row$mcse_emp_se, sd(toy$mmrm_est) / sqrt(2 * 4),
  tolerance = 1e-8, info = "MCSE(emp SE) = emp SE / sqrt(2(n_sim - 1))"
)
expect_equal(
  mmrm_row$power, mean(toy$mmrm_p < 0.05),
  tolerance = 1e-8, info = "power = proportion of p < alpha"
)
expect_equal(
  cox_row$coverage, mean(abs(toy$cox_est - log(0.75)) <= 1.96 * 0.15),
  tolerance = 1e-8, info = "coverage of the stated truth"
)

smy_no_truth <- summarize_results(toy, truth = c())
expect_true(
  all(is.na(smy_no_truth$bias)) && all(is.na(smy_no_truth$coverage)),
  info = "without a truth value, bias and coverage are NA"
)

all_na <- data.frame(
  rep = 1:3,
  mmrm_est = NA_real_, mmrm_se = NA_real_, mmrm_p = NA_real_,
  cox_est = c(-0.2, -0.3, -0.1), cox_se = rep(0.1, 3),
  cox_p = c(0.01, 0.02, 0.5)
)
smy_na <- summarize_results(all_na, truth = c(mmrm = -0.15, cox = 0))
expect_true(
  is.na(smy_na$bias[smy_na$method == method_labels[["mmrm"]]]),
  info = "bias is NA (not an error) when every MMRM fit failed"
)
expect_equal(
  smy_na$convergence[smy_na$method == method_labels[["mmrm"]]], 0,
  info = "convergence rate is 0 when every MMRM fit failed"
)

## ---- wilson_ci() ------------------------------------------------------

ci <- wilson_ci(1, 200)
expect_true(
  ci[, "upper"] > 0.9999 && ci[, "lower"] < 0.99,
  info = "Wilson interval at p = 1 has nonzero width"
)

## ---- large_sample_targets(): latent target matches beta_3 * 2.5 -------

set.seed(5)
tg <- large_sample_targets(n_per_arm = 2e4)
expect_equal(
  unname(tg[["mmrm"]]), -0.075 * 2.5, tolerance = 0.02, scale = 1,
  info = "latent change difference at month 30 is beta_3 * 2.5"
)
expect_true(
  tg[["cox"]] < 0 && tg[["mmrm_cdr"]] < 0,
  info = "stage and Cox targets favor the active arm"
)

## ---- run_simulation(): per-replication RNG streams --------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(6)
r1 <- run_simulation(n_reps = 2, n_per_arm = 30, n_cores = 2L)
r3 <- run_simulation(n_reps = 2, n_per_arm = 30, n_cores = 1L)
set.seed(6)
r2 <- run_simulation(n_reps = 2, n_per_arm = 30, n_cores = 1L)
RNGkind("default")
expect_identical(
  r1, r2,
  info = "same seed gives identical replications for any n_cores"
)
expect_false(
  isTRUE(all.equal(r1$mmrm_est, r3$mmrm_est)),
  info = "consecutive calls do not reuse random numbers"
)
expect_true(
  all(c("mmrm_cdr_est", "cox_est", "ev_pbo", "comp_act") %in% names(r1)),
  info = "run_simulation() records all three methods and descriptors"
)
