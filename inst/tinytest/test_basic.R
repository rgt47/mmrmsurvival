library(tinytest)

# Regression tests for the simulation machinery that produces the
# manuscript's reported numbers (analysis/scripts/sim_study.R). This
# file replaces the placeholder `expect_true(TRUE)` flagged by
# pub_review_whitepaper_2026-08-16.md (Minor issue 1, checklist item
# 9). `sim_study.R` lives under `analysis/scripts/`, not `R/`, so it
# is `source()`d directly rather than tested via the installed
# package namespace.

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
    requireNamespace("dplyr", quietly = TRUE)
})

if (!have_deps) {
  exit_file(
    "MASS/nlme/survival/dplyr not all available; skipping simulation tests"
  )
}

suppressPackageStartupMessages({
  library(MASS)
  library(nlme)
  library(survival)
  library(dplyr)
})
source(sim_study_path)

## ---- sim_trial(): structural checks ---------------------------------

set.seed(1)
dat <- sim_trial(n_per_arm = 20)

expect_true(
  all(c("id", "trt", "visit", "time", "y_star", "cdr") %in% names(dat)),
  info = "sim_trial() output has the columns downstream code expects"
)
expect_true(
  all(dat$trt %in% c(0, 1)),
  info = "treatment indicator is coded 0/1"
)
expect_true(
  all(dat$cdr %in% c(0, 0.5, 1, 2, 3)),
  info = "observed CDR is confined to the staged category set"
)
expect_true(
  nrow(dat) <= 40 * 5,
  info = "row count cannot exceed n_subjects * n_visits (dropout removes rows, never adds)"
)
## explicit monotone-dropout check: for each id, observed visits must
## be a contiguous run starting at visit 1 (dropout is absorbing).
visits_by_id <- split(dat$visit, dat$id)
monotone_ok <- vapply(visits_by_id, function(v) {
  identical(sort(v), seq_len(length(v)))
}, logical(1))
expect_true(
  all(monotone_ok),
  info = "dropout is absorbing: each id's observed visits are 1..k for some k"
)

## ---- fit_mmrm(): estimand-labeling regression test -------------------
## Historical bug (docs/morris-audit-2026-04-17.md): an earlier
## implementation picked the `trt` main-effect coefficient alone,
## which under `trt * visit_f` is the effect at the REFERENCE visit,
## not the last visit. This test fits the model both ways on the same
## data and confirms fit_mmrm() returns the LAST-visit contrast
## (beta_trt + beta_trt:visit_f_last), not the reference-visit
## coefficient alone, and that the two differ materially (so the test
## would have caught the historical bug).

set.seed(42)
dat_est <- sim_trial(n_per_arm = 60)
res_est <- fit_mmrm(dat_est)

dat_manual <- dat_est |>
  group_by(id) |>
  mutate(
    baseline = y_star[visit == min(visit)],
    change = y_star - baseline
  ) |>
  ungroup() |>
  filter(visit > 1) |>
  mutate(visit_f = factor(visit))
mod_manual <- gls(
  change ~ trt * visit_f + baseline,
  data = dat_manual,
  correlation = corSymm(form = ~1 | id),
  weights = varIdent(form = ~1 | visit_f),
  na.action = na.omit,
  control = glsControl(maxIter = 200, opt = "optim")
)
tt_manual <- summary(mod_manual)$tTable
trt_only_est <- unname(tt_manual["trt", "Value"])
last_visit_level <- max(levels(dat_manual$visit_f))
L <- rep(0, nrow(tt_manual))
names(L) <- rownames(tt_manual)
L["trt"] <- 1
L[paste0("trt:visit_f", last_visit_level)] <- 1
expected_est <- sum(L * tt_manual[, "Value"])

expect_equal(
  unname(res_est["est"]), expected_est,
  tolerance = 1e-6,
  info = "fit_mmrm() reproduces the explicit last-visit linear combination"
)
expect_false(
  isTRUE(all.equal(unname(res_est["est"]), trt_only_est, tolerance = 0.01)),
  info = paste(
    "fit_mmrm() estimate must differ materially from the historical",
    "trt-only (reference-visit) bug's estimate on this dataset"
  )
)

## ---- fit_mmrm(): graceful degradation on unfittable data --------------

dat_no_postbaseline <- dat_est[dat_est$visit == 1, ]
res_degraded <- fit_mmrm(dat_no_postbaseline)
expect_true(
  all(is.na(res_degraded)),
  info = "fit_mmrm() returns NA (not an error) when no post-baseline data survive filtering"
)

## ---- summarize_results(): Monte Carlo SE formulas vs. hand calc -------
## Morris, White & Crowther (2019) Table 6 formulas, verified against
## a small hand-computed toy case.

toy <- data.frame(
  rep = 1:5,
  mmrm_est = c(-0.10, -0.14, -0.16, -0.20, -0.12),
  mmrm_se  = rep(0.05, 5),
  mmrm_p   = c(0.04, 0.01, 0.005, 0.001, 0.06),
  cox_est  = c(-0.20, -0.30, -0.25, -0.28, -0.35),
  cox_se   = rep(0.15, 5),
  cox_p    = c(0.20, 0.03, 0.09, 0.04, 0.01)
)
smy <- summarize_results(toy, true_mmrm = -0.15, true_cox = log(0.75))
mmrm_row <- smy[smy$method == "MMRM", ]
cox_row <- smy[smy$method == "Cox", ]

expect_equal(
  mmrm_row$bias, mean(toy$mmrm_est) - (-0.15),
  tolerance = 1e-8, info = "MMRM bias = mean(estimate) - true value"
)
expect_equal(
  mmrm_row$emp_se, sd(toy$mmrm_est),
  tolerance = 1e-8, info = "empirical SE = sd of estimates across reps"
)
expect_equal(
  mmrm_row$mcse_bias, sd(toy$mmrm_est) / sqrt(5),
  tolerance = 1e-8, info = "MCSE(bias) = empirical SE / sqrt(n_sim)"
)
expect_equal(
  mmrm_row$power, mean(toy$mmrm_p < 0.05),
  tolerance = 1e-8, info = "power = proportion of p < alpha"
)
expect_equal(
  cox_row$power, mean(toy$cox_p < 0.05),
  tolerance = 1e-8, info = "Cox power computed with the same formula"
)
expect_equal(
  smy$n_valid, c(5, 5),
  info = "n_valid counts non-NA estimates per method"
)

## ---- summarize_results(): degenerate n handling ------------------------

all_na <- data.frame(
  rep = 1:3,
  mmrm_est = NA_real_, mmrm_se = NA_real_, mmrm_p = NA_real_,
  cox_est = c(-0.2, -0.3, -0.1), cox_se = rep(0.1, 3),
  cox_p = c(0.01, 0.02, 0.5)
)
smy_na <- summarize_results(all_na, true_mmrm = -0.15, true_cox = log(0.75))
expect_true(
  is.na(smy_na$bias[smy_na$method == "MMRM"]),
  info = "bias is NA (not an error) when every MMRM fit failed"
)
expect_equal(
  smy_na$convergence[smy_na$method == "MMRM"], 0,
  info = "convergence rate is 0 when every MMRM fit failed"
)
