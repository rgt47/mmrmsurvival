library(MASS)
library(nlme)
library(survival)
library(tidyverse)

sim_trial <- function(n_per_arm = 200,
                      n_visits = 5,
                      beta = c(0.5, 0.30, 0, -0.075),
                      sigma_b = matrix(
                        c(
                          0.04,
                          0.3 * sqrt(0.04 * 0.02),
                          0.3 * sqrt(0.04 * 0.02),
                          0.02
                        ),
                        nrow = 2
                      ),
                      sigma_e = 0.1,
                      dropout_rate = 0.15) {
  n <- 2 * n_per_arm
  trt <- rep(0:1, each = n_per_arm)
  times <- seq(0.5, 2.5, by = 0.5)

  re <- mvrnorm(n, mu = c(0, 0), Sigma = sigma_b)

  records <- vector("list", n * n_visits)
  idx <- 0L
  for (i in seq_len(n)) {
    dropped <- FALSE
    for (j in seq_len(n_visits)) {
      if (dropped) next
      t_j <- times[j]
      mu_ij <- beta[1] + beta[2] * t_j +
        beta[3] * trt[i] + beta[4] * (t_j * trt[i]) +
        re[i, 1] + re[i, 2] * t_j
      y_star <- mu_ij + rnorm(1, 0, sigma_e)

      cdr_obs <- cut(
        y_star,
        breaks = c(-Inf, 0.25, 0.75, 1.5, 2.5, Inf),
        labels = c(0, 0.5, 1, 2, 3)
      )
      cdr_obs <- as.numeric(as.character(cdr_obs))

      p_drop <- 1 - exp(
        -dropout_rate * 0.5 * (1 + 0.5 * (y_star - 0.5))
      )
      if (runif(1) < p_drop) dropped <- TRUE

      idx <- idx + 1L
      records[[idx]] <- data.frame(
        id = i,
        trt = trt[i],
        visit = j,
        time = t_j,
        y_star = y_star,
        cdr = cdr_obs
      )
    }
  }
  do.call(rbind, records[seq_len(idx)])
}

fit_mmrm <- function(dat) {
  dat <- dat |>
    group_by(id) |>
    mutate(
      baseline = y_star[visit == min(visit)],
      change = y_star - baseline
    ) |>
    ungroup() |>
    filter(visit > 1) |>
    mutate(visit_f = factor(visit))

  mod <- tryCatch(
    gls(
      change ~ trt * visit_f + baseline,
      data = dat,
      correlation = corSymm(form = ~ 1 | id),
      weights = varIdent(form = ~ 1 | visit_f),
      na.action = na.omit,
      control = glsControl(
        maxIter = 200, opt = "optim"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(mod)) return(c(est = NA, se = NA, p = NA))

  # Estimand: treatment effect at the LAST visit, expressed as the
  # linear combination beta_trt + beta_{trt:visit_f_last}. The previous
  # implementation picked `trt` alone, which under the trt * visit_f
  # parameterisation corresponds to the treatment effect at the
  # REFERENCE visit, not the last visit. That mismatched the
  # reported true value (-0.15, i.e. beta[4] * (2.5 - 0.5) = -0.15).
  tt <- summary(mod)$tTable
  coef_names <- rownames(tt)
  last_visit_level <- max(levels(dat$visit_f))
  trt_main <- "trt"
  trt_int <- paste0("trt:visit_f", last_visit_level)

  if (!(trt_main %in% coef_names) ||
      !(trt_int %in% coef_names)) {
    return(c(est = NA, se = NA, p = NA))
  }

  L <- rep(0, length(coef_names))
  names(L) <- coef_names
  L[trt_main] <- 1
  L[trt_int] <- 1

  beta_hat <- tt[, "Value"]
  vcov_mat <- vcov(mod)
  est <- sum(L * beta_hat)
  se <- sqrt(as.numeric(t(L) %*% vcov_mat %*% L))
  z <- est / se
  p <- 2 * stats::pnorm(-abs(z))
  c(est = est, se = se, p = p)
}

fit_cox <- function(dat) {
  surv_dat <- dat |>
    group_by(id) |>
    summarise(
      trt = first(trt),
      event_time = {
        conv <- which(cdr >= 1)
        if (length(conv) > 0) time[conv[1]] else max(time)
      },
      event = as.numeric(any(cdr >= 1)),
      .groups = "drop"
    )

  mod <- tryCatch(
    coxph(Surv(event_time, event) ~ trt, data = surv_dat),
    error = function(e) NULL
  )

  if (is.null(mod)) return(c(est = NA, se = NA, p = NA))

  s <- summary(mod)
  # coef(mod) is a named numeric (name = "trt"), so bare c() would
  # produce `est.trt`, `se.trt`, `p.trt`. Strip names before c() so
  # downstream code can index by the simple names est/se/p.
  c(
    est = unname(coef(mod)),
    se = unname(s$coefficients[, "se(coef)"]),
    p = unname(s$coefficients[, "Pr(>|z|)"])
  )
}

run_simulation <- function(n_reps = 2000,
                           n_per_arm = 200) {
  # Morris et al. (2019) §4.1: the RNG seed is set ONCE by the caller;
  # this function does NOT call `set.seed()`. Per-replicate RNG states
  # are captured and attached as an attribute of the return value for
  # diagnostic reproducibility of any failing rep.
  results <- vector("list", n_reps)
  rng_states <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
    rng_states[[r]] <- .Random.seed
    if (r %% 100 == 0) {
      message(sprintf("Replication %d / %d", r, n_reps))
    }
    dat <- sim_trial(n_per_arm = n_per_arm)
    mmrm_res <- fit_mmrm(dat)
    cox_res <- fit_cox(dat)
    results[[r]] <- data.frame(
      rep = r,
      mmrm_est = mmrm_res["est"],
      mmrm_se = mmrm_res["se"],
      mmrm_p = mmrm_res["p"],
      cox_est = cox_res["est"],
      cox_se = cox_res["se"],
      cox_p = cox_res["p"],
      row.names = NULL
    )
  }
  out <- do.call(rbind, results)
  attr(out, "rng_states") <- rng_states
  out
}

summarize_results <- function(sim_results,
                              true_mmrm = -0.075 * 2,
                              true_cox = log(0.75),
                              alpha = 0.05) {
  # Monte Carlo SE formulas per Morris, White & Crowther (2019)
  # Table 6. Degenerate cases (n < 1 or < 2) yield NA_real_.
  safe_mean <- function(x) if (length(x) >= 1) mean(x) else NA_real_
  safe_sd <- function(x) if (length(x) >= 2) stats::sd(x) else NA_real_
  safe_div <- function(num, den) {
    if (isTRUE(is.finite(num) && is.finite(den) && den > 0)) {
      num / den
    } else NA_real_
  }

  mmrm_valid <- !is.na(sim_results$mmrm_est)
  cox_valid <- !is.na(sim_results$cox_est)
  n_mmrm <- sum(mmrm_valid)
  n_cox <- sum(cox_valid)

  mmrm_est <- sim_results$mmrm_est[mmrm_valid]
  cox_est <- sim_results$cox_est[cox_valid]

  mmrm_cov <- (
    sim_results$mmrm_est - 1.96 * sim_results$mmrm_se <= true_mmrm &
      sim_results$mmrm_est + 1.96 * sim_results$mmrm_se >= true_mmrm
  )
  cox_cov <- (
    sim_results$cox_est - 1.96 * sim_results$cox_se <= true_cox &
      sim_results$cox_est + 1.96 * sim_results$cox_se >= true_cox
  )

  data.frame(
    method = c("MMRM", "Cox"),
    n_valid = c(n_mmrm, n_cox),
    bias = c(
      if (n_mmrm >= 1) safe_mean(mmrm_est) - true_mmrm else NA_real_,
      if (n_cox >= 1) safe_mean(cox_est) - true_cox else NA_real_
    ),
    mcse_bias = c(
      safe_div(safe_sd(mmrm_est), sqrt(n_mmrm)),
      safe_div(safe_sd(cox_est), sqrt(n_cox))
    ),
    emp_se = c(safe_sd(mmrm_est), safe_sd(cox_est)),
    mcse_emp_se = c(
      safe_div(safe_sd(mmrm_est), sqrt(2 * (n_mmrm - 1))),
      safe_div(safe_sd(cox_est), sqrt(2 * (n_cox - 1)))
    ),
    mean_se = c(
      mean(sim_results$mmrm_se[mmrm_valid]),
      mean(sim_results$cox_se[cox_valid])
    ),
    power = c(
      mean(sim_results$mmrm_p[mmrm_valid] < alpha),
      mean(sim_results$cox_p[cox_valid] < alpha)
    ),
    mcse_power = c(
      sqrt(
        mean(sim_results$mmrm_p[mmrm_valid] < alpha) *
          (1 - mean(sim_results$mmrm_p[mmrm_valid] < alpha)) /
          n_mmrm
      ),
      sqrt(
        mean(sim_results$cox_p[cox_valid] < alpha) *
          (1 - mean(sim_results$cox_p[cox_valid] < alpha)) /
          n_cox
      )
    ),
    coverage = c(
      mean(mmrm_cov, na.rm = TRUE),
      mean(cox_cov, na.rm = TRUE)
    ),
    mcse_coverage = c(
      sqrt(
        mean(mmrm_cov, na.rm = TRUE) *
          (1 - mean(mmrm_cov, na.rm = TRUE)) / n_mmrm
      ),
      sqrt(
        mean(cox_cov, na.rm = TRUE) *
          (1 - mean(cox_cov, na.rm = TRUE)) / n_cox
      )
    ),
    convergence = c(n_mmrm / nrow(sim_results),
                    n_cox / nrow(sim_results))
  )
}

# The driver block that previously auto-ran 2000 replications when
# sourced non-interactively was moved to `analysis/scripts/01_run_simulation.R`.
# Sourcing this file now only defines the functions; it does not run
# the simulation. This prevents an accidental 2000-rep run when the
# report sources `sim_study.R` from its setup chunk.
