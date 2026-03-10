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

  tt <- summary(mod)$tTable
  row_nm <- grep("^trt$", rownames(tt), value = TRUE)
  if (length(row_nm) == 0) {
    row_nm <- grep("trt", rownames(tt), value = TRUE)[1]
  }
  c(
    est = tt[row_nm, "Value"],
    se = tt[row_nm, "Std.Error"],
    p = tt[row_nm, "p-value"]
  )
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
  c(
    est = coef(mod),
    se = s$coefficients[, "se(coef)"],
    p = s$coefficients[, "Pr(>|z|)"]
  )
}

run_simulation <- function(n_reps = 2000,
                           n_per_arm = 200,
                           seed = 20260310) {
  set.seed(seed)
  results <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
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
  do.call(rbind, results)
}

if (!interactive()) {
  message("Running simulation with 2000 replications...")
  sim_results <- run_simulation(n_reps = 2000)
  out_path <- file.path(
    "analysis", "data", "sim_results.rds"
  )
  saveRDS(sim_results, out_path)
  message(sprintf("Results saved to %s", out_path))

  true_mmrm <- -0.075 * 2
  true_cox <- log(0.75)

  metrics <- sim_results |>
    summarise(
      mmrm_bias = mean(mmrm_est, na.rm = TRUE) -
        true_mmrm,
      mmrm_emp_se = sd(mmrm_est, na.rm = TRUE),
      mmrm_mean_se = mean(mmrm_se, na.rm = TRUE),
      mmrm_power = mean(mmrm_p < 0.05, na.rm = TRUE),
      mmrm_converged = mean(!is.na(mmrm_est)),
      cox_bias = mean(cox_est, na.rm = TRUE) - true_cox,
      cox_emp_se = sd(cox_est, na.rm = TRUE),
      cox_mean_se = mean(cox_se, na.rm = TRUE),
      cox_power = mean(cox_p < 0.05, na.rm = TRUE),
      cox_converged = mean(!is.na(cox_est))
    )
  print(as.data.frame(metrics))
}
