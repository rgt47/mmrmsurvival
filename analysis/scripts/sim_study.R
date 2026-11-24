library(MASS)
library(nlme)
library(survival)
library(tidyverse)

# Sections delimited by `## ---- label ----` are read into the
# manuscript with knitr::read_chunk(), so the code shown in the paper
# is the code that runs. Sourcing this file only defines functions.

## ---- sim-function ----
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
                      dropout_rate = 0.15,
                      dropout = c("latent", "mcar"),
                      dropout_coupling = 0.5,
                      breaks = c(-Inf, 0.25, 0.75, 1.5, 2.5, Inf),
                      screen = TRUE) {
  dropout <- match.arg(dropout)
  n <- 2 * n_per_arm
  stages <- c(0, 0.5, 1, 2, 3)
  inner <- breaks[-c(1, length(breaks))]
  to_stage <- function(y) {
    stages[findInterval(y, inner, left.open = TRUE) + 1]
  }

  # Month 0 is the screening/randomization visit. With screen = TRUE,
  # candidates are drawn until 2 * n_per_arm are at global CDR 0.5,
  # which is the entry criterion of an MCI-due-to-AD trial. Treatment
  # is assigned after screening, so month 0 carries no treatment term.
  re <- matrix(numeric(0), ncol = 2)
  y0 <- numeric(0)
  while (length(y0) < n) {
    m <- 2 * (n - length(y0)) + 10
    re_c <- matrix(mvrnorm(m, mu = c(0, 0), Sigma = sigma_b), ncol = 2)
    y0_c <- beta[1] + re_c[, 1] + rnorm(m, 0, sigma_e)
    keep <- if (screen) to_stage(y0_c) == 0.5 else rep(TRUE, m)
    re <- rbind(re, re_c[keep, , drop = FALSE])
    y0 <- c(y0, y0_c[keep])
  }
  re <- re[seq_len(n), , drop = FALSE]
  y0 <- y0[seq_len(n)]

  trt <- rep(0:1, each = n_per_arm)
  times <- c(0, seq_len(n_visits) * 0.5)
  y <- matrix(NA_real_, n, n_visits + 1)
  y[, 1] <- y0
  for (j in seq_len(n_visits)) {
    t_j <- times[j + 1]
    y[, j + 1] <- beta[1] + beta[2] * t_j +
      beta[3] * trt + beta[4] * (t_j * trt) +
      re[, 1] + re[, 2] * t_j + rnorm(n, 0, sigma_e)
  }

  # Dropout is monotone. Months 0 and 6 are always observed; after
  # each post-baseline visit a participant leaves with probability
  # p_drop over the next 6 months. Under 'latent' the hazard rises
  # with the latent severity at that visit (MAR given y_star); under
  # 'mcar' it is constant.
  observed <- matrix(TRUE, n, n_visits + 1)
  active <- rep(TRUE, n)
  for (j in seq_len(n_visits)) {
    observed[, j + 1] <- active
    mult <- if (dropout == "latent") {
      pmax(0, 1 + dropout_coupling * (y[, j + 1] - 0.5))
    } else {
      1
    }
    p_drop <- 1 - exp(-dropout_rate * 0.5 * mult)
    active <- active & (runif(n) >= p_drop)
  }

  out <- data.frame(
    id = rep(seq_len(n), times = n_visits + 1),
    trt = rep(trt, times = n_visits + 1),
    visit = rep(0:n_visits, each = n),
    time = rep(times, each = n),
    y_star = as.vector(y),
    cdr = to_stage(as.vector(y))
  )
  out <- out[as.vector(observed), ]
  out <- out[order(out$id, out$visit), ]
  rownames(out) <- NULL
  out
}

## ---- analysis-functions ----
# Estimand: difference between arms in mean change from month 0 to
# the last visit, beta_trt + beta_trt:visit_f_last. The month-0
# value is pre-randomization, so adjusting for it does not change
# the target. Its coefficient varies by visit (baseline * visit_f)
# because the regression of change on baseline differs by visit
# under random slopes; a single pooled slope would misspecify the
# mean, which matters once dropout depends on the outcome. Under
# screening, the observed stage at month 0 is 0.5 for everyone, so
# the baseline terms are dropped for outcome = 'cdr'.
fit_mmrm <- function(dat, outcome = "y_star", max_iter = 200) {
  na_out <- c(est = NA, se = NA, p = NA)
  dat <- dat |>
    mutate(value = .data[[outcome]]) |>
    group_by(id) |>
    mutate(
      baseline = value[match(0, visit)],
      change = value - baseline
    ) |>
    ungroup() |>
    filter(visit > 0, !is.na(baseline)) |>
    mutate(visit_f = factor(visit))
  if (nrow(dat) == 0 || nlevels(dat$visit_f) < 2) return(na_out)

  use_baseline <- isTRUE(stats::var(dat$baseline) > 0)
  form <- if (use_baseline) {
    change ~ trt * visit_f + baseline * visit_f
  } else {
    change ~ trt * visit_f
  }

  mod <- tryCatch(
    gls(
      form,
      data = dat,
      correlation = corSymm(form = ~ 1 | id),
      weights = varIdent(form = ~ 1 | visit_f),
      na.action = na.omit,
      control = glsControl(maxIter = max_iter, opt = "optim")
    ),
    error = function(e) NULL
  )
  if (is.null(mod)) return(na_out)

  # gls() can finish without error at a non-optimum. A character
  # apVar signals a non-positive-definite Hessian.
  apvar_ok <- is.matrix(mod$apVar)
  iter_ok <- is.null(mod$numIter) || mod$numIter < max_iter
  if (!apvar_ok || !iter_ok) return(na_out)

  beta_hat <- coef(mod)
  last_level <- levels(dat$visit_f)[nlevels(dat$visit_f)]
  trt_int <- paste0("trt:visit_f", last_level)
  if (!all(c("trt", trt_int) %in% names(beta_hat))) return(na_out)
  L <- setNames(rep(0, length(beta_hat)), names(beta_hat))
  L[c("trt", trt_int)] <- 1

  est <- sum(L * beta_hat)
  se <- sqrt(as.numeric(t(L) %*% vcov(mod) %*% L))
  # Normal (Wald z) reference distribution; see Methods.
  p <- 2 * stats::pnorm(-abs(est / se))
  c(est = est, se = se, p = p)
}

# Time to the first post-baseline visit with global CDR >= threshold.
# Conversion is known only to lie within a 6-month interval; it is
# recorded at the visit where it is first seen.
make_surv_data <- function(dat, threshold = 1) {
  dat |>
    filter(visit > 0) |>
    group_by(id) |>
    summarise(
      trt = first(trt),
      event_time = {
        conv <- which(cdr >= threshold)
        if (length(conv) > 0) time[conv[1]] else max(time)
      },
      event = as.numeric(any(cdr >= threshold)),
      .groups = "drop"
    )
}

fit_cox <- function(dat, threshold = 1) {
  surv_dat <- make_surv_data(dat, threshold)
  mod <- tryCatch(
    coxph(Surv(event_time, event) ~ trt, data = surv_dat),
    error = function(e) NULL
  )
  if (is.null(mod)) return(c(est = NA, se = NA, p = NA))
  s <- summary(mod)$coefficients
  c(
    est = unname(coef(mod)),
    se = unname(s[, "se(coef)"]),
    p = unname(s[, "Pr(>|z|)"])
  )
}

## ---- run-functions ----
# Per-trial descriptors of the simulated population, by arm (0, 1):
# event fraction, fraction converting at the first post-baseline
# visit, and fraction completing the last visit.
trial_descriptors <- function(dat, threshold = 1) {
  sd_ <- make_surv_data(dat, threshold)
  last <- max(dat$visit)
  completers <- unique(dat$id[dat$visit == last])
  comp <- tapply(sd_$id %in% completers, sd_$trt, mean)
  ev <- tapply(sd_$event, sd_$trt, mean)
  ev1 <- tapply(sd_$event == 1 & sd_$event_time == 0.5, sd_$trt, mean)
  c(
    ev_pbo = ev[["0"]], ev_act = ev[["1"]],
    ev1_pbo = ev1[["0"]], ev1_act = ev1[["1"]],
    comp_pbo = comp[["0"]], comp_act = comp[["1"]]
  )
}

# Each replication gets its own L'Ecuyer-CMRG stream, derived in
# sequence from the caller's seed. Results therefore do not depend on
# n_cores, and the caller's stream is advanced past all of them, so
# consecutive calls do not reuse random numbers.
rep_streams <- function(n_reps) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") {
    stop("run_simulation() requires RNGkind(\"L'Ecuyer-CMRG\")")
  }
  if (!exists(".Random.seed", envir = globalenv())) stats::runif(1)
  s <- get(".Random.seed", envir = globalenv())
  streams <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    s <- parallel::nextRNGStream(s)
    streams[[r]] <- s
  }
  assign(".Random.seed", parallel::nextRNGStream(s), envir = globalenv())
  streams
}

run_simulation <- function(n_reps = 2000,
                           n_per_arm = 200,
                           n_cores = 1L,
                           ...) {
  streams <- rep_streams(n_reps)
  parent_seed <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", parent_seed, envir = globalenv()))
  one_rep <- function(r) {
    assign(".Random.seed", streams[[r]], envir = globalenv())
    dat <- sim_trial(n_per_arm = n_per_arm, ...)
    m_lat <- fit_mmrm(dat, "y_star")
    m_cdr <- fit_mmrm(dat, "cdr")
    cx <- fit_cox(dat)
    data.frame(
      rep = r,
      mmrm_est = m_lat[["est"]], mmrm_se = m_lat[["se"]],
      mmrm_p = m_lat[["p"]],
      mmrm_cdr_est = m_cdr[["est"]], mmrm_cdr_se = m_cdr[["se"]],
      mmrm_cdr_p = m_cdr[["p"]],
      cox_est = cx[["est"]], cox_se = cx[["se"]],
      cox_p = cx[["p"]],
      as.list(trial_descriptors(dat))
    )
  }
  results <- if (n_cores > 1L) {
    parallel::mclapply(seq_len(n_reps), one_rep, mc.cores = n_cores)
  } else {
    lapply(seq_len(n_reps), one_rep)
  }
  # A killed mclapply child returns NULL, which rbind would drop.
  failed <- vapply(results, function(x) {
    is.null(x) || inherits(x, "try-error")
  }, logical(1))
  if (any(failed)) stop(sum(failed), " replications failed")
  do.call(rbind, results)
}

# Power at each per-arm sample size; bias and coverage are not
# needed here, so no truth values are passed.
run_power_sweep <- function(sample_sizes = c(100, 150, 200, 300, 400),
                            n_reps = 200,
                            n_cores = 1L,
                            ...) {
  purrr::map_dfr(sample_sizes, function(n) {
    smy <- summarize_results(
      run_simulation(n_reps, n, n_cores, ...),
      truth = c()
    )
    tibble::tibble(
      n = n, method = smy$method,
      power = smy$power, n_valid = smy$n_valid
    )
  })
}

# Estimand values under the data-generating mechanism, from one very
# large trial with no dropout: the latent MMRM contrast, the mean
# difference in observed global CDR change at the last visit, and the
# large-sample limit of the Cox log hazard ratio.
large_sample_targets <- function(n_per_arm = 2e5, ...) {
  dat <- sim_trial(n_per_arm = n_per_arm, dropout_rate = 0, ...)
  last <- max(dat$visit)
  chg <- dat |>
    group_by(id) |>
    summarise(
      trt = first(trt),
      d_lat = y_star[visit == last] - y_star[visit == 0],
      d_cdr = cdr[visit == last] - cdr[visit == 0],
      .groups = "drop"
    )
  cox_mod <- coxph(
    Surv(event_time, event) ~ trt,
    data = make_surv_data(dat)
  )
  c(
    mmrm = mean(chg$d_lat[chg$trt == 1]) - mean(chg$d_lat[chg$trt == 0]),
    mmrm_cdr = mean(chg$d_cdr[chg$trt == 1]) -
      mean(chg$d_cdr[chg$trt == 0]),
    cox = unname(coef(cox_mod)),
    cox_se = unname(sqrt(diag(vcov(cox_mod)))),
    trial_descriptors(dat)
  )
}

method_labels <- c(
  mmrm = "MMRM, latent score",
  mmrm_cdr = "MMRM, global CDR stage",
  cox = "Cox, time to CDR >= 1"
)

# Monte Carlo SEs follow Morris, White and Crowther (2019), Table 6.
# `truth` is a named vector keyed by method prefix; a method with no
# entry (or NA) gets NA bias and coverage.
summarize_results <- function(sim_results,
                              truth = c(mmrm = -0.075 * 2.5),
                              alpha = 0.05) {
  methods <- intersect(names(method_labels),
                       sub("_est$", "", grep("_est$",
                         names(sim_results), value = TRUE)))
  rows <- lapply(methods, function(m) {
    est_all <- sim_results[[paste0(m, "_est")]]
    se_all <- sim_results[[paste0(m, "_se")]]
    p_all <- sim_results[[paste0(m, "_p")]]
    ok <- !is.na(est_all)
    n_ok <- sum(ok)
    est <- est_all[ok]
    tv <- if (m %in% names(truth)) unname(truth[[m]]) else NA_real_
    sd_est <- if (n_ok >= 2) stats::sd(est) else NA_real_
    pw <- if (n_ok >= 1) mean(p_all[ok] < alpha) else NA_real_
    cov_i <- (est_all - 1.96 * se_all <= tv) &
      (est_all + 1.96 * se_all >= tv)
    cv <- if (n_ok >= 1 && !is.na(tv)) mean(cov_i[ok]) else NA_real_
    data.frame(
      method = unname(method_labels[m]),
      n_valid = n_ok,
      bias = if (n_ok >= 1) mean(est) - tv else NA_real_,
      mcse_bias = if (n_ok >= 2) sd_est / sqrt(n_ok) else NA_real_,
      emp_se = sd_est,
      mcse_emp_se = if (n_ok >= 2) sd_est / sqrt(2 * (n_ok - 1))
                    else NA_real_,
      mean_se = if (n_ok >= 1) mean(se_all[ok]) else NA_real_,
      power = pw,
      mcse_power = if (n_ok >= 1) sqrt(pw * (1 - pw) / n_ok)
                   else NA_real_,
      coverage = cv,
      mcse_coverage = if (!is.na(cv)) sqrt(cv * (1 - cv) / n_ok)
                      else NA_real_,
      convergence = n_ok / nrow(sim_results)
    )
  })
  do.call(rbind, rows)
}

# Wilson score interval for a binomial proportion; unlike the Wald
# interval it does not collapse to zero width at 0 or 1.
wilson_ci <- function(p, n, z = 1.96) {
  mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  cbind(lower = mid - half, upper = mid + half)
}
