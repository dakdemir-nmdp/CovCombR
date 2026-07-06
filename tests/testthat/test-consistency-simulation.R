# Consistency Simulations: Truth-Known Recovery Across the Whole Package
#
# These tests simulate small, known ground-truth covariance structures and check
# that every major piece of CovCombR behaves the way theory predicts:
#   - Sigma_hat gets closer to the true Sigma as total degrees of freedom grow
#     (free model and k-factor model)
#   - AIC/BIC/auto factor selection recover the true number of factors as
#     sample size grows
#   - Plugin and bootstrap standard errors are calibrated: true error falls
#     within a normal-approximation interval at roughly the nominal rate
#   - print/summary/coef/fitted remain consistent with the fitted object
#
# Problem sizes are kept small (p <= 8, few dozen to few hundred df) so the
# whole file runs quickly, but large enough to show the expected asymptotics.

# -- Helpers ------------------------------------------------------------

# True covariance shared by several tests below: a single common factor plus
# idiosyncratic variance, i.e. exactly a rank-1 factor model with p = 6.
.consistency_true_sigma_k1 <- function() {
  p <- 6
  var_names <- paste0("V", 1:p)
  lambda <- c(0.9, 0.8, 0.7, -0.6, 0.5, 0.4)
  psi <- c(0.3, 0.4, 0.5, 0.4, 0.3, 0.6)
  Sigma <- outer(lambda, lambda) + diag(psi)
  dimnames(Sigma) <- list(var_names, var_names)
  list(Sigma = Sigma, lambda = lambda, psi = psi, p = p)
}

# Simulate a fixed missing-data design (3 overlapping studies covering all
# pairs at least once) at a given per-study degrees of freedom, for a given
# true Sigma. Returns S_list/nu ready for fit_covcomb().
.consistency_simulate <- function(true_Sigma, nu_per_study, seed) {
  set.seed(seed)
  var_names <- rownames(true_Sigma)
  p <- length(var_names)
  stopifnot(p == 6)

  patterns <- list(1:4, 3:6, c(1, 2, 5, 6))
  S_list <- list()
  nu_vec <- c()
  for (i in seq_along(patterns)) {
    idx <- patterns[[i]]
    W <- rWishart(1, nu_per_study, true_Sigma[idx, idx])[, , 1]
    S <- W / nu_per_study
    dimnames(S) <- list(var_names[idx], var_names[idx])
    id <- paste0("study", i)
    S_list[[id]] <- S
    nu_vec[id] <- nu_per_study
  }
  list(S_list = S_list, nu = nu_vec)
}

.max_abs_err <- function(A, B) max(abs(A - B))

# -- Free-model consistency: error shrinks as nu grows -------------------

test_that("CONSISTENCY: free model Sigma_hat error shrinks as sample size grows", {
  truth <- .consistency_true_sigma_k1()
  true_Sigma <- truth$Sigma

  nu_levels <- c(30, 120, 480)
  errs <- numeric(length(nu_levels))

  for (i in seq_along(nu_levels)) {
    dat <- .consistency_simulate(true_Sigma, nu_levels[i], seed = 3000 + i)
    fit <- suppressWarnings(suppressMessages(
      fit_covcomb(dat$S_list, dat$nu,
        scale_method = "none", se_method = "none", n_factors = NULL,
        control = list(max_iter = 300)
      )
    ))
    ids <- rownames(true_Sigma)
    errs[i] <- .max_abs_err(fit$Sigma_hat[ids, ids], true_Sigma)
  }

  cat("\nFree model: max-abs-error by nu:", errs, "\n")

  # Error should shrink (allow slight non-monotonicity from sampling noise,
  # but the largest sample size must beat the smallest by a clear margin).
  expect_lt(errs[3], errs[1])
  # Rough root-n scaling sanity check: going from nu=30 to nu=480 is a 16x
  # increase in information, so error should shrink by roughly 4x (sqrt(16));
  # allow generous slack for a single-seed draw.
  expect_lt(errs[3], errs[1] / 1.5)
  # Absolute sanity: at the largest sample size, recovery should be tight.
  expect_lt(errs[3], 0.35)
})

test_that("CONSISTENCY: free model recovers truth within plugin-SE error bars at large nu", {
  truth <- .consistency_true_sigma_k1()
  true_Sigma <- truth$Sigma
  ids <- rownames(true_Sigma)

  dat <- .consistency_simulate(true_Sigma, nu_per_study = 400, seed = 3010)
  fit <- suppressWarnings(suppressMessages(
    fit_covcomb(dat$S_list, dat$nu,
      scale_method = "none", se_method = "plugin", n_factors = NULL,
      control = list(max_iter = 300)
    )
  ))

  z <- (fit$Sigma_hat[ids, ids] - true_Sigma) / fit$Sigma_se[ids, ids]
  z_finite <- z[is.finite(z)]

  cat(sprintf(
    "\nFree model z-scores at nu=400: mean=%.3f, max|z|=%.3f, frac(|z|<3)=%.3f\n",
    mean(z_finite), max(abs(z_finite)), mean(abs(z_finite) < 3)
  ))

  # With a single realization, individual z-scores are noisy, but nearly all
  # entries should fall within a generous +/-3.5 SE band, and the average
  # magnitude should be well below the tail.
  expect_true(mean(abs(z_finite) < 3.5) > 0.9)
  expect_lt(mean(abs(z_finite)), 1.5)
})

# -- Factor-model consistency ---------------------------------------------

test_that("CONSISTENCY: k=1 factor model recovers true Sigma better as sample size grows", {
  truth <- .consistency_true_sigma_k1()
  true_Sigma <- truth$Sigma

  nu_levels <- c(30, 120, 480)
  errs <- numeric(length(nu_levels))

  for (i in seq_along(nu_levels)) {
    dat <- .consistency_simulate(true_Sigma, nu_levels[i], seed = 3100 + i)
    fit <- suppressWarnings(suppressMessages(
      fit_fa_em(dat$S_list, dat$nu, k = 1,
        scale_method = "none",
        control = list(max_iter = 300)
      )
    ))
    ids <- rownames(true_Sigma)
    errs[i] <- .max_abs_err(fit$Sigma_hat[ids, ids], true_Sigma)
  }

  cat("\nk=1 factor model: max-abs-error by nu:", errs, "\n")

  expect_lt(errs[3], errs[1])
  expect_lt(errs[3], 0.25)
  expect_true(all(is.finite(errs)))
})

test_that("CONSISTENCY: factor model reconstructs Sigma from Lambda/Psi consistently", {
  # Structural check: Sigma_hat == Lambda_hat %*% t(Lambda_hat) + diag(Psi_hat)
  # must hold exactly (this is a definitional identity, not a statistical one).
  truth <- .consistency_true_sigma_k1()
  dat <- .consistency_simulate(truth$Sigma, nu_per_study = 100, seed = 3200)

  fit <- suppressWarnings(suppressMessages(
    fit_fa_em(dat$S_list, dat$nu, k = 1, scale_method = "none")
  ))

  recon <- fit$Lambda_hat %*% t(fit$Lambda_hat) + diag(fit$Psi_hat, nrow = length(fit$Psi_hat))
  expect_equal(unname(recon), unname(fit$Sigma_hat), tolerance = 1e-8)

  # Psi (unique variances) must stay strictly positive
  expect_true(all(fit$Psi_hat > 0))
})

test_that("CONSISTENCY: k=2 factor model recovers a true 2-factor structure", {
  set.seed(3300)
  p <- 8
  var_names <- paste0("V", 1:p)
  Lambda_true <- matrix(0, p, 2)
  Lambda_true[1:4, 1] <- c(0.9, 0.8, 0.7, 0.6)
  Lambda_true[5:8, 2] <- c(0.9, -0.8, 0.7, 0.6)
  Psi_true <- rep(0.3, p)
  true_Sigma <- Lambda_true %*% t(Lambda_true) + diag(Psi_true)
  dimnames(true_Sigma) <- list(var_names, var_names)

  patterns <- list(1:6, 3:8, c(1, 2, 4, 5, 7, 8))

  .sim_k2 <- function(nu_per_study, seed) {
    set.seed(seed)
    S_list <- list()
    nu_vec <- c()
    for (i in seq_along(patterns)) {
      idx <- patterns[[i]]
      W <- rWishart(1, nu_per_study, true_Sigma[idx, idx])[, , 1]
      S <- W / nu_per_study
      dimnames(S) <- list(var_names[idx], var_names[idx])
      id <- paste0("study", i)
      S_list[[id]] <- S
      nu_vec[id] <- nu_per_study
    }
    list(S_list = S_list, nu = nu_vec)
  }

  nu_levels <- c(40, 200)
  errs <- numeric(length(nu_levels))
  for (i in seq_along(nu_levels)) {
    dat <- .sim_k2(nu_levels[i], seed = 3300 + i)
    fit <- suppressWarnings(suppressMessages(
      fit_fa_em(dat$S_list, dat$nu, k = 2, scale_method = "none",
                control = list(max_iter = 300))
    ))
    errs[i] <- .max_abs_err(fit$Sigma_hat[var_names, var_names], true_Sigma)
  }

  cat("\nk=2 factor model: max-abs-error by nu:", errs, "\n")
  expect_lt(errs[2], errs[1])
  expect_lt(errs[2], 0.4)
})

# -- Model selection: AIC/BIC/auto pick the true number of factors -------

test_that("CONSISTENCY: BIC selects the true k=1 factor model over free/other-k as sample size grows", {
  truth <- .consistency_true_sigma_k1()
  dat <- .consistency_simulate(truth$Sigma, nu_per_study = 300, seed = 3400)

  fit_auto <- suppressWarnings(suppressMessages(
    fit_covcomb(dat$S_list, dat$nu,
      scale_method = "none", se_method = "none", n_factors = "auto",
      control = list(max_iter = 300, ic = "BIC")
    )
  ))

  expect_true(!is.null(fit_auto$ic_table))
  cat("\nAuto factor-selection IC table (nu=300):\n")
  print(fit_auto$ic_table)

  # The selected model should be the k=1 factor model (true structure), not
  # an over-parameterized free-Sigma fit, once sample size is generous.
  expect_identical(fit_auto$model, "factor")
  expect_identical(fit_auto$n_factors, 1L)
  expect_true(fit_auto$ic_table$selected[fit_auto$ic_table$n_factors == 1 &
    !is.na(fit_auto$ic_table$n_factors)])
})

test_that("CONSISTENCY: compute_aic/compute_bic penalize free model more than true k=1 model", {
  truth <- .consistency_true_sigma_k1()
  dat <- .consistency_simulate(truth$Sigma, nu_per_study = 300, seed = 3410)

  fit_k1 <- suppressWarnings(suppressMessages(
    fit_fa_em(dat$S_list, dat$nu, k = 1, scale_method = "none",
              control = list(max_iter = 300))
  ))
  fit_free <- suppressWarnings(suppressMessages(
    fit_covcomb(dat$S_list, dat$nu, scale_method = "none", se_method = "none",
                n_factors = NULL, control = list(max_iter = 300))
  ))

  aic_k1 <- compute_aic(fit_k1)
  aic_free <- compute_aic(fit_free)
  bic_k1 <- compute_bic(fit_k1, dat$nu)
  bic_free <- compute_bic(fit_free, dat$nu)

  cat(sprintf(
    "\nAIC: k1=%.2f free=%.2f | BIC: k1=%.2f free=%.2f\n",
    aic_k1, aic_free, bic_k1, bic_free
  ))

  # The true, more parsimonious k=1 model should be favored by BIC (which
  # penalizes the free model's larger parameter count more heavily).
  expect_lt(bic_k1, bic_free)

  comparison <- compare_models(fit_k1, fit_free)
  expect_true(is.finite(comparison$p_value) || is.na(comparison$df))
})

# -- Standard error calibration (Monte Carlo over replicate datasets) -----

test_that("CONSISTENCY: plugin SEs are calibrated across repeated simulated datasets", {
  # Simulate many independent datasets from the same truth/design, fit with
  # plugin SEs each time, and check that the empirical coverage of a
  # nominal ~95% interval (Sigma_hat +/- 1.96*SE) is in a reasonable range.
  # This directly checks the plugin SE formula against actual sampling
  # variability, for fully-observed entries only (plugin SEs are known to be
  # anti-conservative for imputed entries, which is documented behavior, not
  # tested here).
  truth <- .consistency_true_sigma_k1()
  true_Sigma <- truth$Sigma
  ids <- rownames(true_Sigma)

  n_rep <- 60
  nu_per_study <- 120
  covered <- logical(0)

  for (r in seq_len(n_rep)) {
    dat <- .consistency_simulate(true_Sigma, nu_per_study, seed = 4000 + r)
    fit <- suppressWarnings(suppressMessages(
      fit_covcomb(dat$S_list, dat$nu,
        scale_method = "none", se_method = "plugin", n_factors = NULL,
        control = list(max_iter = 200)
      )
    ))
    # Entry (3,4) [1-indexed V3,V4] is jointly observed in study1 (1:4) and
    # study2 (3:6), so it has full coverage and a well-defined plugin SE.
    est <- fit$Sigma_hat["V3", "V4"]
    se <- fit$Sigma_se["V3", "V4"]
    truth_val <- true_Sigma["V3", "V4"]
    covered <- c(covered, abs(est - truth_val) < 1.96 * se)
  }

  coverage_rate <- mean(covered)
  cat(sprintf("\nPlugin SE empirical coverage (nominal 95%%): %.3f over %d reps\n",
              coverage_rate, n_rep))

  # Allow a wide band around 95% given only 60 replicates (Monte Carlo SE of
  # a coverage proportion at n=60 is ~2.8%), but coverage should not collapse.
  expect_gt(coverage_rate, 0.80)
})

test_that("CONSISTENCY: bootstrap SEs approximate the empirical sampling SD of Sigma_hat", {
  skip_on_cran()
  # Ground truth: generate many independent datasets, compute the *actual*
  # Monte Carlo SD of Sigma_hat across replicates, then compare to the
  # bootstrap SE reported for a single fit. They should agree in order of
  # magnitude.
  truth <- .consistency_true_sigma_k1()
  true_Sigma <- truth$Sigma
  ids <- rownames(true_Sigma)
  nu_per_study <- 150

  n_rep <- 40
  sigma_hats <- array(NA_real_, dim = c(6, 6, n_rep))
  for (r in seq_len(n_rep)) {
    dat <- .consistency_simulate(true_Sigma, nu_per_study, seed = 4100 + r)
    fit <- suppressWarnings(suppressMessages(
      fit_covcomb(dat$S_list, dat$nu,
        scale_method = "none", se_method = "none", n_factors = NULL,
        control = list(max_iter = 200)
      )
    ))
    sigma_hats[, , r] <- fit$Sigma_hat[ids, ids]
  }
  empirical_sd <- apply(sigma_hats, c(1, 2), sd)
  dimnames(empirical_sd) <- list(ids, ids)

  # Single fit with built-in bootstrap SE
  dat1 <- .consistency_simulate(true_Sigma, nu_per_study, seed = 4100 + 1)
  fit_boot <- suppressWarnings(suppressMessages(
    fit_covcomb(dat1$S_list, dat1$nu,
      scale_method = "none", se_method = "bootstrap", n_factors = NULL,
      control = list(max_iter = 200, bootstrap = list(B = 80, seed = 99))
    )
  ))
  boot_se <- fit_boot$Sigma_se[ids, ids]

  # Compare only entry V3,V4 (fully covered pair) to avoid comparing entries
  # whose empirical SD is dominated by imputation-pattern idiosyncrasies.
  ratio <- boot_se["V3", "V4"] / empirical_sd["V3", "V4"]
  cat(sprintf(
    "\nBootstrap SE vs Monte Carlo empirical SD (V3,V4): boot=%.4f empirical=%.4f ratio=%.3f\n",
    boot_se["V3", "V4"], empirical_sd["V3", "V4"], ratio
  ))

  # Bootstrap SE should be within a factor of 2 of the true Monte Carlo SD
  # (loose bound: 40 MC reps and 80 bootstrap reps both carry substantial
  # Monte Carlo noise on an SD estimate).
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)
})

# -- Log-likelihood / model comparison sanity across sample sizes --------

test_that("CONSISTENCY: normalized log-likelihood per df stabilizes as nu grows", {
  truth <- .consistency_true_sigma_k1()
  nu_levels <- c(30, 300)
  ll_per_df <- numeric(length(nu_levels))

  for (i in seq_along(nu_levels)) {
    dat <- .consistency_simulate(truth$Sigma, nu_levels[i], seed = 3500 + i)
    fit <- suppressWarnings(suppressMessages(
      fit_covcomb(dat$S_list, dat$nu, scale_method = "none", se_method = "none",
                  n_factors = NULL, control = list(max_iter = 300))
    ))
    ll_per_df[i] <- tail(fit$history$log_likelihood_per_df, 1)
  }

  cat("\nNormalized log-likelihood per df by nu:", ll_per_df, "\n")
  # Should not diverge wildly as nu grows by 10x - both values should be
  # finite and within a similar order of magnitude of each other.
  expect_true(all(is.finite(ll_per_df)))
  expect_lt(abs(diff(ll_per_df)), 5)
})

# -- print/summary/coef/fitted consistency with the fitted object --------

test_that("CONSISTENCY: print/summary/coef/fitted stay consistent with fit internals", {
  truth <- .consistency_true_sigma_k1()
  dat <- .consistency_simulate(truth$Sigma, nu_per_study = 100, seed = 3600)

  fit <- suppressWarnings(suppressMessages(
    fit_covcomb(dat$S_list, dat$nu,
      scale_method = "estimate", alpha_normalization = "geometric",
      se_method = "plugin", n_factors = NULL,
      control = list(max_iter = 300)
    )
  ))

  expect_identical(coef(fit), fit$Sigma_hat)
  expect_identical(fitted(fit), fit$S_hat)

  print_out <- capture.output(print(fit))
  expect_true(any(grepl("Combined covariance matrix", print_out)))
  expect_true(any(grepl(sprintf("%d x %d", nrow(fit$Sigma_hat), ncol(fit$Sigma_hat)), print_out)))

  mean_diag_reported <- as.numeric(sub(".*Mean diagonal: ", "", print_out[grepl("Mean diagonal", print_out)][1]))
  expect_equal(mean_diag_reported, mean(diag(fit$Sigma_hat)), tolerance = 1e-3)

  summary_out <- capture.output(summary(fit))
  expect_true(any(grepl("Eigenvalue spectrum", summary_out)))
  expect_true(any(grepl("Standard errors", summary_out)))

  # print/summary must not error for a factor-model fit either, and should
  # reflect the same Sigma_hat via coef()/fitted().
  fit_fa <- suppressWarnings(suppressMessages(
    fit_fa_em(dat$S_list, dat$nu, k = 1, scale_method = "none")
  ))
  expect_identical(coef(fit_fa), fit_fa$Sigma_hat)
  expect_identical(fitted(fit_fa), fit_fa$S_hat)
  expect_silent(capture.output(print(fit_fa)))
})

# -- Scaling of estimated alpha with sample size --------------------------

test_that("CONSISTENCY: alpha_hat scale-factor recovery improves with sample size", {
  truth <- .consistency_true_sigma_k1()
  true_Sigma <- truth$Sigma
  var_names <- rownames(true_Sigma)
  true_alpha <- c(study1 = 0.7, study2 = 1.0, study3 = 1.6)

  .sim_alpha <- function(nu_per_study, seed) {
    set.seed(seed)
    patterns <- list(1:4, 3:6, c(1, 2, 5, 6))
    S_list <- list()
    nu_vec <- c()
    for (i in seq_along(patterns)) {
      idx <- patterns[[i]]
      id <- paste0("study", i)
      W <- rWishart(1, nu_per_study, true_alpha[[id]] * true_Sigma[idx, idx])[, , 1]
      S <- W / nu_per_study
      dimnames(S) <- list(var_names[idx], var_names[idx])
      S_list[[id]] <- S
      nu_vec[id] <- nu_per_study
    }
    list(S_list = S_list, nu = nu_vec)
  }

  # A single simulated dataset has high Monte Carlo variance in alpha_hat
  # (heterogeneous-scale estimation is noisy at small nu), so average the
  # recovery error over a handful of independent replicates at each sample
  # size to get a stable comparison while keeping the problem small.
  nu_levels <- c(30, 300)
  n_rep <- 6
  true_ratio <- true_alpha / true_alpha[["study2"]]
  ratio_errs <- numeric(length(nu_levels))

  for (i in seq_along(nu_levels)) {
    rep_errs <- numeric(n_rep)
    for (r in seq_len(n_rep)) {
      dat <- .sim_alpha(nu_levels[i], seed = 3700 + 100 * i + r)
      fit <- suppressWarnings(suppressMessages(
        fit_covcomb(dat$S_list, dat$nu,
          scale_method = "estimate", alpha_normalization = "arithmetic",
          se_method = "none", n_factors = NULL,
          control = list(max_iter = 300)
        )
      ))
      est_ratio <- fit$alpha_hat / fit$alpha_hat[["study2"]]
      rep_errs[r] <- max(abs(est_ratio - true_ratio))
    }
    ratio_errs[i] <- mean(rep_errs)
  }

  cat("\nAlpha ratio recovery mean error by nu:", ratio_errs, "\n")
  expect_lt(ratio_errs[2], ratio_errs[1])
  expect_lt(ratio_errs[2], 0.25)
})
