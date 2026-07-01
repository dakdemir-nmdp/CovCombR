## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)
suppressPackageStartupMessages(library(CovCombR))

set.seed(2025)
p_example  <- 6
var_names  <- paste0("V", seq_len(p_example))
true_Sigma <- diag(p_example)
dimnames(true_Sigma) <- list(var_names, var_names)

S1_example <- rWishart(1, 50, true_Sigma[1:4, 1:4])[,,1] / 50
dimnames(S1_example) <- list(var_names[1:4], var_names[1:4])

S2_example <- rWishart(1, 60, true_Sigma[3:6, 3:6])[,,1] / 60
dimnames(S2_example) <- list(var_names[3:6], var_names[3:6])

S_list <- list(Study1 = S1_example, Study2 = S2_example)
nu     <- c(Study1 = 50, Study2 = 60)

## ----n_factors----------------------------------------------------------------
# Default: auto-selection by BIC
fit_auto <- fit_covcomb(S_list, nu, se_method = "none")
fit_auto$n_factors
fit_auto$ic_table

# Force a specific number of factors
fit_1f <- fit_covcomb(S_list, nu, n_factors = 1, se_method = "none")
fit_2f <- fit_covcomb(S_list, nu, n_factors = 2, se_method = "none")

# Use the unconstrained free-Sigma model (requires all pairs observed)
fit_free <- fit_covcomb(S_list, nu, n_factors = NULL, se_method = "none")

# AIC instead of BIC for auto-selection
fit_aic <- fit_covcomb(S_list, nu,
                       n_factors = "auto",
                       control   = list(ic = "AIC"),
                       se_method = "none")
fit_aic$ic_table

## ----scale_method-------------------------------------------------------------
library(CovCombR)
set.seed(2025)

# Study 2 has inflated scale (alpha = 2)
S2_scaled <- rWishart(1, 60, 2 * true_Sigma[3:6, 3:6])[,,1] / 60
dimnames(S2_scaled) <- list(var_names[3:6], var_names[3:6])

fit_geom <- fit_covcomb(
  S_list = list(S1 = S1_example, S2 = S2_scaled),
  nu     = c(S1 = 50, S2 = 60),
  scale_method        = "estimate",
  alpha_normalization = "geometric"
)

fit_arith <- fit_covcomb(
  S_list = list(S1 = S1_example, S2 = S2_scaled),
  nu     = c(S1 = 50, S2 = 60),
  scale_method        = "estimate",
  alpha_normalization = "arithmetic"
)

# Geometric: Sigma_hat == S_hat (scale absorbed into alpha)
fit_geom$alpha_hat
all.equal(coef(fit_geom), fitted(fit_geom))

# Arithmetic: S_hat = mean(alpha) * Sigma_hat
fit_arith$alpha_hat
fit_arith$S_hat_scale

## ----init_sigma---------------------------------------------------------------
set.seed(2025)
var_names5 <- paste0("V", 1:5)
Sigma5     <- diag(5)
dimnames(Sigma5) <- list(var_names5, var_names5)

S1b <- rWishart(1, 100, Sigma5[1:3, 1:3])[,,1] / 100
S2b <- rWishart(1, 100, Sigma5[3:5, 3:5])[,,1] / 100
dimnames(S1b) <- list(var_names5[1:3], var_names5[1:3])
dimnames(S2b) <- list(var_names5[3:5], var_names5[3:5])

fit_avg <- fit_covcomb(
  list(S1 = S1b, S2 = S2b), c(S1 = 100, S2 = 100),
  init_sigma = "avg_padded", se_method = "none"
)
fit_id <- fit_covcomb(
  list(S1 = S1b, S2 = S2b), c(S1 = 100, S2 = 100),
  init_sigma = "identity", se_method = "none"
)

fit_avg$convergence$iterations
fit_id$convergence$iterations
max(abs(coef(fit_avg) - coef(fit_id)))

## ----control_params-----------------------------------------------------------
# Quick exploratory fit
fit_fast <- fit_covcomb(S_list, nu,
                        control   = list(max_iter = 100, tol = 1e-4),
                        se_method = "none")

# High-precision fit
fit_precise <- fit_covcomb(S_list, nu,
                           control   = list(max_iter = 2000, tol = 1e-9),
                           se_method = "none")

# Increased regularisation for near-singular problems
fit_stable <- fit_covcomb(S_list, nu,
                          control   = list(ridge = 1e-5, min_eigen = 1e-4),
                          se_method = "none")

## ----bootstrap_workflow, eval=FALSE-------------------------------------------
# Publication workflow: bootstrap SEs with reproducibility
fit_pub <- fit_covcomb(
  S_list = S_list,
  nu     = nu,
  se_method = "bootstrap",
  control   = list(
    bootstrap = list(
      B              = 200,
      seed           = 2025,
      progress       = FALSE,
      retain_samples = TRUE
    )
  )
)

# Access bootstrap samples (p x p x B array)
boot_samples <- fit_pub$bootstrap_samples

# Compute 95% quantile confidence intervals
CI_lower <- apply(boot_samples, 1:2, quantile, probs = 0.025)
CI_upper <- apply(boot_samples, 1:2, quantile, probs = 0.975)

# Hypothesis test: Is Sigma[1,2] significantly different from 0?
Sigma_hat <- coef(fit_pub)
Sigma_se  <- fit_pub$Sigma_se
z_stat    <- Sigma_hat[1, 2] / Sigma_se[1, 2]
p_value   <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
z_stat
p_value

## ----sem_workflow, eval=FALSE-------------------------------------------------
# SEM for large problems — always validate against bootstrap on a subset
fit_sem <- fit_covcomb(
  S_list    = S_list,
  nu        = nu,
  se_method = "sem",
  control   = list(sem = list(h = 1e-6, ridge = 1e-8))
)

fit_sem$sem$condition_number  # Should be < 1e6
fit_sem$sem$min_eigenvalue    # Should be > 0

# Validate on same data with bootstrap
fit_check <- fit_covcomb(S_list, nu, se_method = "bootstrap",
                         control = list(bootstrap = list(B = 100)))
cor(as.vector(fit_sem$Sigma_se),
    as.vector(fit_check$Sigma_se),
    use = "complete.obs")  # Should be > 0.95
