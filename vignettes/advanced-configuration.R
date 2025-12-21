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

# Load package
suppressPackageStartupMessages(library(CovCombR))

# Create reusable example data for demonstrations throughout
set.seed(2025)
p_example <- 6
var_names_example <- paste0("V", 1:p_example)
true_Sigma_example <- diag(p_example)
dimnames(true_Sigma_example) <- list(var_names_example, var_names_example)

# Create overlapping study data
S1_example <- rWishart(1, 50, true_Sigma_example[1:4, 1:4])[,,1] / 50
dimnames(S1_example) <- list(var_names_example[1:4], var_names_example[1:4])

S2_example <- rWishart(1, 60, true_Sigma_example[3:6, 3:6])[,,1] / 60
dimnames(S2_example) <- list(var_names_example[3:6], var_names_example[3:6])

# Standard S_list and nu for examples
S_list <- list(Study1 = S1_example, Study2 = S2_example)
nu <- c(Study1 = 50, Study2 = 60)

## ----eval=TRUE----------------------------------------------------------------
library(CovCombR)

# Simulate data with heterogeneous scales
set.seed(2025)
p <- 6
var_names <- paste0("V", 1:p)
true_Sigma <- diag(p)
dimnames(true_Sigma) <- list(var_names, var_names)

# Study 1: standard scale (α = 1)
S1 <- rWishart(1, 50, true_Sigma[1:4, 1:4])[,,1] / 50
dimnames(S1) <- list(var_names[1:4], var_names[1:4])

# Study 2: inflated scale (α = 2)
S2 <- rWishart(1, 60, 2 * true_Sigma[3:6, 3:6])[,,1] / 60
dimnames(S2) <- list(var_names[3:6], var_names[3:6])

# Fit with geometric normalization (default)
fit_geom <- fit_covcomb(
  S_list = list(S1 = S1, S2 = S2),
  nu = c(S1 = 50, S2 = 60),
  scale_method = "estimate",
  alpha_normalization = "geometric"
)

# Fit with arithmetic normalization
fit_arith <- fit_covcomb(
  S_list = list(S1 = S1, S2 = S2),
  nu = c(S1 = 50, S2 = 60),
  scale_method = "estimate",
  alpha_normalization = "arithmetic"
)

# Compare
cat("Geometric normalization:\n")
cat("  Alpha:", fit_geom$alpha_hat, "\n")
cat("  Sigma_hat[1,1]:", coef(fit_geom)[1,1], "\n")
cat("  S_hat[1,1]:", fitted(fit_geom)[1,1], "\n")
cat("  Identical?", all.equal(coef(fit_geom), fitted(fit_geom)), "\n\n")

cat("Arithmetic normalization:\n")
cat("  Alpha:", fit_arith$alpha_hat, "\n")
cat("  Sigma_hat[1,1]:", coef(fit_arith)[1,1], "\n")
cat("  S_hat[1,1]:", fitted(fit_arith)[1,1], "\n")
cat("  Identical?", all.equal(coef(fit_arith), fitted(fit_arith)), "\n")

## ----eval=TRUE----------------------------------------------------------------
# Use the example data from setup
fit <- fit_covcomb(S_list, nu, scale_method = "estimate")
cat("Alpha parameters:\n")
print(fit$alpha_hat)

## ----eval=FALSE---------------------------------------------------------------
# # Example workflow (hypothetical - requires actual genotype data)
# # GRM from 50K SNP array
# GRM_50K <- compute_grm(genotypes_50K)  # Your GRM computation function
# 
# # GRM from 500K SNP array
# GRM_500K <- compute_grm(genotypes_500K)
# 
# # Combine with scale estimation
# fit <- fit_covcomb(
#   S_list = list(Array50K = GRM_50K, Array500K = GRM_500K),
#   nu = c(Array50K = 50000, Array500K = 500000),
#   scale_method = "estimate"
# )
# 
# # Check relative scales
# print(fit$alpha_hat)
# # Expected: alpha for 500K might be larger due to finer resolution

## ----eval=TRUE----------------------------------------------------------------
# Difficult case: chain topology with small p
set.seed(2025)
var_names <- paste0("V", 1:5)
true_Sigma <- diag(5)
dimnames(true_Sigma) <- list(var_names, var_names)

S1 <- rWishart(1, 100, true_Sigma[1:3, 1:3])[,,1] / 100
dimnames(S1) <- list(var_names[1:3], var_names[1:3])

S2 <- rWishart(1, 100, true_Sigma[3:5, 3:5])[,,1] / 100
dimnames(S2) <- list(var_names[3:5], var_names[3:5])

# Default initialization
fit_default <- fit_covcomb(
  S_list = list(S1 = S1, S2 = S2),
  nu = c(S1 = 100, S2 = 100),
  init_sigma = "avg_padded"
)

# Identity initialization
fit_identity <- fit_covcomb(
  S_list = list(S1 = S1, S2 = S2),
  nu = c(S1 = 100, S2 = 100),
  init_sigma = "identity"
)

# Compare convergence
cat("Avg padded iterations:", fit_default$convergence$iterations, "\n")
cat("Identity iterations:", fit_identity$convergence$iterations, "\n")

# Compare estimates
cat("Max difference:", max(abs(coef(fit_default) - coef(fit_identity))), "\n")

## ----eval=TRUE----------------------------------------------------------------
# Quick and dirty fit
fit_fast <- fit_covcomb(S_list, nu,
                        control = list(max_iter = 100, tol = 1e-4))

# High precision fit
fit_precise <- fit_covcomb(S_list, nu,
                           control = list(max_iter = 2000, tol = 1e-8))

## ----eval=TRUE----------------------------------------------------------------
# Problematic dataset with near-singularity
fit_stable <- fit_covcomb(S_list, nu,
                          control = list(ridge = 1e-5, min_eigen = 1e-4))

## ----eval=FALSE---------------------------------------------------------------
# # Note: In practice use B=200 for publication; here we show workflow
# fit_pub <- fit_covcomb(
#   S_list = S_list,
#   nu = nu,
#   se_method = "bootstrap",
#   control = list(
#     bootstrap = list(
#       B = 200,              # 200 replicates for publication
#       seed = 2025,          # Reproducibility
#       parallel = FALSE,     # Serial execution (safer)
#       retain_samples = TRUE # Keep all bootstrap estimates
#     )
#   )
# )
# 
# # Access bootstrap samples
# boot_samples <- fit_pub$bootstrap_samples  # Array: p x p x B
# 
# # Compute 95% confidence intervals (quantile method)
# Sigma_hat <- coef(fit_pub)
# CI_lower <- apply(boot_samples, 1:2, quantile, probs = 0.025)
# CI_upper <- apply(boot_samples, 1:2, quantile, probs = 0.975)
# 
# # Hypothesis test: Is Sigma[1,2] significantly different from 0?
# cat("Estimate:", Sigma_hat[1,2], "\n")
# cat("95% CI: [", CI_lower[1,2], ",", CI_upper[1,2], "]\n")

## ----eval=FALSE---------------------------------------------------------------
# # Bootstrap convergence monitoring (use smaller B for speed)
# fit <- fit_covcomb(S_list, nu, se_method = "bootstrap",
#                   control = list(bootstrap = list(B = 200, seed = 2025)))
# 
# # Check convergence status (stored internally)
# # If many bootstrap replicates fail, SEs may be unreliable
# # Increase max_iter or relax tol if this occurs

## ----eval=FALSE---------------------------------------------------------------
# # Large problem where bootstrap is slow (hypothetical workflow)
# set.seed(2025)
# p <- 50
# var_names <- paste0("V", 1:p)
# 
# # Create large study data (simulation code omitted for brevity)
# # S_list_large <- create_large_study_data(p = 50, K = 3)
# # nu_large <- c(Study1 = 500, Study2 = 600, Study3 = 550)
# 
# # Fit with SEM
# fit_sem <- fit_covcomb(
#   S_list = S_list_large,
#   nu = nu_large,
#   se_method = "sem",
#   control = list(
#     sem = list(h = 1e-6, ridge = 1e-6)
#   )
# )
# 
# # CRITICAL: Validate SEM against bootstrap on a subset
# p_subset <- 10
# indices <- 1:p_subset
# 
# S_list_subset <- lapply(S_list_large, function(S) {
#   S[indices, indices]
# })
# 
# fit_boot_subset <- fit_covcomb(
#   S_list = S_list_subset,
#   nu = nu_large,
#   se_method = "bootstrap",
#   control = list(bootstrap = list(B = 200, seed = 2025))
# )
# 
# fit_sem_subset <- fit_covcomb(
#   S_list = S_list_subset,
#   nu = nu_large,
#   se_method = "sem",
#   control = list(sem = list(h = 1e-6))
# )
# 
# # Compare SEs
# se_boot <- fit_boot_subset$Sigma_se
# se_sem <- fit_sem_subset$Sigma_se
# 
# cat("Mean SE (bootstrap):", mean(se_boot, na.rm = TRUE), "\n")
# cat("Mean SE (SEM):", mean(se_sem, na.rm = TRUE), "\n")
# cat("Correlation:", cor(as.vector(se_boot), as.vector(se_sem),
#                         use = "complete.obs"), "\n")

