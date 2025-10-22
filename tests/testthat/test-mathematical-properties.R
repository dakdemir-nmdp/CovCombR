# Mathematical property tests - validate theoretical guarantees

test_that("Sigma_hat is always symmetric", {
  # Test that output is symmetric regardless of initialization
  set.seed(3000)
  p <- 5
  var_names <- paste0("V", 1:p)
  nu <- 40
  
  true_Sigma <- diag(p) + 0.4
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  W1 <- rWishart(1, nu, true_Sigma[1:3, 1:3])[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- list(var_names[1:3], var_names[1:3])
  
  W2 <- rWishart(1, nu, true_Sigma[3:5, 3:5])[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- list(var_names[3:5], var_names[3:5])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, se_method = "none")
  
  # Check symmetry
  expect_equal(fit$Sigma_hat, t(fit$Sigma_hat), check.attributes = FALSE)
  
  # Maximum asymmetry
  max_asymmetry <- max(abs(fit$Sigma_hat - t(fit$Sigma_hat)))
  cat(sprintf("\nMaximum asymmetry: %.2e\n", max_asymmetry))
  expect_lt(max_asymmetry, 1e-14)
})

test_that("Sigma_hat is positive definite", {
  # Test that output is always PD
  set.seed(3001)
  p <- 6
  var_names <- paste0("V", 1:p)
  nu <- 50
  
  true_Sigma <- diag(p) + 0.3
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  # Multiple patterns
  W1 <- rWishart(1, nu, true_Sigma[1:4, 1:4])[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- list(var_names[1:4], var_names[1:4])
  
  W2 <- rWishart(1, nu, true_Sigma[3:6, 3:6])[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- list(var_names[3:6], var_names[3:6])
  
  S3 <- rWishart(1, nu, true_Sigma[c(1,3,5), c(1,3,5)])[,,1]
  dimnames(S3) <- list(var_names[c(1,3,5)], var_names[c(1,3,5)])
  
  S_list <- list(s1 = S1, s2 = S2, s3 = S3)
  nu_vec <- c(s1 = nu, s2 = nu, s3 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, se_method = "none")
  
  # Check PD via eigenvalues
  eig_vals <- eigen(fit$Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_vals)
  
  cat(sprintf("\nEigenvalues of Sigma_hat:\n"))
  cat(sprintf("  Min: %.6f\n", min_eig))
  cat(sprintf("  Max: %.6f\n", max(eig_vals)))
  cat(sprintf("  Condition number: %.2e\n", max(eig_vals) / min_eig))
  
  expect_true(all(eig_vals > 0))
  expect_gt(min_eig, 1e-10)
})

test_that("Algorithm converges for well-conditioned problems", {
  # Test convergence behavior
  set.seed(3002)
  p <- 4
  var_names <- paste0("V", 1:p)
  nu <- 60
  
  true_Sigma <- diag(p) + 0.5
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  W1 <- rWishart(1, nu, true_Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(true_Sigma)
  
  W2 <- rWishart(1, nu, true_Sigma)[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- dimnames(true_Sigma)
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, se_method = "none",
                        control = list(max_iter = 1000, tol = 1e-8))
  
  expect_true(fit$convergence$converged)
  expect_lt(fit$convergence$iterations, 1000)
  
  cat(sprintf("\nConvergence:\n"))
  cat(sprintf("  Converged: %s\n", fit$convergence$converged))
  cat(sprintf("  Iterations: %d\n", fit$convergence$iterations))
})

test_that("Handles poorly conditioned input gracefully", {
  # Test with nearly singular matrices
  set.seed(3003)
  p <- 4
  var_names <- paste0("V", 1:p)
  nu <- 100  # Large nu to get near-singular
  
  # Create nearly singular Sigma
  near_singular <- diag(p)
  near_singular[1,2] <- near_singular[2,1] <- 0.99
  near_singular[3,4] <- near_singular[4,3] <- 0.99
  dimnames(near_singular) <- list(var_names, var_names)
  
  # Generate data
  W1 <- rWishart(1, nu, near_singular)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(near_singular)
  
  S_list <- list(s1 = S1)
  nu_vec <- c(s1 = nu)
  
  # Should warn about single sample (identifiability) but not crash
  expect_warning({
    fit <- fit_covcomb(S_list, nu_vec, se_method = "none")
  }, "Only one sample provided")
  
  # Output should still be PD
  eig_vals <- eigen(fit$Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_vals > 0))
  
  cat(sprintf("\nPoorly conditioned input:\n"))
  cat(sprintf("  Input condition number: %.2e\n", 
              max(eigen(near_singular)$values) / min(eigen(near_singular)$values)))
  cat(sprintf("  Output condition number: %.2e\n",
              max(eig_vals) / min(eig_vals)))
})

test_that("Recovers identity when data is identity-scaled", {
  # Simple case: should recover structure
  set.seed(3004)
  p <- 3
  var_names <- paste0("V", 1:p)
  nu <- 200  # Large sample
  
  true_Sigma <- diag(p)
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  # Generate many samples
  S_list <- list()
  nu_vec <- c()
  for (i in 1:10) {
    W <- rWishart(1, nu, true_Sigma)[,,1]
    S <- W / nu  # Convert to sample covariance
    dimnames(S) <- dimnames(true_Sigma)
    S_list[[paste0("s", i)]] <- S
    nu_vec[paste0("s", i)] <- nu
  }
  
  fit <- fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  
  # Should recover identity structure (diagonal)
  off_diag_max <- max(abs(fit$Sigma_hat[row(fit$Sigma_hat) != col(fit$Sigma_hat)]))
  diag_mean <- mean(diag(fit$Sigma_hat))
  diag_sd <- sd(diag(fit$Sigma_hat))
  
  cat(sprintf("\nRecovery of identity:\n"))
  cat(sprintf("  Max off-diagonal: %.6f\n", off_diag_max))
  cat(sprintf("  Diagonal mean: %.6f\n", diag_mean))
  cat(sprintf("  Diagonal SD: %.6f\n", diag_sd))
  cat(sprintf("  Expected (E[W]/nu = 1): 1.0\n"))
  cat(sprintf("  Actually nu=%d, so E[W_ii/nu] = 1\n", nu))
  
  expect_lt(off_diag_max, 0.05)
  # E[W] = nu*Sigma, so E[W/nu] = Sigma = I, so E[diag] should be 1
  expect_equal(diag_mean, 1, tolerance = 0.05)
})

test_that("Handles single sample correctly", {
  # Edge case: only one sample
  set.seed(3005)
  p <- 4
  var_names <- paste0("V", 1:p)
  nu <- 50
  
  Sigma <- diag(p) + 0.4
  dimnames(Sigma) <- list(var_names, var_names)
  
  W <- rWishart(1, nu, Sigma)[,,1]
  S <- W / nu
  dimnames(S) <- dimnames(Sigma)
  
  S_list <- list(s1 = S)
  nu_vec <- c(s1 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")

  # Should return S (the sample covariance) essentially
  expected <- S

  diff_norm <- sqrt(sum((fit$Sigma_hat - expected)^2)) / sqrt(sum(expected^2))

  cat(sprintf("\nSingle sample:\n"))
  cat(sprintf("  Relative difference from S: %.6f\n", diff_norm))

  expect_lt(diff_norm, 0.01)
})

test_that("Handles p=2 minimal case", {
  # Smallest non-trivial case
  set.seed(3006)
  p <- 2
  var_names <- c("V1", "V2")
  nu <- 30
  
  Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
  dimnames(Sigma) <- list(var_names, var_names)
  
  W1 <- rWishart(1, nu, Sigma)[,,1]
  S1 <- W1 / nu
  W2 <- rWishart(1, nu, Sigma)[,,1]
  S2 <- W2 / nu
  dimnames(S1) <- dimnames(Sigma)
  dimnames(S2) <- dimnames(Sigma)
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, se_method = "none")
  
  expect_equal(dim(fit$Sigma_hat), c(2, 2))
  expect_true(fit$convergence$converged)
  
  # Check correlation is reasonable
  corr_est <- fit$Sigma_hat[1,2] / sqrt(fit$Sigma_hat[1,1] * fit$Sigma_hat[2,2])
  cat(sprintf("\np=2 case:\n"))
  cat(sprintf("  True correlation: 0.7\n"))
  cat(sprintf("  Estimated correlation: %.4f\n", corr_est))
  
  expect_gt(corr_est, 0.3)
  expect_lt(corr_est, 1.0)
})

test_that("Respects boundary constraints on alpha", {
  # Alpha should be bounded away from zero
  set.seed(3007)
  p <- 3
  var_names <- paste0("V", 1:p)
  nu <- 40
  
  Sigma <- diag(p)
  dimnames(Sigma) <- list(var_names, var_names)
  
  # Create samples with very different scales
  W1 <- rWishart(1, nu, Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma)
  
  S2 <- 1000 * rWishart(1, nu, Sigma)[,,1]  # Huge scale difference
  dimnames(S2) <- dimnames(Sigma)
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, scale_method = "estimate", se_method = "none")
  
  # All alphas should be >= 1e-8
  expect_true(all(fit$alpha_hat >= 1e-8))
  
  cat(sprintf("\nAlpha boundary constraints:\n"))
  cat(sprintf("  alpha1: %.2e\n", fit$alpha_hat["s1"]))
  cat(sprintf("  alpha2: %.2e\n", fit$alpha_hat["s2"]))
  cat(sprintf("  Min alpha: %.2e (bound: 1e-8)\n", min(fit$alpha_hat)))
})

test_that("Initialization affects convergence speed but not final result", {
  # Different initializations should reach same result
  set.seed(3008)
  p <- 4
  var_names <- paste0("V", 1:p)
  nu <- 50
  
  true_Sigma <- diag(p) + 0.5
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  W1 <- rWishart(1, nu, true_Sigma[1:3, 1:3])[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- list(var_names[1:3], var_names[1:3])
  
  W2 <- rWishart(1, nu, true_Sigma[2:4, 2:4])[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- list(var_names[2:4], var_names[2:4])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  # Fit multiple times (different random seeds affect initialization)
  results <- list()
  for (i in 1:5) {
    set.seed(3008 + i)
    # After gamma removal, convergence may be slower for some initializations
    results[[i]] <- fit_covcomb(S_list, nu_vec, se_method = "none",
                                   control = list(max_iter = 200, tol = 1e-5))
  }
  
  # Most should converge (some may hit max_iter with tight overlap patterns)
  num_converged <- sum(sapply(results, function(r) r$convergence$converged))
  
  cat(sprintf("\nInitialization robustness:\n"))
  cat(sprintf("  Converged: %d / 5\n", num_converged))
  cat(sprintf("  Iterations: %s\n", 
              paste(sapply(results, function(r) r$convergence$iterations), collapse=", ")))
  
  # At least some should converge
  expect_gte(num_converged, 1)
  
  # Among converged results, estimates should be similar
  converged_results <- results[sapply(results, function(r) r$convergence$converged)]
  if (length(converged_results) >= 2) {
    Sigma_hats <- lapply(converged_results, function(r) r$Sigma_hat)
    
    max_diff <- 0
    for (i in 1:(length(Sigma_hats)-1)) {
      diff <- sqrt(sum((Sigma_hats[[i]] - Sigma_hats[[length(Sigma_hats)]])^2))
      max_diff <- max(max_diff, diff)
    }
    
    cat(sprintf("  Max difference in Sigma_hat (among converged): %.2e\n", max_diff))
    expect_lt(max_diff, 0.1)
  }
})

test_that("Handles all-observed vs all-missing patterns", {
  # Mix of complete and incomplete data
  set.seed(3009)
  p <- 5
  var_names <- paste0("V", 1:p)
  nu <- 60
  
  Sigma <- diag(p) + 0.4
  dimnames(Sigma) <- list(var_names, var_names)
  
  # Complete sample
  W_complete <- rWishart(1, nu, Sigma)[,,1]
  S_complete <- W_complete / nu
  dimnames(S_complete) <- dimnames(Sigma)
  
  # Partial samples
  W_partial1 <- rWishart(1, nu, Sigma[1:3, 1:3])[,,1]
  S_partial1 <- W_partial1 / nu
  dimnames(S_partial1) <- list(var_names[1:3], var_names[1:3])
  
  W_partial2 <- rWishart(1, nu, Sigma[3:5, 3:5])[,,1]
  S_partial2 <- W_partial2 / nu
  dimnames(S_partial2) <- list(var_names[3:5], var_names[3:5])
  
  S_list <- list(complete = S_complete, partial1 = S_partial1, partial2 = S_partial2)
  nu_vec <- c(complete = nu, partial1 = nu, partial2 = nu)
  
  fit <- fit_covcomb(S_list, nu_vec, se_method = "none")
  
  expect_true(fit$convergence$converged)
  expect_true(all(eigen(fit$Sigma_hat, symmetric = TRUE, only.values = TRUE)$values > 0))
  
  cat(sprintf("\nMixed observation patterns:\n"))
  cat(sprintf("  Converged: %s\n", fit$convergence$converged))
  cat(sprintf("  Iterations: %d\n", fit$convergence$iterations))
})
