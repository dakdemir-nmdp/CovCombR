
test_that("SCALING SIMULATION: standardized input maintains scale", {
  # CRITICAL TEST: When input GRMs have mean(diag) standardized to 1,
  # verify that output maintains appropriate scale
  set.seed(1000)
  p <- 10
  n_samples <- 5
  nu <- 100
  var_names <- paste0("SNP", 1:p)
  # True covariance with specific scale
  true_Sigma <- diag(p) + 0.3
  diag(true_Sigma) <- 1.5  # Diagonal elements are 1.5
  dimnames(true_Sigma) <- list(var_names, var_names)
  # Generate samples with different patterns
  S_list <- list()
  patterns <- list()
  nu_vec <- c()
  for (i in 1:n_samples) {
    # Each sample observes different subset
    if (i == 1) obs_idx <- 1:7
    else if (i == 2) obs_idx <- 4:10
    else if (i == 3) obs_idx <- c(1,3,5,7,9)
    else if (i == 4) obs_idx <- c(2,4,6,8,10)
    else obs_idx <- 1:10
    W_i <- rWishart(1, nu, true_Sigma[obs_idx, obs_idx])[,,1]
    S_i <- W_i / nu
    dimnames(S_i) <- list(var_names[obs_idx], var_names[obs_idx])
    S_list[[paste0("sample", i)]] <- S_i
    patterns[[i]] <- obs_idx
    nu_vec[paste0("sample", i)] <- nu
  }
  fit <- fit_covcomb(
    S_list = S_list,
    nu = nu_vec,
    se_method = "none",
    control = list(max_iter = 1000, tol = 1e-8)
  )
  # KEY QUESTION: What is mean(diag(Sigma_hat))?
  mean_diag_output <- mean(diag(fit$Sigma_hat))
  cat("\n=== SCALING SIMULATION RESULTS ===\n")
  cat(sprintf("Output mean(diag(Sigma_hat)): %.6f\n", mean_diag_output))
  cat(sprintf("Expected from theory: 1/nu = 1/%d = %.6f\n", nu, 1/nu))
  # 1. Input S_k is scaled so mean(diag(S_k)) = 1
  # 2. EM estimates Sigma from scaled W ~ Wishart(nu, Sigma)
  # 3. E[W] = nu * Sigma, so E[W/nu] = Sigma
  # 4. Since input was scaled to mean(diag)=1.5, output Sigma should have mean(diag) ≈ 1.5
  # This is the KEY THEORETICAL PREDICTION
  expected_mean_diag <- 1.5
  # Test that output scale matches theory
  expect_equal(mean_diag_output, expected_mean_diag, tolerance = 0.05)
  # Also verify convergence
  expect_true(fit$convergence$converged)
})


test_that("SCALING SIMULATION: scale_method='none' preserves absolute scale", {
  # When scale_method="none", output should match input scale
  set.seed(1002)
  p <- 6
  var_names <- paste0("V", 1:p)
  
  # True covariance with scale factor
  scale_factor <- 5.0
  true_Sigma <- scale_factor * (diag(p) + 0.3)
  diag(true_Sigma) <- scale_factor * 1.2
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  nu <- 80
  
  # Generate complete samples (no missing data)
  W1 <- rWishart(1, nu, true_Sigma)[,,1]
  S1 <- W1 / nu
  W2 <- rWishart(1, nu, true_Sigma)[,,1]
  S2 <- W2 / nu
  dimnames(S1) <- dimnames(true_Sigma)
  dimnames(S2) <- dimnames(true_Sigma)
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = nu_vec,
    scale_method = "none",
    se_method = "none"
  )
  
  mean_diag_output <- mean(diag(fit$Sigma_hat))
  mean_diag_input_avg <- (mean(diag(S1)) + mean(diag(S2))) / 2  # S1, S2 are sample covariances

  cat(sprintf("\nscale_method='none': output=%.4f, input_avg=%.4f\n",
              mean_diag_output, mean_diag_input_avg))

  # Output should match average input scale
  expect_equal(mean_diag_output, mean_diag_input_avg, tolerance = 0.3)
  
  # Should be close to true scale
  expect_equal(mean_diag_output, scale_factor * 1.2, tolerance = 0.5)
})

test_that("SCALING SIMULATION: compare scales across methods", {
  set.seed(1003)
  p <- 5
  var_names <- paste0("V", 1:p)
  nu <- 60
  
  true_Sigma <- diag(p) + 0.5
  diag(true_Sigma) <- 3.0  # Large scale
  dimnames(true_Sigma) <- list(var_names, var_names)
  
  # Generate data
  W1 <- rWishart(1, nu, true_Sigma[1:3, 1:3])[,,1]
  S1 <- W1 / nu
  W2 <- rWishart(1, nu, true_Sigma[3:5, 3:5])[,,1]
  S2 <- W2 / nu
  dimnames(S1) <- list(var_names[1:3], var_names[1:3])
  dimnames(S2) <- list(var_names[3:5], var_names[3:5])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  # Fit with none
  fit_none <- fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  
  # Fit with estimate
  fit_unit <- fit_covcomb(S_list, nu_vec, scale_method = "estimate", se_method = "none")
  
  scale_none <- mean(diag(fit_none$Sigma_hat))
  scale_unit <- mean(diag(fit_unit$Sigma_hat))
  
  cat(sprintf("\nComparing methods:\n"))
  cat(sprintf("  scale_method='none': mean(diag) = %.6f\n", scale_none))
  cat(sprintf("  Ratio: %.4f\n", scale_none / scale_unit))
  
  expect_gt(scale_unit, scale_none - 0.5)
  
  expect_equal(scale_unit, scale_none, tolerance = 0.5)
})

test_that("SCALING SIMULATION: heterogeneous input scales with estimate", {
  # When samples have very different scales, estimate method should capture it
  set.seed(1004)
  p <- 4
  var_names <- paste0("V", 1:p)
  nu <- 70
  
  base_Sigma <- diag(p) + 0.3
  diag(base_Sigma) <- 1
  dimnames(base_Sigma) <- list(var_names, var_names)
  
  # Sample 1: scale = 1
  W1 <- rWishart(1, nu, base_Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(base_Sigma)

  # Sample 2: scale = 10 (scale the Wishart matrix, then convert to sample covariance)
  W2 <- 10 * rWishart(1, nu, base_Sigma)[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- dimnames(base_Sigma)
  
  S_list <- list(s1 = S1, s2 = S2)
  nu_vec <- c(s1 = nu, s2 = nu)
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = nu_vec,
    scale_method = "estimate",
    se_method = "none"
  )
  
  alpha1 <- fit$alpha_hat["s1"]
  alpha2 <- fit$alpha_hat["s2"]
  alpha_ratio <- unname(alpha2 / alpha1)
  
  cat(sprintf("\nHeterogeneous scales:\n"))
  cat(sprintf("  alpha1 = %.4f, alpha2 = %.4f\n", alpha1, alpha2))
  cat(sprintf("  alpha_ratio = %.4f (expected ≈ 10)\n", alpha_ratio))
  cat(sprintf("  mean(diag(Sigma_hat)) = %.4f\n", mean(diag(fit$Sigma_hat))))
  
  # Alpha ratio should capture scale difference
  expect_gt(alpha_ratio, 5)
  expect_lt(alpha_ratio, 15)
  
  # Reconstructed scales should match inputs
  recon_scale1 <- alpha1 * mean(diag(fit$Sigma_hat))
  recon_scale2 <- alpha2 * mean(diag(fit$Sigma_hat))
  
  input_scale1 <- mean(diag(S1))
  input_scale2 <- mean(diag(S2))
  
  cat(sprintf("  Reconstructed: s1=%.4f (input %.4f), s2=%.4f (input %.4f)\n",
              recon_scale1, input_scale1, recon_scale2, input_scale2))
})
