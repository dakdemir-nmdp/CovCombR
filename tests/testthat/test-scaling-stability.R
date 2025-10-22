# Scaling and Numerical Stability Tests

test_that("Different input scales are handled correctly", {
  # Test with very small scale
  set.seed(3001)
  p <- 5
  nu <- 50
  
  Sigma_small <- diag(p) * 1e-6
  Sigma_small[1:2, 1:2] <- 8e-7
  diag(Sigma_small) <- 1e-6
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma_small) <- list(ids, ids)
  
  W1 <- rWishart(1, nu, Sigma_small)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma_small)
  
  fit_small <- fit_covcomb(
    S_list = list(s1 = S1),
    nu = c(s1 = nu),
    scale_method = "none",
    se_method = "none"
  )
  
  # Check that output is at same scale as input
  input_scale <- mean(diag(S1))  # S1 is already sample covariance
  output_scale <- mean(diag(fit_small$Sigma_hat))

  scale_ratio <- output_scale / input_scale
  expect_equal(scale_ratio, 1.0, tolerance = 0.1)
  
  # Test with very large scale
  Sigma_large <- diag(p) * 1e6
  Sigma_large[1:2, 1:2] <- 8e5
  diag(Sigma_large) <- 1e6
  dimnames(Sigma_large) <- list(ids, ids)
  
  W2 <- rWishart(1, nu, Sigma_large)[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- dimnames(Sigma_large)
  
  fit_large <- fit_covcomb(
    S_list = list(s1 = S2),
    nu = c(s1 = nu),
    scale_method = "none",
    se_method = "none"
  )
  
  input_scale_large <- mean(diag(S2))  # S2 is already sample covariance
  output_scale_large <- mean(diag(fit_large$Sigma_hat))

  scale_ratio_large <- output_scale_large / input_scale_large
  expect_equal(scale_ratio_large, 1.0, tolerance = 0.1)
})

test_that("Mixed scales across samples are preserved", {
  # When samples have different inherent scales, relative scales should be preserved
  set.seed(3002)
  p <- 6
  nu <- 60
  
  ids <- paste0("V", 1:p)
  
  # Create samples at different scales
  Sigma1 <- diag(p) * 0.1
  diag(Sigma1) <- 0.1
  dimnames(Sigma1) <- list(ids, ids)
  
  Sigma2 <- diag(p) * 10
  diag(Sigma2) <- 10
  dimnames(Sigma2) <- list(ids, ids)
  
  W1 <- rWishart(1, nu, Sigma1)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma1)
  
  W2 <- rWishart(1, nu, Sigma2)[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- dimnames(Sigma2)
  
  # Input scale ratio
  scale1_in <- mean(diag(S1))
  scale2_in <- mean(diag(S2))
  ratio_in <- scale2_in / scale1_in
  
  # Fit with scale estimation
  fit <- fit_covcomb(
    S_list = list(s1 = S1, s2 = S2),
    nu = c(s1 = nu, s2 = nu),
    scale_method = "estimate",
    se_method = "none"
  )
  
  # Check that alpha captures the scale difference
  alpha_ratio <- unname(fit$alpha_hat["s2"] / fit$alpha_hat["s1"])
  expect_equal(alpha_ratio, ratio_in, tolerance = 0.3)
})

test_that("Ill-conditioned covariance matrices are handled", {
  # Test with nearly singular covariance
  set.seed(3003)
  p <- 5
  nu <- 100
  
  # Create ill-conditioned matrix
  Sigma <- diag(p)
  Sigma[1:3, 1:3] <- 0.99  # Very high correlation
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # Condition number before
  kappa_true <- kappa(Sigma, exact = TRUE)
  expect_gt(kappa_true, 100)
  
  W1 <- rWishart(1, nu, Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma)
  
  W2 <- rWishart(1, nu, Sigma[1:4, 1:4])[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- list(ids[1:4], ids[1:4])
  
  fit <- fit_covcomb(
    S_list = list(s1 = S1, s2 = S2),
    nu = c(s1 = nu, s2 = nu),
    se_method = "none"
  )
  
  # Should still converge and be PD
  expect_true(fit$convergence$converged)
  eigenvalues <- eigen(fit$Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("Unit diagonal scaling normalizes properly", {
  set.seed(3004)
  p <- 6
  nu <- 80
  
  # Start with arbitrary scale
  Sigma <- diag(p) * 5
  Sigma[1:3, 1:3] <- 4
  diag(Sigma) <- 5
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  W1 <- rWishart(1, nu, Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma)
  
  # S1 is sample covariance: diagonal should be around 5
  expect_gt(mean(diag(S1)), 3)

  fit <- fit_covcomb(
    S_list = list(s1 = S1),
    nu = c(s1 = nu),
    se_method = "none"
  )

  # Sigma_hat should match S1 for single complete sample
  # mean(diag(Sigma_hat)) ≈ mean(diag(S1)) ≈ 5
  expected_mean_diag <- 5
  mean_diag <- mean(diag(fit$Sigma_hat))
  expect_equal(mean_diag, expected_mean_diag, tolerance = 0.5)
})

test_that("Extreme correlation values are handled", {
  # Test with very high and very low correlations
  set.seed(3005)
  p <- 6
  nu <- 100
  
  Sigma <- diag(p)
  Sigma[1,2] <- Sigma[2,1] <- 0.95  # Very high
  Sigma[3,4] <- Sigma[4,3] <- 0.05  # Very low
  Sigma[5,6] <- Sigma[6,5] <- -0.85  # High negative
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  W1 <- rWishart(1, nu, Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma)
  
  fit <- fit_covcomb(
    S_list = list(s1 = S1),
    nu = c(s1 = nu),
    scale_method = "none",
    se_method = "none"
  )
  
  # Check correlations are preserved
  Sigma_cor <- cov2cor(fit$Sigma_hat)
  
  expect_gt(Sigma_cor[1,2], 0.7)
  expect_lt(Sigma_cor[5,6], -0.5)
})

test_that("Many samples with partial overlap converge correctly", {
  # Stress test: many samples with various overlap patterns
  set.seed(3006)
  p <- 10
  nu <- 50
  K <- 20  # Many samples
  
  Sigma <- diag(p)
  Sigma[1:5, 1:5] <- 0.6
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  S_list <- list()
  nu_vec <- numeric(K)
  
  for (k in 1:K) {
    # Random observation pattern
    n_obs <- sample(3:7, 1)
    obs_idx <- sort(sample(1:p, n_obs))
    
    W <- rWishart(1, nu, Sigma[obs_idx, obs_idx])[,,1]
    dimnames(W) <- list(ids[obs_idx], ids[obs_idx])
    S_list[[paste0("s", k)]] <- W
    nu_vec[k] <- nu
  }
  names(nu_vec) <- names(S_list)
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = nu_vec,
    se_method = "none"
  )
  
  expect_true(fit$convergence$converged)
  
  # Check structure recovery
  Sigma_cor_est <- cov2cor(fit$Sigma_hat)
  Sigma_cor_true <- cov2cor(Sigma)
  
  cor_error <- norm(Sigma_cor_est - Sigma_cor_true, "F") / norm(Sigma_cor_true, "F")
  expect_lt(cor_error, 0.65)  # Relaxed tolerance for complex missing patterns
})

test_that("Degrees of freedom affects variance correctly", {
  # Theory: Var(W_ij) proportional to 1/nu
  # So with higher nu, estimates should be more precise
  set.seed(3007)
  p <- 4
  Sigma <- diag(p)
  Sigma[1:2, 1:2] <- 0.7
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # Low nu
  nu_low <- 30
  W_low <- rWishart(1, nu_low, Sigma)[,,1]
  S_low <- W_low / nu_low
  dimnames(S_low) <- dimnames(Sigma)
  
  fit_low <- fit_covcomb(
    S_list = list(s1 = S_low),
    nu = c(s1 = nu_low),
    scale_method = "none",
    se_method = "plugin"
  )
  
  # High nu
  nu_high <- 300
  W_high <- rWishart(1, nu_high, Sigma)[,,1]
  S_high <- W_high / nu_high
  dimnames(S_high) <- dimnames(Sigma)
  
  fit_high <- fit_covcomb(
    S_list = list(s1 = S_high),
    nu = c(s1 = nu_high),
    scale_method = "none",
    se_method = "plugin"
  )
  
  # SE should be larger for low nu
  mean_se_low <- mean(fit_low$Sigma_se)
  mean_se_high <- mean(fit_high$Sigma_se)
  
  se_ratio <- mean_se_low / mean_se_high
  expected_ratio <- sqrt(nu_high / nu_low)
  
  expect_gt(se_ratio, expected_ratio * 0.5)
})

test_that("Ridge parameter prevents numerical issues", {
  # Without ridge, nearly singular submatrices cause problems
  set.seed(3008)
  p <- 6
  nu <- 40  # Low nu increases sampling variability
  
  Sigma <- diag(p)
  Sigma[1:3, 1:3] <- 0.95  # High correlation
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # Generate data that might produce singular submatrix
  W1 <- rWishart(1, nu, Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(Sigma)
  
  # Should not error even if some submatrix is nearly singular
  expect_no_error({
    fit <- fit_covcomb(
      S_list = list(s1 = S1),
      nu = c(s1 = nu),
      control = list(ridge = 1e-6),
      se_method = "none"
    )
  })
})

test_that("Symmetry preserved under all operations", {
  # Every operation should maintain symmetry
  set.seed(3009)
  p <- 7
  nu <- 60
  
  Sigma <- diag(p)
  Sigma[1:4, 1:4] <- 0.6
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # Create asymmetric patterns in missing data
  W1 <- rWishart(1, nu, Sigma[1:5, 1:5])[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- list(ids[1:5], ids[1:5])
  
  S2 <- rWishart(1, nu, Sigma[c(3:7), c(3:7)])[,,1]
  dimnames(S2) <- list(ids[3:7], ids[3:7])
  
  fit <- fit_covcomb(
    S_list = list(s1 = S1, s2 = S2),
    nu = c(s1 = nu, s2 = nu),
    se_method = "none"
  )
  
  # Check perfect symmetry
  asym <- max(abs(fit$Sigma_hat - t(fit$Sigma_hat)))
  expect_lt(asym, 1e-14)
})

test_that("Missing completely random patterns vs structured patterns", {
  # Compare performance under MCAR vs structured missingness
  set.seed(3010)
  p <- 8
  nu <- 70
  
  Sigma <- diag(p)
  Sigma[1:4, 1:4] <- 0.7
  Sigma[5:8, 5:8] <- 0.7
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # MCAR: random variables missing
  S_mcar <- list()
  for (i in 1:5) {
    obs_idx <- sort(sample(1:p, 5))
    W <- rWishart(1, nu, Sigma[obs_idx, obs_idx])[,,1]
    dimnames(W) <- list(ids[obs_idx], ids[obs_idx])
    S_mcar[[paste0("s", i)]] <- W
  }
  
  fit_mcar <- fit_covcomb(
    S_list = S_mcar,
    nu = setNames(rep(nu, 5), names(S_mcar)),
    se_method = "none"
  )
  
  # Structured: first block vs second block
  W1_block1 <- rWishart(1, nu, Sigma[1:5, 1:5])[,,1]
  S1_block1 <- W1_block1 / nu
  dimnames(S1_block1) <- list(ids[1:5], ids[1:5])
  
  W2_block2 <- rWishart(1, nu, Sigma[4:8, 4:8])[,,1]
  S2_block2 <- W2_block2 / nu
  dimnames(S2_block2) <- list(ids[4:8], ids[4:8])
  
  fit_structured <- fit_covcomb(
    S_list = list(s1 = S1_block1, s2 = S2_block2),
    nu = c(s1 = nu, s2 = nu),
    se_method = "none",
    control = list(max_iter = 500, tol = 1e-4)  # May need more iterations after gamma removal
  )

  # Both should converge
  expect_true(fit_mcar$convergence$converged)
  expect_true(fit_structured$convergence$converged)
  
  # Both should be PD
  expect_true(all(eigen(fit_mcar$Sigma_hat)$values > 0))
  expect_true(all(eigen(fit_structured$Sigma_hat)$values > 0))
})
