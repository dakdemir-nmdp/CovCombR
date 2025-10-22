test_that("geometric mean normalization produces product of 1", {
  # Test that with geometric mean normalization,
  # the product of all alpha_k equals 1

  set.seed(5001)
  p <- 4
  Sigma_true <- diag(p)

  # Create 3 samples with different coverage
  S1 <- rWishart(1, 50, Sigma_true[1:3, 1:3])[,,1]
  dimnames(S1) <- list(paste0("v", 1:3), paste0("v", 1:3))

  S2 <- rWishart(1, 60, Sigma_true[2:4, 2:4])[,,1]
  dimnames(S2) <- list(paste0("v", 2:4), paste0("v", 2:4))

  S3 <- rWishart(1, 40, Sigma_true[c(1,3,4), c(1,3,4)])[,,1]
  dimnames(S3) <- list(paste0("v", c(1,3,4)), paste0("v", c(1,3,4)))

  S_list <- list(s1 = S1, s2 = S2, s3 = S3)
  nu <- c(s1 = 50, s2 = 60, s3 = 40)

  # Fit with geometric mean normalization (default)
  fit_geom <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "geometric",
    se_method = "none"
  )

  # Check that product equals 1 (geometric mean constraint)
  product_alpha <- prod(fit_geom$alpha_hat)
  expect_equal(product_alpha, 1.0, tolerance = 1e-6)

  # Also verify using geometric mean directly
  K <- length(fit_geom$alpha_hat)
  geom_mean <- exp(mean(log(fit_geom$alpha_hat)))
  expect_equal(geom_mean, 1.0, tolerance = 1e-6)
})

test_that("arithmetic mean normalization produces mean of 1", {
  # Test that with arithmetic mean normalization,
  # the mean of all alpha_k equals 1

  set.seed(5002)
  p <- 4
  Sigma_true <- diag(p)

  # Create 3 samples with different coverage
  S1 <- rWishart(1, 50, Sigma_true[1:3, 1:3])[,,1]
  dimnames(S1) <- list(paste0("v", 1:3), paste0("v", 1:3))

  S2 <- rWishart(1, 60, Sigma_true[2:4, 2:4])[,,1]
  dimnames(S2) <- list(paste0("v", 2:4), paste0("v", 2:4))

  S3 <- rWishart(1, 40, Sigma_true[c(1,3,4), c(1,3,4)])[,,1]
  dimnames(S3) <- list(paste0("v", c(1,3,4)), paste0("v", c(1,3,4)))

  S_list <- list(s1 = S1, s2 = S2, s3 = S3)
  nu <- c(s1 = 50, s2 = 60, s3 = 40)

  # Fit with arithmetic mean normalization
  fit_arith <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "arithmetic",
    se_method = "none"
  )

  # Check that mean equals 1 (arithmetic mean constraint)
  mean_alpha <- mean(fit_arith$alpha_hat)
  expect_equal(mean_alpha, 1.0, tolerance = 1e-6)
})

test_that("geometric mean is default normalization", {
  # Test that geometric mean is used when alpha_normalization is not specified

  set.seed(5003)
  p <- 3
  Sigma_true <- diag(p)

  S1 <- rWishart(1, 30, Sigma_true[1:2, 1:2])[,,1]
  dimnames(S1) <- list(paste0("v", 1:2), paste0("v", 1:2))

  S2 <- rWishart(1, 40, Sigma_true[2:3, 2:3])[,,1]
  dimnames(S2) <- list(paste0("v", 2:3), paste0("v", 2:3))

  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 30, s2 = 40)

  # Fit without specifying alpha_normalization (should default to geometric)
  fit_default <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    se_method = "none"
  )

  # Should satisfy geometric mean constraint
  product_alpha <- prod(fit_default$alpha_hat)
  expect_equal(product_alpha, 1.0, tolerance = 1e-6)
})

test_that("both normalizations produce valid covariance estimates", {
  # Both normalization methods should produce valid results
  # that differ only by a global scale factor

  set.seed(5004)
  p <- 4
  Sigma_true <- matrix(c(
    1.0, 0.3, 0.2, 0.1,
    0.3, 1.0, 0.4, 0.2,
    0.2, 0.4, 1.0, 0.3,
    0.1, 0.2, 0.3, 1.0
  ), 4, 4)

  S1 <- rWishart(1, 50, Sigma_true[1:3, 1:3])[,,1]
  dimnames(S1) <- list(paste0("v", 1:3), paste0("v", 1:3))

  S2 <- rWishart(1, 60, Sigma_true[2:4, 2:4])[,,1]
  dimnames(S2) <- list(paste0("v", 2:4), paste0("v", 2:4))

  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 50, s2 = 60)

  fit_geom <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "geometric",
    se_method = "none"
  )

  fit_arith <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "arithmetic",
    se_method = "none"
  )

  # Both should converge
  expect_true(fit_geom$convergence$converged)
  expect_true(fit_arith$convergence$converged)

  # Both covariance estimates should be positive definite
  expect_true(all(eigen(fit_geom$Sigma_hat)$values > 0))
  expect_true(all(eigen(fit_arith$Sigma_hat)$values > 0))

  # The ratio between the two Sigma estimates should be constant
  # (because they differ only by a global scale factor)
  ratio_matrix <- fit_geom$Sigma_hat / fit_arith$Sigma_hat
  ratio_values <- ratio_matrix[!is.na(ratio_matrix)]
  expect_true(sd(ratio_values) / mean(ratio_values) < 0.01)
})

test_that("normalization with scale_method='none' keeps alphas at 1", {
  # When scale_method = "none", alphas should remain 1
  # regardless of normalization setting

  set.seed(5005)
  p <- 3
  Sigma_true <- diag(p)

  S1 <- rWishart(1, 30, Sigma_true[1:2, 1:2])[,,1]
  dimnames(S1) <- list(paste0("v", 1:2), paste0("v", 1:2))

  S2 <- rWishart(1, 40, Sigma_true[2:3, 2:3])[,,1]
  dimnames(S2) <- list(paste0("v", 2:3), paste0("v", 2:3))

  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 30, s2 = 40)

  fit_none <- fit_covcomb(
    S_list, nu,
    scale_method = "none",
    alpha_normalization = "geometric",
    se_method = "none"
  )

  # All alphas should be 1
  expect_true(all(abs(fit_none$alpha_hat - 1) < 1e-10))
})

test_that("invalid alpha_normalization raises error", {
  set.seed(5006)
  S1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  dimnames(S1) <- list(c("v1", "v2"), c("v1", "v2"))

  S_list <- list(s1 = S1)
  nu <- c(s1 = 10)

  expect_error(
    fit_covcomb(
      S_list, nu,
      scale_method = "estimate",
      alpha_normalization = "invalid_method"
    ),
    "should be one of"
  )
})

test_that("geometric mean handles edge case with 2 samples", {
  # With 2 samples, geometric mean should work correctly

  set.seed(5007)
  p <- 3
  Sigma_true <- diag(p)

  S1 <- rWishart(1, 50, Sigma_true[1:2, 1:2])[,,1]
  dimnames(S1) <- list(paste0("v", 1:2), paste0("v", 1:2))

  S2 <- rWishart(1, 50, Sigma_true[2:3, 2:3])[,,1]
  dimnames(S2) <- list(paste0("v", 2:3), paste0("v", 2:3))

  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 50, s2 = 50)

  fit <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "geometric",
    se_method = "none"
  )

  # Product should equal 1
  expect_equal(prod(fit$alpha_hat), 1.0, tolerance = 1e-6)

  # With symmetric setup, alphas might be similar
  # but constraint should still hold
  expect_equal(as.numeric(fit$alpha_hat[1] * fit$alpha_hat[2]), 1.0, tolerance = 1e-6)
})

test_that("geometric mean handles many samples", {
  # Test with K=5 samples

  set.seed(5008)
  p <- 6
  Sigma_true <- diag(p)

  # Create overlapping samples
  S_list <- list()
  nu_vec <- c()

  for (k in 1:5) {
    # Each sample observes 3 consecutive variables (with wraparound)
    idx <- ((k-1):(k+1)) %% p + 1
    idx <- sort(unique(idx))

    S_k <- rWishart(1, 30 + k*5, Sigma_true[idx, idx])[,,1]
    dimnames(S_k) <- list(paste0("v", idx), paste0("v", idx))

    S_list[[paste0("s", k)]] <- S_k
    nu_vec[paste0("s", k)] <- 30 + k*5
  }

  fit <- fit_covcomb(
    S_list, nu_vec,
    scale_method = "estimate",
    alpha_normalization = "geometric",
    se_method = "none"
  )

  # Product of all 5 alphas should equal 1
  expect_equal(prod(fit$alpha_hat), 1.0, tolerance = 1e-6)

  # Geometric mean should equal 1
  K <- length(fit$alpha_hat)
  geom_mean <- exp(mean(log(fit$alpha_hat)))
  expect_equal(geom_mean, 1.0, tolerance = 1e-6)
})

test_that("normalization formulas match specification", {
  # Verify the exact formulas from the specification

  set.seed(5009)
  p <- 3
  Sigma_true <- diag(p)

  S1 <- rWishart(1, 40, Sigma_true[1:2, 1:2])[,,1]
  dimnames(S1) <- list(paste0("v", 1:2), paste0("v", 1:2))

  S2 <- rWishart(1, 50, Sigma_true[2:3, 2:3])[,,1]
  dimnames(S2) <- list(paste0("v", 2:3), paste0("v", 2:3))

  S3 <- rWishart(1, 60, Sigma_true[c(1,3), c(1,3)])[,,1]
  dimnames(S3) <- list(paste0("v", c(1,3)), paste0("v", c(1,3)))

  S_list <- list(s1 = S1, s2 = S2, s3 = S3)
  nu <- c(s1 = 40, s2 = 50, s3 = 60)

  # Test geometric mean formula
  fit_geom <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "geometric",
    se_method = "none"
  )

  K <- length(fit_geom$alpha_hat)
  # Geometric mean: (∏ α_k)^(1/K) = 1
  geom_mean_actual <- prod(fit_geom$alpha_hat)^(1/K)
  expect_equal(geom_mean_actual, 1.0, tolerance = 1e-6)

  # Test arithmetic mean formula
  fit_arith <- fit_covcomb(
    S_list, nu,
    scale_method = "estimate",
    alpha_normalization = "arithmetic",
    se_method = "none"
  )

  # Arithmetic mean: (1/K) ∑ α_k = 1
  arith_mean_actual <- sum(fit_arith$alpha_hat) / K
  expect_equal(arith_mean_actual, 1.0, tolerance = 1e-6)
})
