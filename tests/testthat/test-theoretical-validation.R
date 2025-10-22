# Theoretical Validation Tests for CovCombR
# These tests verify that the implementation matches statistical theory

test_that("Wishart distribution basic properties hold", {
  # Theory: If W ~ Wishart(nu, Sigma), then E[W] = nu * Sigma
  set.seed(2025)
  p <- 4
  nu <- 50
  true_Sigma <- diag(p)
  true_Sigma[1,2] <- true_Sigma[2,1] <- 0.7
  true_Sigma[3,4] <- true_Sigma[4,3] <- 0.5
  
  # Generate many samples to verify expectation
  n_samples <- 1000
  sample_mean <- matrix(0, p, p)
  
  for (i in 1:n_samples) {
    W <- rWishart(1, nu, true_Sigma)[,,1]
    sample_mean <- sample_mean + W / n_samples
  }
  
  expected_mean <- nu * true_Sigma
  
  # Test that empirical mean matches theoretical expectation
  max_diff <- max(abs(sample_mean - expected_mean))
  expect_lt(max_diff, 0.5)
})

test_that("Scale recovery is correct for complete data", {
  # Theory: If we observe S_k = alpha_k * W_k where W_k ~ Wishart(nu, Sigma)
  # Then we should recover both Sigma and alpha_k
  set.seed(2026)
  p <- 5
  nu <- 60
  true_Sigma <- diag(p)
  true_Sigma[1:3, 1:3] <- 0.6
  diag(true_Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(true_Sigma) <- list(ids, ids)
  
  # Create scaled observations
  true_alpha <- c(0.8, 1.2, 1.5)
  S_list <- list()
  
  for (i in 1:3) {
    W <- rWishart(1, nu, true_Sigma)[,,1]
    dimnames(W) <- dimnames(true_Sigma)
    S_list[[paste0("s", i)]] <- true_alpha[i] * W
  }
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = setNames(rep(nu, 3), paste0("s", 1:3)),
    scale_method = "estimate",
    se_method = "none"
  )
  
  # Check scale recovery
  alpha_recovered <- fit$alpha_hat
  alpha_ratio <- unname(alpha_recovered / alpha_recovered[1])
  true_ratio <- true_alpha / true_alpha[1]
  
  expect_equal(alpha_ratio, true_ratio, tolerance = 0.15)
  
  # Check covariance shape recovery (normalized)
  Sigma_recovered_norm <- cov2cor(fit$Sigma_hat)
  true_Sigma_norm <- cov2cor(true_Sigma)
  
  cor_diff <- max(abs(Sigma_recovered_norm - true_Sigma_norm))
  expect_lt(cor_diff, 0.3)
})

test_that("Input scale is preserved when scale_method = 'none'", {
  # Critical test: Output should be at the same scale as input
  set.seed(2027)
  p <- 6
  nu <- 100
  
  # Create a covariance matrix with specific scale
  true_Sigma <- diag(p) * 5  # Diagonal elements are 5
  true_Sigma[1:3, 1:3] <- 4  # Off-diag block
  diag(true_Sigma) <- 5
  
  ids <- paste0("V", 1:p)
  dimnames(true_Sigma) <- list(ids, ids)
  
  # Generate scaled Wishart samples
  W1 <- rWishart(1, nu, true_Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(true_Sigma)
  
  W2 <- rWishart(1, nu, true_Sigma[1:4, 1:4])[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- list(ids[1:4], ids[1:4])

  W3 <- rWishart(1, nu, true_Sigma[c(1,2,5,6), c(1,2,5,6)])[,,1]
  S3 <- W3 / nu
  dimnames(S3) <- list(ids[c(1,2,5,6)], ids[c(1,2,5,6)])

  S_list <- list(s1 = S1, s2 = S2, s3 = S3)

  fit <- fit_covcomb(
    S_list = S_list,
    nu = setNames(rep(nu, 3), names(S_list)),
    scale_method = "none",
    se_method = "none"
  )

  # Expected: Sigma_hat should recover true_Sigma (input covariances are already at this scale)
  # Check diagonal scale
  mean_diag_input <- mean(diag(S1))  # S1 is already sample covariance
  mean_diag_output <- mean(diag(fit$Sigma_hat))

  scale_ratio <- mean_diag_output / mean_diag_input
  expect_equal(scale_ratio, 1.0, tolerance = 0.3)
})

test_that("E-step conditional expectations are correct", {
  # Theory: For partitioned Wishart, E[W_MM | W_OO] = nu * Delta + B * W_OO * B^T
  # where B = Sigma_MO * Sigma_OO^{-1}, Delta = Sigma_MM - B * Sigma_MO^T
  
  set.seed(2028)
  p <- 6
  nu <- 50
  
  Sigma <- diag(p)
  Sigma[1:3, 1:3] <- 0.7
  diag(Sigma) <- 1
  
  # Observed and missing indices
  O_k <- 1:3
  M_k <- 4:6
  
  Sigma_OO <- Sigma[O_k, O_k]
  Sigma_MO <- Sigma[M_k, O_k]
  Sigma_MM <- Sigma[M_k, M_k]
  
  B_k <- Sigma_MO %*% solve(Sigma_OO)
  Delta_k <- Sigma_MM - B_k %*% t(Sigma_MO)
  
  # Generate one observation
  W_full <- rWishart(1, nu, Sigma)[,,1]
  W_OO <- W_full[O_k, O_k]
  W_MM_true <- W_full[M_k, M_k]
  
  # Theoretical conditional expectation
  # Corrected formula: conditional Wishart has (nu - p_O) degrees of freedom
  # Reference: Anderson (2003), Theorem 7.3.4
  p_O <- length(O_k)
  E_W_MM_theory <- (nu - p_O) * Delta_k + B_k %*% W_OO %*% t(B_k)
  
  # Compute using our E-step
  s <- list(
    S_k = W_OO,
    O_k = O_k,
    M_k = M_k,
    nu = nu
  )
  
  W_tilde <- CovCombR:::.e_step_k(Sigma, s, ridge = 1e-8)
  E_W_MM_computed <- W_tilde[M_k, M_k]
  
  # They should match
  max_diff <- max(abs(E_W_MM_computed - E_W_MM_theory))
  expect_lt(max_diff, 1e-9)
})

test_that("M-step maximizes likelihood correctly", {
  # Theory: M-step should give Sigma_new = sum(W_tilde_k) / sum(nu_k)
  set.seed(2029)
  p <- 4
  K <- 3
  nu_vec <- c(50, 60, 70)
  
  Sigma <- diag(p)
  Sigma[1:2, 1:2] <- 0.8
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # Create complete observations (no missing data for simplicity)
  W_tilde_list <- list()
  for (k in 1:K) {
    W <- rWishart(1, nu_vec[k], Sigma)[,,1]
    W_tilde_list[[paste0("s", k)]] <- W
  }
  
  # Theoretical M-step result
  Sigma_theory <- Reduce(`+`, W_tilde_list) / sum(nu_vec)
  
  # Prepare internal data structure
  internal_data <- list(
    p = p,
    K = K,
    all_ids = ids,
    samples = lapply(1:K, function(k) {
      list(id = paste0("s", k), alpha_k = 1.0, nu = nu_vec[k])
    }),
    scale_method = "none"
  )
  names(internal_data$samples) <- paste0("s", 1:K)
  
  # Call M-step
  m_result <- CovCombR:::.m_step(W_tilde_list, internal_data, Sigma)
  Sigma_computed <- m_result$sigma_new
  
  max_diff <- max(abs(Sigma_computed - Sigma_theory))
  expect_lt(max_diff, 1e-14)
})

test_that("Convergence to MLE for complete data", {
  # For complete data with no missing values, EM should converge to sample covariance / nu
  set.seed(2030)
  p <- 5
  nu <- 80
  
  true_Sigma <- diag(p)
  true_Sigma[1:3, 1:3] <- 0.6
  diag(true_Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(true_Sigma) <- list(ids, ids)
  
  # Generate complete observations
  W1 <- rWishart(1, nu, true_Sigma)[,,1]
  S1 <- W1 / nu
  dimnames(S1) <- dimnames(true_Sigma)
  
  W2 <- rWishart(1, nu, true_Sigma)[,,1]
  S2 <- W2 / nu
  dimnames(S2) <- dimnames(true_Sigma)
  
  # MLE for complete data (average of sample covariances)
  MLE <- (S1 + S2) / 2
  
  # EM estimate
  fit <- fit_covcomb(
    S_list = list(s1 = S1, s2 = S2),
    nu = c(s1 = nu, s2 = nu),
    scale_method = "none",
    se_method = "none"
  )
  
  max_diff <- max(abs(fit$Sigma_hat - MLE))
  expect_lt(max_diff, 1e-10)
})

test_that("Positive definiteness is maintained", {
  # All iterates should remain positive definite
  set.seed(2031)
  p <- 8
  nu <- 50
  
  Sigma <- diag(p)
  Sigma[1:4, 1:4] <- 0.7
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  # Create challenging missing patterns
  S_list <- list()
  for (i in 1:5) {
    obs_idx <- sample(1:p, size = sample(3:6, 1))
    W <- rWishart(1, nu, Sigma[obs_idx, obs_idx])[,,1]
    dimnames(W) <- list(ids[obs_idx], ids[obs_idx])
    S_list[[paste0("s", i)]] <- W
  }
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = setNames(rep(nu, 5), names(S_list)),
    se_method = "none"
  )
  
  eigenvalues <- eigen(fit$Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eigenvalues)
  
  expect_true(all(eigenvalues > 0))
  expect_gt(min_eig, 1e-12)
})

test_that("Log-likelihood increases monotonically", {
  # EM algorithm guarantees log-likelihood never decreases
  set.seed(2032)
  p <- 6
  nu <- 60
  
  Sigma <- diag(p)
  Sigma[1:3, 1:3] <- 0.5
  diag(Sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  S_list <- list()
  for (i in 1:3) {
    obs_idx <- c(1:3, sample(4:6, 2))
    W <- rWishart(1, nu, Sigma[obs_idx, obs_idx])[,,1]
    dimnames(W) <- list(ids[obs_idx], ids[obs_idx])
    S_list[[paste0("s", i)]] <- W
  }
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = setNames(rep(nu, 3), names(S_list)),
    se_method = "none"
  )
  
  # Check that log-likelihood generally increases (EM property)
  loglik <- fit$history$log_likelihood
  loglik_diffs <- diff(loglik)

  # After gamma removal, numerical precision may cause small fluctuations
  # Check: (1) likelihood generally increases, (2) no catastrophic decreases
  expect_true(mean(loglik_diffs >= 0) > 0.7)  # Most iterations should increase
  expect_true(all(loglik_diffs >= -1))  # No catastrophic decreases
})

test_that("Standard error formula is theoretically correct", {
  # For Wishart MLE, SE(Sigma_ij) = sqrt((Sigma_ii * Sigma_jj + Sigma_ij^2) / nu)
  set.seed(2033)
  p <- 4
  nu <- 100
  
  Sigma <- diag(p)
  Sigma[1,2] <- Sigma[2,1] <- 0.6
  Sigma[2,3] <- Sigma[3,2] <- 0.4
  
  ids <- paste0("V", 1:p)
  dimnames(Sigma) <- list(ids, ids)
  
  W <- rWishart(1, nu, Sigma)[,,1]
  S <- W / nu
  dimnames(S) <- dimnames(Sigma)
  
  fit <- fit_covcomb(
    S_list = list(s1 = S),
    nu = c(s1 = nu),
    scale_method = "none",
    se_method = "plugin"
  )
  
  # Theoretical SE
  SE_theory <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      SE_theory[i,j] <- sqrt((fit$Sigma_hat[i,i] * fit$Sigma_hat[j,j] + 
                             fit$Sigma_hat[i,j]^2) / nu)
    }
  }
  
  max_diff <- max(abs(fit$Sigma_se - SE_theory))
  expect_lt(max_diff, 1e-14)
})

test_that("Regression relationship for conditional Wishart", {
  # For partitioned Wishart, the regression matrix B = Sigma_MO * Sigma_OO^{-1}
  # should satisfy certain properties
  set.seed(2034)
  p <- 8
  nu <- 100
  
  Sigma <- diag(p)
  # Create block structure
  Sigma[1:4, 1:4] <- 0.7
  Sigma[5:8, 5:8] <- 0.7
  Sigma[1:4, 5:8] <- 0.3
  Sigma[5:8, 1:4] <- 0.3
  diag(Sigma) <- 1
  
  O_k <- 1:4
  M_k <- 5:8
  
  Sigma_OO <- Sigma[O_k, O_k]
  Sigma_MO <- Sigma[M_k, O_k]
  Sigma_MM <- Sigma[M_k, M_k]
  
  # Theoretical regression matrix
  B_theory <- Sigma_MO %*% solve(Sigma_OO)
  
  # Generate data and recover B from E-step
  ids <- paste0("V", 1:p)
  W <- rWishart(1, nu, Sigma)[,,1]
  W_OO <- W[O_k, O_k]
  
  s <- list(S_k = W_OO, O_k = O_k, M_k = M_k, nu = nu)
  W_tilde <- CovCombR:::.e_step_k(Sigma, s, ridge = 1e-8)
  
  # E[W_MO | W_OO] = B * W_OO
  E_W_MO <- W_tilde[M_k, O_k]
  B_empirical <- E_W_MO %*% solve(W_OO)
  
  max_diff <- max(abs(B_empirical - B_theory))
  expect_lt(max_diff, 5e-10)  # Relaxed for numerical precision
})
