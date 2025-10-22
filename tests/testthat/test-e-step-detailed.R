# E-step detailed tests - testing .e_step_k() function

test_that(".e_step_k() complete data returns observed block only", {
  # When all data is observed, E-step should just return the observed matrix
  set.seed(200)
  p <- 4
  sigma <- diag(p) + 0.3
  diag(sigma) <- 1
  
  # Complete observation
  nu <- 50
  S_k <- rWishart(1, nu, sigma)[,,1]
  
  sample_info <- list(
    O_k = 1:p,
    M_k = integer(0),  # No missing
    S_k = S_k,
    nu = nu
  )
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Should match S_k exactly (no imputation needed)
  expect_equal(W_tilde, S_k, tolerance = 1e-10)
})

test_that(".e_step_k() returns symmetric matrix", {
  set.seed(201)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1
  
  # Partial observation
  O_k <- 1:3
  M_k <- 4:5
  nu <- 40
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Must be symmetric
  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-12)
})

test_that(".e_step_k() observed block matches input", {
  set.seed(202)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1
  
  O_k <- 1:3
  M_k <- 4:5
  nu <- 40
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Observed block should match S_k
  expect_equal(W_tilde[O_k, O_k], S_k, tolerance = 1e-10)
})

test_that(".e_step_k() conditional expectation formula for E[W_MO|W_OO]", {
  # Test: E[W_MO | W_OO] = B_k * W_OO where B_k = Sigma_MO * Sigma_OO^{-1}
  set.seed(203)
  p <- 4
  sigma <- diag(p) + 0.5
  diag(sigma) <- 1
  
  O_k <- 1:2
  M_k <- 3:4
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Compute B_k manually
  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))
  
  # E[W_MO] should equal B_k * S_k
  E_W_MO_expected <- B_k %*% S_k
  E_W_MO_actual <- W_tilde[M_k, O_k]
  
  expect_equal(E_W_MO_actual, E_W_MO_expected, tolerance = 1e-8)
  
  # And by symmetry
  expect_equal(W_tilde[O_k, M_k], t(E_W_MO_expected), tolerance = 1e-8)
})

test_that(".e_step_k() conditional expectation formula for E[W_MM|W_OO]", {
  # Test: E[W_MM | W_OO] = nu * Delta_k + B_k * W_OO * B_k^T
  # where Delta_k = Sigma_MM - B_k * Sigma_MO^T
  set.seed(204)
  p <- 4
  sigma <- diag(p) + 0.5
  diag(sigma) <- 1
  
  O_k <- 1:2
  M_k <- 3:4
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Compute manually
  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  Sigma_MM <- sigma[M_k, M_k]
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))
  Delta_k <- Sigma_MM - B_k %*% t(Sigma_MO)

  # Corrected formula: conditional Wishart has (nu - p_O) degrees of freedom
  # Reference: Anderson (2003), Theorem 7.3.4
  p_O <- length(O_k)
  E_W_MM_expected <- (nu - p_O) * Delta_k + B_k %*% S_k %*% t(B_k)
  E_W_MM_actual <- W_tilde[M_k, M_k]

  expect_equal(E_W_MM_actual, E_W_MM_expected, tolerance = 1e-6)
})

test_that(".e_step_k() B_k matrix is regression coefficient", {
  # B_k = Sigma_MO * Sigma_OO^{-1} is the matrix of regression coefficients
  # Verify: B_k * Sigma_OO = Sigma_MO (definition of matrix inverse)
  set.seed(205)
  p <- 6
  sigma <- diag(p) + 0.3
  diag(sigma) <- 1
  
  O_k <- c(1, 3, 5)
  M_k <- c(2, 4, 6)
  nu <- 60
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  # Extract components
  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))
  
  # B_k * Sigma_OO should approximate Sigma_MO
  expect_equal(B_k %*% Sigma_OO, Sigma_MO, tolerance = 1e-8)
})

test_that(".e_step_k() Delta_k is conditional variance", {
  # Delta_k = Sigma_MM - Sigma_MO * Sigma_OO^{-1} * Sigma_OM
  # This is the Schur complement, must be PD
  set.seed(206)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1
  
  O_k <- 1:3
  M_k <- 4:5
  
  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  Sigma_MM <- sigma[M_k, M_k]
  
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))
  Delta_k <- Sigma_MM - B_k %*% t(Sigma_MO)
  
  # Delta_k must be symmetric
  expect_equal(Delta_k, t(Delta_k), tolerance = 1e-12)
  
  # Delta_k must be PD (all eigenvalues positive)
  eigs <- eigen(Delta_k, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that(".e_step_k() handles single missing variable", {
  # M_k has only one element
  set.seed(207)
  p <- 3
  sigma <- matrix(c(1, 0.6, 0.4,
                    0.6, 1, 0.5,
                    0.4, 0.5, 1), 3, 3)
  
  O_k <- 1:2
  M_k <- 3
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Should be symmetric
  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-12)
  
  # Observed block should match
  expect_equal(W_tilde[O_k, O_k], S_k, tolerance = 1e-10)
  
  # Missing block should be positive
  expect_gt(W_tilde[M_k, M_k], 0)
})

test_that(".e_step_k() handles single observed variable", {
  # O_k has only one element
  set.seed(208)
  p <- 3
  sigma <- matrix(c(1, 0.6, 0.4,
                    0.6, 1, 0.5,
                    0.4, 0.5, 1), 3, 3)
  
  O_k <- 1
  M_k <- 2:3
  nu <- 50
  S_k <- matrix(rWishart(1, nu, sigma[O_k, O_k, drop=FALSE])[,,1])
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Should be symmetric
  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-12)
  
  # Observed scalar should match
  expect_equal(W_tilde[O_k, O_k], S_k[1,1], tolerance = 1e-10)
})

test_that(".e_step_k() imputation increases with correlation", {
  # Higher correlation → larger imputed values
  set.seed(209)
  p <- 3
  nu <- 100
  O_k <- 1
  M_k <- 2:3
  
  # Low correlation
  sigma_low <- diag(3)
  sigma_low[1, 2:3] <- 0.1
  sigma_low[2:3, 1] <- 0.1
  S_k_low <- matrix(2.0)  # Observed variance
sample_info_low <- list(O_k = O_k, M_k = M_k, S_k = S_k_low, nu = nu)
  W_low <- CovCombR:::.e_step_k(sigma_low, sample_info_low, ridge = 1e-10)
  
  # High correlation
  sigma_high <- diag(3)
  sigma_high[1, 2:3] <- 0.8
  sigma_high[2:3, 1] <- 0.8
  S_k_high <- matrix(2.0)  # Same observed variance
sample_info_high <- list(O_k = O_k, M_k = M_k, S_k = S_k_high, nu = nu)
  W_high <- CovCombR:::.e_step_k(sigma_high, sample_info_high, ridge = 1e-10)
  
  # Imputed cross-covariance should be larger with high correlation
  expect_gt(abs(W_high[2, 1]), abs(W_low[2, 1]))
  expect_gt(abs(W_high[3, 1]), abs(W_low[3, 1]))
})

test_that(".e_step_k() ridge regularization prevents singularity", {
  # Create nearly singular Sigma_OO
  set.seed(210)
  p <- 4
  
  # Create ill-conditioned covariance
  eig_vals <- c(1, 0.5, 0.1, 1e-15)  # One very small eigenvalue
  eig_vecs <- qr.Q(qr(matrix(rnorm(16), 4, 4)))
  sigma <- eig_vecs %*% diag(eig_vals) %*% t(eig_vecs)
  
  O_k <- 1:3
  M_k <- 4
  nu <- 50
  
  # Generate S_k from healthy submatrix
  S_k <- rWishart(1, nu, diag(3))[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  # Should not crash with ridge
  expect_no_error({
    W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-6)
  })
  
  # Result should be symmetric and finite
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-6)
  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-10)
  expect_true(all(is.finite(W_tilde)))
})

test_that(".e_step_k() dimensions are correct", {
  set.seed(211)
  p <- 6
  sigma <- diag(p) + 0.3
  diag(sigma) <- 1
  
  O_k <- c(1, 3, 5)
  M_k <- c(2, 4, 6)
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # Should be p x p
  expect_equal(dim(W_tilde), c(p, p))
  
  # Subblocks should have correct dimensions
  expect_equal(dim(W_tilde[O_k, O_k]), c(length(O_k), length(O_k)))
  expect_equal(dim(W_tilde[M_k, M_k]), c(length(M_k), length(M_k)))
  expect_equal(dim(W_tilde[M_k, O_k]), c(length(M_k), length(O_k)))
  expect_equal(dim(W_tilde[O_k, M_k]), c(length(O_k), length(M_k)))
})

test_that(".e_step_k() result is positive definite", {
  # W_tilde should be PD (it's an expected Wishart matrix)
  set.seed(212)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1
  
  O_k <- 1:3
  M_k <- 4:5
  nu <- 60
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  
  sample_info <- list(O_k = O_k, M_k = M_k, S_k = S_k, nu = nu)
  
  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)
  
  # All eigenvalues should be positive
  eigs <- eigen(W_tilde, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})
