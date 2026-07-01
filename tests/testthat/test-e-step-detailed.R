# E-step detailed tests - testing .e_step_k() function
# NOTE: .e_step_k() uses W_k = nu * S_k internally (Wishart scale, not sample
# covariance scale). All sample_info structs must supply W_k, not S_k.

test_that(".e_step_k() complete data returns observed block only", {
  # When all data is observed, E-step should return W_k (the Wishart matrix)
  set.seed(200)
  p <- 4
  sigma <- diag(p) + 0.3
  diag(sigma) <- 1

  nu <- 50
  S_k <- rWishart(1, nu, sigma)[,,1]
  W_k <- nu * S_k  # Wishart matrix

  sample_info <- list(
    O_k = 1:p,
    M_k = integer(0),  # No missing
    W_k = W_k,
    nu = nu
  )

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  # Should match W_k exactly (no imputation needed)
  expect_equal(W_tilde, W_k, tolerance = 1e-10)
})

test_that(".e_step_k() returns symmetric matrix", {
  set.seed(201)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1

  O_k <- 1:3
  M_k <- 4:5
  nu <- 40
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-12)
})

test_that(".e_step_k() observed block matches input W_k", {
  set.seed(202)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1

  O_k <- 1:3
  M_k <- 4:5
  nu <- 40
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  # Observed block should match W_k (Wishart matrix, not sample covariance)
  expect_equal(W_tilde[O_k, O_k], W_k, tolerance = 1e-10)
})

test_that(".e_step_k() conditional expectation formula for E[W_MO|W_OO]", {
  # E[W_MO | W_OO] = B_k * W_OO where B_k = Sigma_MO * Sigma_OO^{-1}
  set.seed(203)
  p <- 4
  sigma <- diag(p) + 0.5
  diag(sigma) <- 1

  O_k <- 1:2
  M_k <- 3:4
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))

  # E[W_MO] = B_k * W_OO (Wishart scale)
  E_W_MO_expected <- B_k %*% W_k
  E_W_MO_actual <- W_tilde[M_k, O_k]

  expect_equal(E_W_MO_actual, E_W_MO_expected, tolerance = 1e-8)
  expect_equal(W_tilde[O_k, M_k], t(E_W_MO_expected), tolerance = 1e-8)
})

test_that(".e_step_k() conditional expectation formula for E[W_MM|W_OO]", {
  # E[W_MM | W_OO] = (nu - p_O) * Delta_k + B_k * W_OO * B_k^T
  set.seed(204)
  p <- 4
  sigma <- diag(p) + 0.5
  diag(sigma) <- 1

  O_k <- 1:2
  M_k <- 3:4
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  Sigma_MM <- sigma[M_k, M_k]
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))
  Delta_k <- Sigma_MM - B_k %*% t(Sigma_MO)

  # Correct formula: (nu - p_O) * Delta_k + B_k * W_OO * B_k^T
  p_O <- length(O_k)
  E_W_MM_expected <- (nu - p_O) * Delta_k + B_k %*% W_k %*% t(B_k)
  E_W_MM_actual <- W_tilde[M_k, M_k]

  expect_equal(E_W_MM_actual, E_W_MM_expected, tolerance = 1e-6)
})

test_that(".e_step_k() B_k matrix is regression coefficient", {
  # B_k = Sigma_MO * Sigma_OO^{-1}
  set.seed(205)
  p <- 6
  sigma <- diag(p) + 0.3
  diag(sigma) <- 1

  O_k <- c(1, 3, 5)
  M_k <- c(2, 4, 6)
  nu <- 60
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  Sigma_OO <- sigma[O_k, O_k]
  Sigma_MO <- sigma[M_k, O_k]
  B_k <- Sigma_MO %*% solve(Sigma_OO + diag(1e-10, length(O_k)))

  expect_equal(B_k %*% Sigma_OO, Sigma_MO, tolerance = 1e-8)
})

test_that(".e_step_k() Delta_k is conditional variance", {
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

  expect_equal(Delta_k, t(Delta_k), tolerance = 1e-12)

  eigs <- eigen(Delta_k, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that(".e_step_k() handles single missing variable", {
  set.seed(207)
  p <- 3
  sigma <- matrix(c(1, 0.6, 0.4,
                    0.6, 1, 0.5,
                    0.4, 0.5, 1), 3, 3)

  O_k <- 1:2
  M_k <- 3
  nu <- 50
  S_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-12)
  expect_equal(W_tilde[O_k, O_k], W_k, tolerance = 1e-10)
  expect_gt(W_tilde[M_k, M_k], 0)
})

test_that(".e_step_k() handles single observed variable", {
  set.seed(208)
  p <- 3
  sigma <- matrix(c(1, 0.6, 0.4,
                    0.6, 1, 0.5,
                    0.4, 0.5, 1), 3, 3)

  O_k <- 1
  M_k <- 2:3
  nu <- 50
  S_k <- matrix(rWishart(1, nu, sigma[O_k, O_k, drop=FALSE])[,,1])
  W_k <- nu * S_k

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  expect_equal(W_tilde, t(W_tilde), tolerance = 1e-12)
  expect_equal(W_tilde[O_k, O_k], W_k[1, 1], tolerance = 1e-10)
})

test_that(".e_step_k() imputation increases with correlation", {
  set.seed(209)
  p <- 3
  nu <- 100
  O_k <- 1
  M_k <- 2:3

  # Low correlation
  sigma_low <- diag(3)
  sigma_low[1, 2:3] <- 0.1
  sigma_low[2:3, 1] <- 0.1
  W_k_low <- matrix(2.0 * nu)  # W_k = nu * S_k; scalar observed Wishart
  sample_info_low <- list(O_k = O_k, M_k = M_k, W_k = W_k_low, nu = nu)
  W_low <- CovCombR:::.e_step_k(sigma_low, sample_info_low, ridge = 1e-10)

  # High correlation
  sigma_high <- diag(3)
  sigma_high[1, 2:3] <- 0.8
  sigma_high[2:3, 1] <- 0.8
  W_k_high <- matrix(2.0 * nu)
  sample_info_high <- list(O_k = O_k, M_k = M_k, W_k = W_k_high, nu = nu)
  W_high <- CovCombR:::.e_step_k(sigma_high, sample_info_high, ridge = 1e-10)

  # Imputed cross-covariance should be larger with high correlation
  expect_gt(abs(W_high[2, 1]), abs(W_low[2, 1]))
  expect_gt(abs(W_high[3, 1]), abs(W_low[3, 1]))
})

test_that(".e_step_k() ridge regularization prevents singularity", {
  set.seed(210)
  p <- 4

  eig_vals <- c(1, 0.5, 0.1, 1e-15)
  eig_vecs <- qr.Q(qr(matrix(rnorm(16), 4, 4)))
  sigma <- eig_vecs %*% diag(eig_vals) %*% t(eig_vecs)

  O_k <- 1:3
  M_k <- 4
  nu <- 50
  W_k <- rWishart(1, nu, diag(3))[,,1]

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  expect_no_error({
    W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-6)
  })

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
  W_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  expect_equal(dim(W_tilde), c(p, p))
  expect_equal(dim(W_tilde[O_k, O_k]), c(length(O_k), length(O_k)))
  expect_equal(dim(W_tilde[M_k, M_k]), c(length(M_k), length(M_k)))
  expect_equal(dim(W_tilde[M_k, O_k]), c(length(M_k), length(O_k)))
  expect_equal(dim(W_tilde[O_k, M_k]), c(length(O_k), length(M_k)))
})

test_that(".e_step_k() result is positive definite", {
  set.seed(212)
  p <- 5
  sigma <- diag(p) + 0.4
  diag(sigma) <- 1

  O_k <- 1:3
  M_k <- 4:5
  nu <- 60
  W_k <- rWishart(1, nu, sigma[O_k, O_k])[,,1]

  sample_info <- list(O_k = O_k, M_k = M_k, W_k = W_k, nu = nu)

  W_tilde <- CovCombR:::.e_step_k(sigma, sample_info, ridge = 1e-10)

  eigs <- eigen(W_tilde, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})
