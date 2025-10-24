test_that("log-cholesky parameterization is bijective", {
  # Test with a simple 3x3 covariance matrix
  set.seed(123)
  p <- 3
  L_true <- matrix(c(
    1.5, 0, 0,
    0.3, 2.0, 0,
    -0.5, 0.8, 1.2
  ), nrow = 3, byrow = TRUE)
  Sigma_true <- L_true %*% t(L_true)

  # Convert to log-cholesky parameters
  theta <- CovCombR:::sigma_to_logchol(Sigma_true)

  # Expected parameter count: p + p*(p-1)/2 = 3 + 3 = 6
  expect_equal(length(theta), 6)

  # Convert back to covariance matrix
  Sigma_recon <- CovCombR:::logchol_to_sigma(theta)

  # Should recover the original matrix
  expect_equal(Sigma_recon, Sigma_true, tolerance = 1e-10)
})

test_that("log-cholesky works for 1x1 matrix", {
  Sigma <- matrix(2.5, 1, 1)
  theta <- CovCombR:::sigma_to_logchol(Sigma)

  expect_equal(length(theta), 1)
  expect_equal(theta[1], log(sqrt(2.5)))

  Sigma_recon <- CovCombR:::logchol_to_sigma(theta)
  expect_equal(Sigma_recon, Sigma, tolerance = 1e-10)
})

test_that("log-cholesky works for 2x2 matrix", {
  Sigma <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
  theta <- CovCombR:::sigma_to_logchol(Sigma)

  # Should have 2 + 1 = 3 parameters
  expect_equal(length(theta), 3)

  Sigma_recon <- CovCombR:::logchol_to_sigma(theta)
  expect_equal(Sigma_recon, Sigma, tolerance = 1e-10)
})

test_that("log-cholesky diagonal elements are log-transformed", {
  # Create diagonal matrix
  Sigma <- diag(c(1, 4, 9))
  theta <- CovCombR:::sigma_to_logchol(Sigma)

  # First 3 parameters should be log(sqrt(diagonal))
  expect_equal(theta[1:3], c(log(1), log(2), log(3)), tolerance = 1e-10)

  # Remaining parameters (off-diagonals) should be zero
  expect_equal(theta[4:6], rep(0, 3), tolerance = 1e-10)
})

test_that("jacobian_logchol_to_sigma has correct dimensions", {
  set.seed(456)
  p <- 3
  # Create a proper positive definite matrix
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- A %*% t(A) + diag(p)  # Ensure positive definite

  theta <- CovCombR:::sigma_to_logchol(Sigma)
  J <- CovCombR:::jacobian_logchol_to_sigma(theta)

  # Should be p^2 x (p + p*(p-1)/2)
  expect_equal(nrow(J), p^2)
  expect_equal(ncol(J), length(theta))
})

test_that("jacobian is numerically accurate", {
  set.seed(789)
  p <- 2
  Sigma <- matrix(c(2, 0.8, 0.8, 3), 2, 2)

  theta <- CovCombR:::sigma_to_logchol(Sigma)
  J <- CovCombR:::jacobian_logchol_to_sigma(theta)

  # Verify with manual finite differences
  h <- 1e-7
  for (k in seq_along(theta)) {
    theta_plus <- theta
    theta_plus[k] <- theta_plus[k] + h
    Sigma_plus <- CovCombR:::logchol_to_sigma(theta_plus)

    theta_minus <- theta
    theta_minus[k] <- theta_minus[k] - h
    Sigma_minus <- CovCombR:::logchol_to_sigma(theta_minus)

    fd_deriv <- as.vector((Sigma_plus - Sigma_minus) / (2 * h))

    # Compare with Jacobian column
    expect_equal(J[, k], fd_deriv, tolerance = 1e-4)
  }
})
