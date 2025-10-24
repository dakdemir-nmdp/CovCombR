test_that("SEM method computes without error for simple case", {
  # Skip if too slow - this is just a basic smoke test
  skip_on_cran()

  # Create simple overlapping data
  set.seed(2025)
  p <- 3
  true_Sigma <- diag(p) * 2
  true_Sigma[1, 2] <- true_Sigma[2, 1] <- 0.5
  true_Sigma[2, 3] <- true_Sigma[3, 2] <- 0.7

  # Generate two samples with overlap
  S1 <- rWishart(1, 50, true_Sigma[1:2, 1:2])[, , 1] / 50
  dimnames(S1) <- list(paste0("V", 1:2), paste0("V", 1:2))

  S2 <- rWishart(1, 60, true_Sigma[2:3, 2:3])[, , 1] / 60
  dimnames(S2) <- list(paste0("V", 2:3), paste0("V", 2:3))

  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 50, s2 = 60)

  # Fit model first
  fit <- fit_covcomb(S_list, nu, se_method = "none")

  # Test SEM directly
  sem_result <- CovCombR:::compute_se_sem(
    fit_result = fit,
    S_list = S_list,
    nu = nu,
    scale_method = "none",
    alpha_normalization = "geometric"
  )

  # Check that SEM returns expected components
  expect_true(is.list(sem_result))
  expect_true("Sigma_se" %in% names(sem_result))
  expect_true("I_obs" %in% names(sem_result))
  expect_true("I_com" %in% names(sem_result))
  expect_true("R" %in% names(sem_result))

  # Check dimensions
  expect_equal(dim(sem_result$Sigma_se), c(p, p))
  expect_true(all(!is.na(sem_result$Sigma_se)))
  expect_true(all(sem_result$Sigma_se >= 0))

  # Check condition number is reported
  expect_true(is.numeric(sem_result$condition_number))
  expect_true(is.finite(sem_result$condition_number))
})

test_that("log-cholesky parameterization works in EM context", {
  # Test that we can convert back and forth
  set.seed(123)
  p <- 2
  Sigma <- matrix(c(2, 0.5, 0.5, 3), 2, 2)

  theta <- CovCombR:::sigma_to_logchol(Sigma)
  Sigma_recon <- CovCombR:::logchol_to_sigma(theta)

  expect_equal(Sigma, Sigma_recon, tolerance = 1e-10)
})

test_that("Q function can be evaluated", {
  skip_on_cran()

  # Create minimal data
  set.seed(456)
  p <- 2
  Sigma_true <- matrix(c(1, 0.3, 0.3, 1), 2, 2)

  S1 <- rWishart(1, 30, Sigma_true)[, , 1] / 30
  dimnames(S1) <- list(paste0("V", 1:2), paste0("V", 1:2))

  S_list <- list(s1 = S1)
  nu <- c(s1 = 30)

  # Preprocess
  internal_data <- CovCombR:::.preprocess_data(
    S_list, nu, "none", "geometric"
  )

  theta <- CovCombR:::sigma_to_logchol(Sigma_true)

  # Evaluate Q function
  Q_val <- CovCombR:::compute_Q_function(
    theta, theta, internal_data, ridge = 1e-8
  )

  expect_true(is.numeric(Q_val))
  expect_true(is.finite(Q_val))
  expect_equal(length(Q_val), 1)
})
