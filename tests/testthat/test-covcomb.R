test_that("fit_covcomb basic functionality works", {
  set.seed(123)
  p <- 5
  true_sigma <- diag(p)
  true_sigma[1,2] <- true_sigma[2,1] <- 0.8
  
  ids <- paste0("V", 1:p)
  dimnames(true_sigma) <- list(ids, ids)
  
  # Full observation case
  # Generate Wishart, then convert to sample covariance
  W_full <- rWishart(1, 50, true_sigma)[,,1]
  S_full <- W_full / 50  # Convert to sample covariance
  dimnames(S_full) <- list(ids, ids)

  fit <- fit_covcomb(
    S_list = list(full = S_full),
    nu = c(full = 50),
    scale_method = "none",
    se_method = "plugin"
  )
  
  expect_s3_class(fit, "covcomb")
  expect_true(fit$convergence$converged)
  expect_equal(dim(fit$Sigma_hat), c(p, p))
  expect_true(all(abs(fit$alpha_hat - 1) < 1e-12))
  expect_false(is.null(fit$S_hat_se))
  expect_equal(dim(fit$S_hat_se), c(p, p))
})

test_that("Missing data patterns handled correctly", {
  set.seed(456)
  p <- 8
  sigma <- diag(p)
  sigma[1:4, 1:4] <- 0.6
  diag(sigma) <- 1

  ids <- paste0("V", 1:p)
  dimnames(sigma) <- list(ids, ids)

  S_list <- list()
  for (i in 1:3) {
    W <- rWishart(1, 50, sigma)[,,1]
    S <- W / 50  # Convert to sample covariance
    dimnames(S) <- dimnames(sigma)
    obs_vars <- ids[sample(1:p, 5)]
    S_list[[paste0("s", i)]] <- S[obs_vars, obs_vars]
  }
  
  fit <- fit_covcomb(
    S_list = S_list,
    nu = setNames(rep(50, 3), paste0("s", 1:3)),
    se_method = "none"
  )
  
  expect_true(fit$convergence$converged)
  expect_true(all(eigen(fit$Sigma_hat)$values > 0))
})

test_that("Scale estimation works", {
  set.seed(789)
  p <- 6
  sigma <- diag(p) * 0.8
  diag(sigma) <- 1
  
  ids <- paste0("V", 1:p)
  dimnames(sigma) <- list(ids, ids)
  
  alpha_true <- c(1.0, 1.5, 0.8)
  S_list_scaled <- list()

  for (i in 1:3) {
    W <- rWishart(1, 60, alpha_true[i] * sigma)[,,1]  # Wishart with scaled Sigma
    S <- W / 60  # Convert to sample covariance
    dimnames(S) <- dimnames(sigma)
    obs_vars <- ids[1:4]
    S_list_scaled[[paste0("sample", i)]] <- S[obs_vars, obs_vars]
  }
  
  fit <- fit_covcomb(
    S_list = S_list_scaled,
    nu = setNames(rep(60, 3), paste0("sample", 1:3)),
    scale_method = "estimate",
    se_method = "none"
  )
  
  expect_true(fit$convergence$converged)
  expect_false(is.null(fit$alpha_hat))
  expect_length(fit$alpha_hat, 3)
  expect_true(any(abs(fit$alpha_hat - 1) > 1e-3))
})

test_that("Plugin SE computation is correct", {
  sigma_test <- diag(3)
  sigma_test[1,2] <- sigma_test[2,1] <- 0.5

  # Create coverage matrix (all entries observed with 100 df)
  coverage_mat <- matrix(100, 3, 3)

  se_result <- CovCombR:::compute_se_plugin(sigma_test, coverage_mat)

  se_11_expected <- sqrt((1*1 + 1^2) / 100)
  se_12_expected <- sqrt((0.5^2 + 1*1) / 100)

  expect_equal(se_result[1,1], se_11_expected, tolerance = 1e-10)
  expect_equal(se_result[1,2], se_12_expected, tolerance = 1e-10)
})

test_that("Symmetry is preserved", {
  set.seed(999)
  p <- 8
  sigma <- diag(p)
  sigma[1:4, 1:4] <- 0.6
  diag(sigma) <- 1

  ids <- paste0("V", 1:p)
  dimnames(sigma) <- list(ids, ids)

  S_list <- list()
  for (i in 1:3) {
    W <- rWishart(1, 50, sigma)[,,1]
    S <- W / 50  # Convert to sample covariance
    dimnames(S) <- dimnames(sigma)
    obs_vars <- ids[sample(1:p, 5)]
    S_list[[paste0("s", i)]] <- S[obs_vars, obs_vars]
  }
  
  fit <- fit_covcomb(S_list, setNames(rep(50, 3), paste0("s", 1:3)))
  
  max_asym <- max(abs(fit$Sigma_hat - t(fit$Sigma_hat)))
  expect_lt(max_asym, 1e-14)
})

test_that("S3 methods work", {
  set.seed(111)
  p <- 5
  sigma <- diag(p)

  ids <- paste0("V", 1:p)
  dimnames(sigma) <- list(ids, ids)

  W <- rWishart(1, 50, sigma)[,,1]
  S <- W / 50  # Convert to sample covariance
  dimnames(S) <- list(ids, ids)
  
  fit <- fit_covcomb(list(s1 = S), c(s1 = 50), se_method = "plugin")
  
  expect_s3_class(fit, "covcomb")
  expect_output(print(fit), "Wishart EM Covariance Combination")
  expect_output(summary(fit), "Eigenvalue spectrum \\(Sigma_hat\\)")
  expect_equal(dim(coef(fit)), c(p, p))
  expect_false(is.null(fit$Sigma_se))
  expect_false(is.null(fit$S_hat_se))
})
