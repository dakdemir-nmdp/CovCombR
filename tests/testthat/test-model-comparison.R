# Tests for model comparison functions and normalized log-likelihood

test_that("Normalized log-likelihood is computed and added to history", {
  set.seed(123)
  p <- 5
  nu <- 50
  var_names <- paste0("V", 1:p)

  # Create true Sigma
  Sigma <- diag(p) + 0.3
  diag(Sigma) <- 1.5
  dimnames(Sigma) <- list(var_names, var_names)

  # Create two overlapping samples
  W1 <- rWishart(1, nu, Sigma[1:3, 1:3])[,,1]
  W2 <- rWishart(1, nu, Sigma[c(2,3,4,5), c(2,3,4,5)])[,,1]
  dimnames(W1) <- list(var_names[1:3], var_names[1:3])
  dimnames(W2) <- list(var_names[c(2,3,4,5)], var_names[c(2,3,4,5)])

  S_list <- list(s1 = W1, s2 = W2)
  nu_vec <- c(s1 = nu, s2 = nu)

  # Fit model with history tracking
  fit <- fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none",
                     control = list(track_history = TRUE))

  # Check that history has the new column
  expect_true("log_likelihood_per_df" %in% names(fit$history))

  # Check that normalized log-likelihood is computed correctly
  total_nu <- sum(nu_vec)
  expected_ll_per_df <- fit$history$log_likelihood / total_nu
  expect_equal(fit$history$log_likelihood_per_df, expected_ll_per_df)

  # Check that normalized values are reasonable (typically -2 to -10)
  final_ll_per_df <- tail(fit$history$log_likelihood_per_df, 1)
  expect_true(final_ll_per_df < 0)  # Should be negative
  expect_true(final_ll_per_df > -20)  # Should not be too large

  cat("\nFinal log-likelihood:", tail(fit$history$log_likelihood, 1), "\n")
  cat("Final log-likelihood per df:", final_ll_per_df, "\n")
})

test_that("Normalized log-likelihood is stable across problem sizes", {
  set.seed(456)

  results <- data.frame(
    p = integer(),
    nu = integer(),
    ll_per_df = numeric()
  )

  for (p in c(5, 10, 20)) {
    nu <- p * 10
    var_names <- paste0("V", 1:p)

    Sigma <- diag(p) + 0.3
    diag(Sigma) <- 1.5
    dimnames(Sigma) <- list(var_names, var_names)

    # Create overlapping samples
    idx1 <- 1:ceiling(2*p/3)
    idx2 <- floor(p/3):p

    W1 <- rWishart(1, nu, Sigma[idx1, idx1])[,,1]
    W2 <- rWishart(1, nu, Sigma[idx2, idx2])[,,1]
    dimnames(W1) <- list(var_names[idx1], var_names[idx1])
    dimnames(W2) <- list(var_names[idx2], var_names[idx2])

    S_list <- list(s1 = W1, s2 = W2)
    nu_vec <- c(s1 = nu, s2 = nu)

    fit <- suppressWarnings(
      fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none",
                  control = list(max_iter = 50))
    )

    final_ll_per_df <- tail(fit$history$log_likelihood_per_df, 1)
    results <- rbind(results, data.frame(p = p, nu = nu, ll_per_df = final_ll_per_df))
  }

  # Check that normalized log-likelihood doesn't grow dramatically with problem size
  # (should be relatively stable)
  ll_range <- diff(range(results$ll_per_df))
  expect_true(ll_range < 15)  # Should not vary too much

  cat("\nNormalized log-likelihood across problem sizes:\n")
  print(results)
})

test_that("compute_aic works correctly", {
  set.seed(789)
  p <- 4
  nu <- 30
  var_names <- paste0("V", 1:p)

  Sigma <- diag(p) + 0.3
  diag(Sigma) <- 1.5
  dimnames(Sigma) <- list(var_names, var_names)

  W1 <- rWishart(1, nu, Sigma[1:3, 1:3])[,,1]
  W2 <- rWishart(1, nu, Sigma[2:4, 2:4])[,,1]
  dimnames(W1) <- list(var_names[1:3], var_names[1:3])
  dimnames(W2) <- list(var_names[2:4], var_names[2:4])

  S_list <- list(s1 = W1, s2 = W2)
  nu_vec <- c(s1 = nu, s2 = nu)

  fit <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  )

  aic <- compute_aic(fit)

  # Check that AIC is computed correctly
  ll <- tail(fit$history$log_likelihood, 1)
  k <- p * (p + 1) / 2
  expected_aic <- -2 * ll + 2 * k

  expect_equal(aic, expected_aic)
  expect_true(is.finite(aic))
  expect_true(aic > 0)  # AIC should be positive for negative log-likelihoods

  cat("\nAIC:", aic, "\n")
})

test_that("compute_bic works correctly", {
  set.seed(101112)
  p <- 4
  nu <- 30
  var_names <- paste0("V", 1:p)

  Sigma <- diag(p) + 0.3
  diag(Sigma) <- 1.5
  dimnames(Sigma) <- list(var_names, var_names)

  W1 <- rWishart(1, nu, Sigma[1:3, 1:3])[,,1]
  W2 <- rWishart(1, nu, Sigma[2:4, 2:4])[,,1]
  dimnames(W1) <- list(var_names[1:3], var_names[1:3])
  dimnames(W2) <- list(var_names[2:4], var_names[2:4])

  S_list <- list(s1 = W1, s2 = W2)
  nu_vec <- c(s1 = nu, s2 = nu)

  fit <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  )

  bic <- compute_bic(fit, nu_vec)

  # Check that BIC is computed correctly
  ll <- tail(fit$history$log_likelihood, 1)
  k <- p * (p + 1) / 2
  total_nu <- sum(nu_vec)
  expected_bic <- -2 * ll + log(total_nu) * k

  expect_equal(bic, expected_bic)
  expect_true(is.finite(bic))
  expect_true(bic > 0)

  cat("\nBIC:", bic, "\n")
})

test_that("compute_aic and compute_bic require covcomb objects", {
  expect_error(compute_aic(list()), "must be an object of class 'covcomb'")
  expect_error(compute_bic(list(), c(s1 = 30)), "must be an object of class 'covcomb'")
})

test_that("compute_bic requires valid nu vector", {
  set.seed(131415)
  p <- 3
  nu <- 20
  var_names <- paste0("V", 1:p)

  Sigma <- diag(p)
  dimnames(Sigma) <- list(var_names, var_names)

  W1 <- rWishart(1, nu, Sigma)[,,1]
  dimnames(W1) <- list(var_names, var_names)

  S_list <- list(s1 = W1)
  nu_vec <- c(s1 = nu)

  fit <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  )

  # Test with invalid nu
  expect_error(compute_bic(fit, c(30)), "must be a named numeric vector")
  expect_error(compute_bic(fit, "invalid"), "must be a named numeric vector")
})

test_that("compare_models works correctly", {
  set.seed(161718)
  p <- 5
  nu <- 40
  var_names <- paste0("V", 1:p)

  Sigma <- diag(p) + 0.3
  diag(Sigma) <- 1.5
  dimnames(Sigma) <- list(var_names, var_names)

  W1 <- rWishart(1, nu, Sigma[1:3, 1:3])[,,1]
  W2 <- rWishart(1, nu, Sigma[c(2,3,4,5), c(2,3,4,5)])[,,1]
  dimnames(W1) <- list(var_names[1:3], var_names[1:3])
  dimnames(W2) <- list(var_names[c(2,3,4,5)], var_names[c(2,3,4,5)])

  S_list <- list(s1 = W1, s2 = W2)
  nu_vec <- c(s1 = nu, s2 = nu)

  # Fit two models
  fit1 <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  )
  fit2 <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "estimate", se_method = "none")
  )

  comparison <- compare_models(fit1, fit2)

  # Check structure
  expect_true(is.list(comparison))
  expect_true("log_likelihood_diff" %in% names(comparison))
  expect_true("favored_model" %in% names(comparison))
  expect_true("strength" %in% names(comparison))
  expect_true("chi_square_stat" %in% names(comparison))

  # Check values
  ll1 <- tail(fit1$history$log_likelihood, 1)
  ll2 <- tail(fit2$history$log_likelihood, 1)
  expected_diff <- ll1 - ll2

  expect_equal(comparison$log_likelihood_diff, expected_diff)
  expect_equal(comparison$chi_square_stat, -2 * expected_diff)
  expect_true(comparison$favored_model %in% c("Model 1", "Model 2"))
  expect_true(comparison$strength %in% c("Very strong", "Strong", "Moderate", "Weak"))

  cat("\nModel comparison:\n")
  cat("  Favored model:", comparison$favored_model, "\n")
  cat("  Strength:", comparison$strength, "\n")
  cat("  Log-likelihood difference:", comparison$log_likelihood_diff, "\n")
  cat("  Chi-square statistic:", comparison$chi_square_stat, "\n")
})

test_that("compare_models strength categories are correct", {
  set.seed(192021)
  p <- 4
  nu <- 30
  var_names <- paste0("V", 1:p)

  Sigma <- diag(p)
  dimnames(Sigma) <- list(var_names, var_names)

  W1 <- rWishart(1, nu, Sigma[1:3, 1:3])[,,1]
  W2 <- rWishart(1, nu, Sigma[2:4, 2:4])[,,1]
  dimnames(W1) <- list(var_names[1:3], var_names[1:3])
  dimnames(W2) <- list(var_names[2:4], var_names[2:4])

  S_list <- list(s1 = W1, s2 = W2)
  nu_vec <- c(s1 = nu, s2 = nu)

  fit <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none")
  )

  # Manually set different log-likelihood values to test strength categories
  fit1 <- fit
  fit2 <- fit
  fit3 <- fit
  fit4 <- fit

  # Create fits with known log-likelihood differences
  fit1$history$log_likelihood[nrow(fit1$history)] <- -100
  fit2$history$log_likelihood[nrow(fit2$history)] <- -101  # diff = 1 (Weak)

  comparison_weak <- compare_models(fit1, fit2)
  expect_equal(comparison_weak$strength, "Weak")

  fit2$history$log_likelihood[nrow(fit2$history)] <- -103  # diff = 3 (Moderate)
  comparison_moderate <- compare_models(fit1, fit2)
  expect_equal(comparison_moderate$strength, "Moderate")

  fit2$history$log_likelihood[nrow(fit2$history)] <- -107  # diff = 7 (Strong)
  comparison_strong <- compare_models(fit1, fit2)
  expect_equal(comparison_strong$strength, "Strong")

  fit2$history$log_likelihood[nrow(fit2$history)] <- -115  # diff = 15 (Very strong)
  comparison_vstrong <- compare_models(fit1, fit2)
  expect_equal(comparison_vstrong$strength, "Very strong")
})

test_that("compare_models requires covcomb objects", {
  expect_error(compare_models(list(), list()),
               "Both fit1 and fit2 must be objects of class 'covcomb'")
})

test_that("Model comparison helpers avoid numerical underflow", {
  # Test that we can safely compare models even with large log-likelihood values
  set.seed(222324)
  p <- 20
  nu <- 200
  var_names <- paste0("V", 1:p)

  Sigma <- diag(p) + 0.3
  diag(Sigma) <- 1.5
  dimnames(Sigma) <- list(var_names, var_names)

  # Create overlapping samples
  idx1 <- 1:14
  idx2 <- 7:20

  W1 <- rWishart(1, nu, Sigma[idx1, idx1])[,,1]
  W2 <- rWishart(1, nu, Sigma[idx2, idx2])[,,1]
  dimnames(W1) <- list(var_names[idx1], var_names[idx1])
  dimnames(W2) <- list(var_names[idx2], var_names[idx2])

  S_list <- list(s1 = W1, s2 = W2)
  nu_vec <- c(s1 = nu, s2 = nu)

  fit1 <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "none", se_method = "none",
                control = list(max_iter = 50))
  )
  fit2 <- suppressWarnings(
    fit_covcomb(S_list, nu_vec, scale_method = "estimate", se_method = "none",
                control = list(max_iter = 50))
  )

  # These should have large negative log-likelihoods
  ll1 <- tail(fit1$history$log_likelihood, 1)
  ll2 <- tail(fit2$history$log_likelihood, 1)

  expect_true(ll1 < -1000)  # Should be very negative
  expect_true(ll2 < -1000)

  # But comparison should still work
  comparison <- compare_models(fit1, fit2)

  expect_true(is.finite(comparison$log_likelihood_diff))
  expect_true(is.finite(comparison$chi_square_stat))
  expect_true(comparison$favored_model %in% c("Model 1", "Model 2"))

  # AIC and BIC should also be computable
  aic1 <- compute_aic(fit1)
  aic2 <- compute_aic(fit2)
  bic1 <- compute_bic(fit1, nu_vec)
  bic2 <- compute_bic(fit2, nu_vec)

  expect_true(is.finite(aic1) && is.finite(aic2))
  expect_true(is.finite(bic1) && is.finite(bic2))

  cat("\nLarge problem (p=20):\n")
  cat("  Log-likelihood 1:", ll1, "\n")
  cat("  Log-likelihood 2:", ll2, "\n")
  cat("  Comparison still works: difference =", comparison$log_likelihood_diff, "\n")
  cat("  AIC 1:", aic1, "  AIC 2:", aic2, "\n")
  cat("  BIC 1:", bic1, "  BIC 2:", bic2, "\n")
})
