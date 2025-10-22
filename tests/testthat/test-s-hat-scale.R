library(testthat)

test_that("S_hat mean diagonal matches observed mean diagonal", {
  set.seed(123)
  p <- 6
  var_names <- paste0("V", seq_len(p))
  true_Sigma <- diag(p)
  dimnames(true_Sigma) <- list(var_names, var_names)
  true_Sigma[1:3, 1:3] <- 0.6
  diag(true_Sigma) <- 1

  W1 <- rWishart(1, 40, true_Sigma[1:4, 1:4])[,,1]
  S1 <- W1 / 40
  dimnames(S1) <- list(var_names[1:4], var_names[1:4])
  W2 <- rWishart(1, 50, true_Sigma[3:6, 3:6])[,,1]
  S2 <- W2 / 50
  dimnames(S2) <- list(var_names[3:6], var_names[3:6])

  S_list <- list(A = S1, B = S2)
  nu <- c(A = 40, B = 50)

  fit <- fit_covcomb(S_list, nu, se_method = "none", control = list(max_iter = 200))

  # S1 and S2 are already sample covariances
  observed_mean_s1 <- mean(diag(S1))
  observed_mean_s2 <- mean(diag(S2))

  S_hat <- fit[["S_hat"]]
  Sigma_hat <- fit[["Sigma_hat"]]

  s_hat_mean_diag <- mean(diag(S_hat))
  sigma_mean_diag <- mean(diag(Sigma_hat))

  # With gamma=1 fixed, S_hat should equal Sigma_hat
  expect_equal(S_hat, Sigma_hat, tolerance = 1e-10)
  expect_equal(s_hat_mean_diag, sigma_mean_diag, tolerance = 1e-10)

  # Verify the mean diagonal is within reasonable range of observed sample covariances
  expect_true(s_hat_mean_diag > min(observed_mean_s1, observed_mean_s2) * 0.8)
  expect_true(s_hat_mean_diag < max(observed_mean_s1, observed_mean_s2) * 1.2)
})
