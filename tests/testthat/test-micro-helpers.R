# Micro-level tests for internal helper functions

test_that(".symmetrize() produces exactly symmetric matrix", {
  # Random asymmetric matrix
  A <- matrix(rnorm(16), 4, 4)
  
  A_sym <- CovCombR:::.symmetrize(A)
  
  # Test exact symmetry
  expect_equal(A_sym, t(A_sym))
  
  # Test average formula: (A + A^T)/2
  expect_equal(A_sym, (A + t(A))/2)
  
  # Test individual elements
  for (i in 1:3) {
    for (j in (i+1):4) {
      expect_equal(A_sym[i,j], A_sym[j,i])
      expect_equal(A_sym[i,j], (A[i,j] + A[j,i])/2)
    }
  }

test_that(".symmetrize() preserves already symmetric matrices", {
  A <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  A_sym <- CovCombR:::.symmetrize(A)
  
  expect_equal(A_sym, A)
})

test_that(".symmetrize() handles edge cases", {
  # 1x1 matrix
  A1 <- matrix(5)
  expect_equal(CovCombR:::.symmetrize(A1), A1)
  
  # Diagonal matrix (already symmetric)
  D <- diag(1:5)
  expect_equal(CovCombR:::.symmetrize(D), D)
  
  # Zero matrix
  Z <- matrix(0, 3, 3)
  expect_equal(CovCombR:::.symmetrize(Z), Z)
})

test_that(".project_to_pd() makes negative eigenvalues positive", {
  # Create matrix with one negative eigenvalue
  eig_vals <- c(3, 2, -0.5)  # One negative
  eig_vecs <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  A <- eig_vecs %*% diag(eig_vals) %*% t(eig_vecs)
  
  A_pd <- CovCombR:::.project_to_pd(A, min_eigen = 1e-10)
  
  # All eigenvalues should be positive
  eigs <- eigen(A_pd, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that(".project_to_pd() preserves already PD matrices", {
  # Create PD matrix
  A <- diag(3) + 0.3
  diag(A) <- 1
  
  A_pd <- CovCombR:::.project_to_pd(A, min_eigen = 1e-10)
  
  # Should be nearly identical (only symmetrization)
  expect_equal(A_pd, A, tolerance = 1e-10)
})

test_that(".project_to_pd() enforces minimum eigenvalue", {
  # Create matrix with very small positive eigenvalues
  eig_vals <- c(1, 1e-12, 1e-15)
  eig_vecs <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  A <- eig_vecs %*% diag(eig_vals) %*% t(eig_vecs)
  
  min_eigen <- 1e-8
  A_pd <- CovCombR:::.project_to_pd(A, min_eigen = min_eigen)
  
  eigs <- eigen(A_pd, symmetric = TRUE, only.values = TRUE)$values
  # Allow floating-point tolerance
  expect_true(all(eigs >= min_eigen - 1e-14))
})

test_that(".project_to_pd() mathematical properties", {
  set.seed(100)
  A <- matrix(rnorm(16), 4, 4)
  A <- A %*% t(A)  # Make symmetric PSD (but might have small eigenvalues)
  
  A_pd <- CovCombR:::.project_to_pd(A, min_eigen = 0.1)
  
  # Should be symmetric
  expect_equal(A_pd, t(A_pd), tolerance = 1e-10)
  
  # Should be PD
  eigs <- eigen(A_pd, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= 0.1))
  
  # Eigenspace should be preserved (eigenvectors same)
  eig_A <- eigen(A, symmetric = TRUE)
  eig_A_pd <- eigen(A_pd, symmetric = TRUE)
  
  # Eigenvectors should match (up to sign)
  for (i in 1:4) {
    # Check if v_i or -v_i matches
    dot_prod <- abs(sum(eig_A$vectors[,i] * eig_A_pd$vectors[,i]))
    expect_gt(dot_prod, 0.99)  # Should be close to 1
  }
})

test_that(".project_to_pd() handles edge cases", {
  # Empty matrix
  expect_equal(CovCombR:::.project_to_pd(matrix(0, 0, 0)), matrix(0, 0, 0))
  
  # 1x1 matrix (now should work!)
  A1 <- matrix(-0.5)
  A1_pd <- CovCombR:::.project_to_pd(A1, min_eigen = 0.1)
  expect_equal(A1_pd, matrix(0.1), tolerance = 1e-10)
  
  # Already good matrix
  A_good <- diag(5)
  expect_equal(CovCombR:::.project_to_pd(A_good, 0.01), A_good)
})

test_that(".preprocess_data() correctly identifies observed/missing indices", {
  p <- 6
  var_names <- paste0("V", 1:p)
  
  S1 <- diag(3)
  dimnames(S1) <- list(var_names[1:3], var_names[1:3])
  
  S2 <- diag(3)
  dimnames(S2) <- list(var_names[4:6], var_names[4:6])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 10, s2 = 15)
  
  data <- CovCombR:::.preprocess_data(S_list, nu, "none")
  
  # Check dimensions
  expect_equal(data$p, 6)
  expect_equal(data$K, 2)
  expect_equal(data$all_ids, var_names)
  
  # Check sample 1: observes V1-V3, missing V4-V6
  expect_equal(unname(data$samples$s1$O_k), 1:3)
  expect_equal(unname(data$samples$s1$M_k), 4:6)
  
  # Check sample 2: observes V4-V6, missing V1-V3
  expect_equal(unname(data$samples$s2$O_k), 4:6)
  expect_equal(unname(data$samples$s2$M_k), 1:3)
})

})

test_that(".preprocess_data() handles overlapping patterns", {
  var_names <- paste0("V", 1:5)
  
  S1 <- diag(3)
  dimnames(S1) <- list(var_names[1:3], var_names[1:3])
  
  S2 <- diag(3)
  dimnames(S2) <- list(var_names[3:5], var_names[3:5])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 10, s2 = 15)
  
  data <- CovCombR:::.preprocess_data(S_list, nu, "none")
  
  # Sample 1: V1-V3 observed, V4-V5 missing
  expect_equal(unname(data$samples$s1$O_k), 1:3)
  expect_equal(unname(data$samples$s1$M_k), 4:5)
  
  # Sample 2: V3-V5 observed, V1-V2 missing
  expect_equal(unname(data$samples$s2$O_k), 3:5)
  expect_equal(unname(data$samples$s2$M_k), 1:2)
  
  # V3 is observed in both
  expect_true(3 %in% data$samples$s1$O_k)
  expect_true(3 %in% data$samples$s2$O_k)
})

test_that(".initialize_sigma() identity method", {
  S1 <- diag(3)
  dimnames(S1) <- list(paste0("V", 1:3), paste0("V", 1:3))
  
  data <- CovCombR:::.preprocess_data(list(s1 = S1), c(s1 = 10), "none")
  
  sigma_init <- CovCombR:::.initialize_sigma(data, "identity")
  
  expect_equal(sigma_init, diag(3))
})

test_that(".initialize_sigma() avg_padded handles complete data", {
  set.seed(101)
  p <- 4
  var_names <- paste0("V", 1:p)
  
  # Two complete samples
  S1 <- diag(p) + 0.2
  diag(S1) <- 1
  dimnames(S1) <- list(var_names, var_names)
  
  S2 <- diag(p) + 0.3
  diag(S2) <- 1
  dimnames(S2) <- list(var_names, var_names)
  
  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 10, s2 = 15)
  
  data <- CovCombR:::.preprocess_data(S_list, nu, "none")
  sigma_init <- CovCombR:::.initialize_sigma(data, "avg_padded")
  
  # Should be average of the sample covariances
  expected <- (S1 + S2) / 2
  expect_equal(sigma_init, expected, tolerance = 1e-10)
})

test_that(".initialize_sigma() avg_padded handles missing data", {
  p <- 4
  var_names <- paste0("V", 1:p)
  
  # Sample 1: observes V1-V2
  S1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  dimnames(S1) <- list(var_names[1:2], var_names[1:2])
  
  # Sample 2: observes V3-V4
  S2 <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
  dimnames(S2) <- list(var_names[3:4], var_names[3:4])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 10, s2 = 15)
  
  data <- CovCombR:::.preprocess_data(S_list, nu, "none")
  sigma_init <- CovCombR:::.initialize_sigma(data, "avg_padded")
  
  # V1-V2 block should match S1 (sample covariance)
  expect_equal(sigma_init[1:2, 1:2], S1, tolerance = 1e-10)

  # V3-V4 block should match S2 (sample covariance)
  expect_equal(sigma_init[3:4, 3:4], S2, tolerance = 1e-10)
  
  # Cross-blocks V1-V2 vs V3-V4 should be 0 (no info)
  expect_equal(unname(sigma_init[1:2, 3:4]), matrix(0, 2, 2), tolerance = 1e-10)
  expect_equal(unname(sigma_init[3:4, 1:2]), matrix(0, 2, 2), tolerance = 1e-10)
  
  # Should be symmetric
  expect_equal(sigma_init, t(sigma_init), tolerance = 1e-10)
})

test_that(".initialize_sigma() avg_padded with overlapping patterns", {
  p <- 4
  var_names <- paste0("V", 1:p)
  
  # Sample 1: observes V1-V3
  S1 <- diag(3) + 0.2
  diag(S1) <- 1
  dimnames(S1) <- list(var_names[1:3], var_names[1:3])
  
  # Sample 2: observes V2-V4
  S2 <- diag(3) + 0.4
  diag(S2) <- 1
  dimnames(S2) <- list(var_names[2:4], var_names[2:4])
  
  S_list <- list(s1 = S1, s2 = S2)
  nu <- c(s1 = 10, s2 = 15)
  
  data <- CovCombR:::.preprocess_data(S_list, nu, "none")
  sigma_init <- CovCombR:::.initialize_sigma(data, "avg_padded")
  
  # V2-V3 overlap: should be averaged sample covariances
  # From S1: [[1, 0.2], [0.2, 1]]
  # From S2: [[1, 0.4], [0.4, 1]]
  # Average: [[1, 0.3], [0.3, 1]]
  expect_equal(unname(sigma_init[2:3, 2:3]),
               matrix(c(1, 0.3, 0.3, 1), 2, 2),
               tolerance = 1e-10)
})

test_that(".initialize_sigma() custom matrix", {
  custom <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  
  S1 <- diag(2)
  dimnames(S1) <- list(paste0("V", 1:2), paste0("V", 1:2))
  data <- CovCombR:::.preprocess_data(list(s1 = S1), c(s1 = 10), "none")
  
  sigma_init <- CovCombR:::.initialize_sigma(data, custom)
  
  expect_equal(sigma_init, custom)
})
