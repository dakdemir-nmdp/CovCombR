#' Standard Error Calculation Methods for CovCombR
#'
#' This file contains utility functions for computing standard errors of the
#' combined covariance estimate from the EM algorithm. It includes parametric
#' bootstrap, plugin standard errors, and asymptotic methods based on observed
#' information (SEM/Oakes).
#'
#' @keywords internal

#' Parametric Bootstrap Standard Errors
#'
#' Computes standard errors for the combined covariance matrix using parametric
#' bootstrap that preserves the missing data pattern and Wishart sampling mechanism.
#'
#' @param S_list List of sample covariance matrices (with row/column names)
#' @param nu Named vector of degrees of freedom for each sample
#' @param scale_method Scaling method ("none" or "estimate")
#' @param alpha_normalization Normalization for alpha ("geometric" or "arithmetic")
#' @param init_sigma Initial covariance matrix for bootstrap fits
#' @param control Control parameters for the EM algorithm
#' @param bootstrap_ctrl Bootstrap control parameters (B, seed, progress, etc.)
#' @param Sigma_hat Converged estimate of the covariance matrix
#' @param alpha_hat Converged estimates of scale factors (if scale_method = "estimate")
#'
#' @return A list with components:
#'   \item{Sigma_se}{Matrix of standard errors for Sigma_hat}
#'   \item{samples}{Array of bootstrap covariance estimates (if retain_samples = TRUE)}
#'   \item{meta}{Bootstrap metadata (B, successes, failures, etc.)}
#'
#' @keywords internal
compute_se_bootstrap <- function(S_list, nu, scale_method = "none",
                                 alpha_normalization = "geometric", init_sigma = NULL,
                                 control = list(), bootstrap_ctrl = list(),
                                 Sigma_hat, alpha_hat = NULL) {
  if (length(S_list) == 0L) {
    stop("Bootstrap requires at least one sample.", call. = FALSE)
  }

  if (is.null(bootstrap_ctrl)) bootstrap_ctrl <- list()

  B <- bootstrap_ctrl$B
  if (is.null(B)) B <- 200
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) {
    stop("bootstrap_ctrl$B must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(round(B))

  seed <- bootstrap_ctrl$seed
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("bootstrap_ctrl$seed must be a single numeric value.", call. = FALSE)
    }
  }

  progress <- isTRUE(bootstrap_ctrl$progress)
  verbose <- isTRUE(bootstrap_ctrl$verbose)
  retain_samples <- isTRUE(bootstrap_ctrl$retain_samples)

  init_sigma_boot <- if (!is.null(bootstrap_ctrl$init_sigma)) {
    bootstrap_ctrl$init_sigma
  } else {
    Sigma_hat
  }

  control_boot <- control
  if (!is.list(control_boot)) control_boot <- list()
  if (is.list(control_boot)) control_boot$bootstrap <- NULL

  p <- nrow(Sigma_hat)
  dim_names <- dimnames(Sigma_hat)
  row_names <- if (!is.null(dim_names)) dim_names[[1]] else NULL
  col_names <- if (!is.null(dim_names)) dim_names[[2]] else NULL
  sample_count <- length(S_list)
  sample_names <- names(S_list)
  if (is.null(sample_names)) {
    sample_names <- as.character(seq_len(sample_count))
    names(S_list) <- sample_names
  }
  if (is.null(names(nu))) {
    names(nu) <- sample_names
  }
  nu <- nu[sample_names]
  if (is.null(alpha_hat)) {
    alpha_hat <- setNames(rep(1, sample_count), sample_names)
  }
  if (is.null(names(alpha_hat))) {
    alpha_hat <- setNames(alpha_hat, sample_names)
  }
  alpha_hat <- alpha_hat[sample_names]
  sample_info <- lapply(sample_names, function(id) {
    S_k <- S_list[[id]]
    if (is.null(rownames(S_k)) || is.null(colnames(S_k))) {
      stop("Each covariance matrix in S_list must have row and column names.", call. = FALSE)
    }
    obs_vars <- rownames(S_k)
    if (!setequal(obs_vars, colnames(S_k))) {
      stop("Row and column names must match for each covariance matrix.", call. = FALSE)
    }
    nu_k <- nu[[id]]
    if (is.null(nu_k) || is.na(nu_k)) {
      stop(sprintf("Degrees of freedom missing for sample '%s'.", id), call. = FALSE)
    }
    list(id = id, vars = obs_vars, nu = nu_k)
  })

  successes <- 0L
  sigma_samples <- vector("list", B)
  failed_indices <- integer(0)

  if (!is.null(seed)) {
    has_old_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_old_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    set.seed(seed)
    on.exit({
      if (has_old_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  show_progress <- progress && B > 1L
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (b in seq_len(B)) {
    S_boot <- vector("list", sample_count)
    names(S_boot) <- sample_names
    for (k in seq_len(sample_count)) {
      info <- sample_info[[k]]
      vars <- info$vars
      sigma_sub <- Sigma_hat[vars, vars, drop = FALSE]
      alpha_k <- alpha_hat[[info$id]]
      if (is.null(alpha_k) || is.na(alpha_k)) alpha_k <- 1
      # Generate Wishart matrix
      W_sim <- stats::rWishart(1, info$nu, alpha_k * sigma_sub)[,,1]
      # Convert to sample covariance (to match input format)
      sim_mat <- W_sim / info$nu
      dimnames(sim_mat) <- list(vars, vars)
      S_boot[[info$id]] <- sim_mat
    }
    nu_boot <- nu[sample_names]

    fit_boot <- tryCatch(
      fit_covcomb(
        S_list = S_boot,
        nu = nu_boot,
        scale_method = scale_method,
        alpha_normalization = alpha_normalization,
        init_sigma = init_sigma_boot,
        control = control_boot,
        se_method = "none"
      ),
      error = function(e) {
        if (verbose) {
          message(sprintf("Bootstrap replicate %d failed: %s", b, e$message))
        }
        NULL
      }
    )

    if (is.null(fit_boot) || !isTRUE(fit_boot$convergence$converged)) {
      failed_indices <- c(failed_indices, b)
      if (!is.null(fit_boot) && verbose) {
        message(sprintf(
          "Bootstrap replicate %d did not converge (iterations: %d, final change: %.3e).",
          b, fit_boot$convergence$iterations, fit_boot$convergence$final_rel_change
        ))
      }
    } else {
      successes <- successes + 1L
      sigma_boot <- fit_boot$Sigma_hat
      if (!is.null(row_names) && !is.null(col_names)) {
        sigma_boot <- sigma_boot[row_names, col_names, drop = FALSE]
      }
      sigma_samples[[successes]] <- sigma_boot
    }

    if (show_progress) utils::setTxtProgressBar(pb, b)
  }

  if (successes == 0L) {
    stop("All bootstrap replicates failed; cannot compute standard errors.", call. = FALSE)
  }

  sigma_samples <- sigma_samples[seq_len(successes)]
  sample_array <- array(NA_real_, dim = c(p, p, successes))
  for (s in seq_len(successes)) {
    sample_array[,,s] <- sigma_samples[[s]]
  }
  if (retain_samples) {
    dimnames(sample_array) <- list(row_names, col_names, paste0("b", seq_len(successes)))
  }

  sigma_se <- apply(sample_array, c(1, 2), stats::sd)
  dimnames(sigma_se) <- dim_names

  meta <- list(
    B = B,
    successes = successes,
    failures = B - successes,
    failed_indices = if (length(failed_indices)) failed_indices else NULL,
    seed = seed
  )

  list(
    Sigma_se = sigma_se,
    samples = if (retain_samples) sample_array else NULL,
    meta = meta
  )
}

#' Plugin Standard Errors
#'
#' Computes "plugin" standard errors assuming complete data and using the
#' closed-form Wishart variance with coverage weights.
#'
#' @param sigma_hat Estimated covariance matrix
#' @param coverage_mat Matrix of coverage weights (sum of degrees of freedom
#'   for pairs of variables observed together)
#'
#' @return Matrix of standard errors with the same dimensions as sigma_hat
#'
#' @keywords internal
compute_se_plugin <- function(sigma_hat, coverage_mat) {
  if (is.null(coverage_mat)) {
    stop("Coverage matrix required for plugin SE computation.", call. = FALSE)
  }
  dim_names <- dimnames(sigma_hat)
  variance_mat <- sigma_hat^2 + outer(diag(sigma_hat), diag(sigma_hat))
  se_mat <- matrix(NA_real_, nrow = nrow(sigma_hat), ncol = ncol(sigma_hat))
  positive <- coverage_mat > 0
  se_mat[positive] <- sqrt(variance_mat[positive] / coverage_mat[positive])
  if (!is.null(dim_names)) dimnames(se_mat) <- dim_names
  se_mat
}


#' Log-Cholesky Parameterization
#'
#' Converts a covariance matrix to log-Cholesky parameters.
#' Uses lower-triangular Cholesky decomposition: Sigma = L L^T where
#' L is lower triangular.
#'
#' @param Sigma Positive definite covariance matrix
#'
#' @return Vector of parameters: log(diagonals) followed by free
#'   off-diagonals (below diagonal)
#'
#' @keywords internal
sigma_to_logchol <- function(Sigma) {
  p <- nrow(Sigma)
  L <- t(chol(Sigma))  # Lower triangular Cholesky
  theta <- numeric(p + p * (p - 1) / 2)
  # Log of diagonal elements
  theta[1:p] <- log(diag(L))
  # Free off-diagonal elements (lower triangle, column-major order)
  if (p > 1) {
    idx <- p + 1
    for (j in 1:(p - 1)) {  # columns
      for (i in (j + 1):p) {  # rows below diagonal
        theta[idx] <- L[i, j]
        idx <- idx + 1
      }
    }
  }
  theta
}

#' Inverse Log-Cholesky Parameterization
#'
#' Converts log-Cholesky parameters back to a covariance matrix.
#' Reconstructs lower-triangular L such that Sigma = L L^T.
#'
#' @param theta Vector of parameters: log(diagonals) followed by
#'   free off-diagonals (below diagonal)
#'
#' @return Positive definite covariance matrix
#'
#' @keywords internal
logchol_to_sigma <- function(theta) {
  # Solve p + p*(p-1)/2 = length(theta) for p
  # This gives p^2 + p - 2*length(theta) = 0
  # Using quadratic formula: p = (-1 + sqrt(1 + 8*length(theta))) / 2
  p <- (sqrt(1 + 8 * length(theta)) - 1) / 2
  p <- as.integer(round(p))

  L <- matrix(0, p, p)
  # Diagonal elements (log-transformed)
  diag(L) <- exp(theta[1:p])
  # Off-diagonal elements (below diagonal, column-major order)
  if (p > 1) {
    idx <- p + 1
    for (j in 1:(p - 1)) {  # columns
      for (i in (j + 1):p) {  # rows below diagonal
        L[i, j] <- theta[idx]
        idx <- idx + 1
      }
    }
  }
  Sigma <- L %*% t(L)
  Sigma
}

#' Jacobian of Sigma w.r.t. log-Cholesky parameters
#'
#' Computes the Jacobian matrix d(vec(Sigma))/d(theta) for the
#' log-Cholesky parameterization using numerical differentiation.
#'
#' @param theta Vector of log-Cholesky parameters
#'
#' @return Jacobian matrix with dimensions (p^2) x length(theta)
#'
#' @keywords internal
jacobian_logchol_to_sigma <- function(theta) {
  p <- (sqrt(1 + 8 * length(theta)) - 1) / 2
  p <- as.integer(round(p))
  Sigma <- logchol_to_sigma(theta)

  # Compute Jacobian using numerical differentiation
  # This is more robust than deriving analytical formulas
  J <- matrix(0, p^2, length(theta))

  h <- 1e-8
  for (k in seq_along(theta)) {
    theta_plus <- theta
    theta_plus[k] <- theta_plus[k] + h
    Sigma_plus <- logchol_to_sigma(theta_plus)
    J[, k] <- as.vector((Sigma_plus - Sigma) / h)
  }

  J
}

#' Supplemented EM (SEM) Standard Errors
#'
#' Computes asymptotic standard errors using the Supplemented EM method
#' (Meng & Rubin, 1991). This implements the observed information matrix
#' for EM without requiring conditional second moments.
#'
#' The method works by:
#' 1. Parameterizing Sigma via log-Cholesky to ensure positive definiteness
#' 2. Computing the complete-data information I_com from the Q function
#' 3. Estimating the EM rate matrix R via finite differences of the EM map
#' 4. Forming I_obs = I_com - I_com^(1/2) R I_com^(1/2)
#' 5. Computing SEs via delta method back to Sigma scale
#'
#' @param fit_result Result from fit_covcomb containing Sigma_hat, alpha_hat
#' @param S_list Original list of sample covariance matrices
#' @param nu Original degrees of freedom vector
#' @param scale_method Scaling method used in the fit ("none" or "estimate")
#' @param alpha_normalization Alpha normalization method ("geometric" or
#'   "arithmetic")
#' @param h Finite difference step size for computing R matrix (default 1e-6)
#' @param ridge Ridge parameter for numerical stability (default 1e-8)
#'
#' @return List with components:
#'   \item{Sigma_se}{Matrix of standard errors for Sigma_hat}
#'   \item{I_obs}{Observed information matrix in log-cholesky parameterization}
#'   \item{I_com}{Complete-data information matrix}
#'   \item{R}{EM rate matrix}
#'   \item{condition_number}{Condition number of I_obs (for diagnostics)}
#'   \item{min_eigenvalue}{Smallest eigenvalue of I_obs (for diagnostics)}
#'
#' @references
#' Meng, X. L., & Rubin, D. B. (1991). Using EM to obtain asymptotic
#' variance-covariance matrices: The SEM algorithm. Journal of the American
#' Statistical Association, 86(416), 899-909.
#'
#' @keywords internal
compute_se_sem <- function(fit_result, S_list, nu,
                           scale_method = "none",
                           alpha_normalization = "geometric",
                           h = 1e-6, ridge = 1e-8) {
  Sigma_hat <- fit_result$Sigma_hat
  alpha_hat <- fit_result$alpha_hat
  p <- nrow(Sigma_hat)
  dim_names <- dimnames(Sigma_hat)

  # Step 1: Convert to log-cholesky parameterization
  theta_hat <- sigma_to_logchol(Sigma_hat)
  d_theta <- length(theta_hat)

  # Step 2: Compute complete-data information I_com
  # This requires preprocessing the data in the same way as fit_covcomb
  internal_data <- .preprocess_data(S_list, nu, scale_method,
                                    alpha_normalization)

  # Update internal_data with converged alpha values
  for (id in names(S_list)) {
    internal_data$samples[[id]]$alpha_k <- alpha_hat[[id]]
  }

  I_com <- compute_complete_data_info(theta_hat, internal_data, ridge)

  # Step 3: Estimate EM rate matrix R via finite differences
  R <- compute_em_rate_matrix(theta_hat, S_list, nu, scale_method,
                              alpha_normalization, h, ridge)

  # Step 4: Compute observed information
  # I_obs = I_com^{1/2} (I - R) I_com^{1/2}
  # Check spectral radius of R first
  eig_R <- eigen(R, symmetric = FALSE, only.values = TRUE)$values
  spectral_radius_R <- max(Mod(eig_R))

  if (spectral_radius_R >= 0.99) {
    warning(
      "EM rate matrix R has spectral radius >= 0.99 (",
      sprintf("%.3f", spectral_radius_R), "), ",
      "indicating slow EM convergence. SEM SEs may be unstable."
    )
  }

  I_com_sqrt <- matrix_sqrt(I_com, ridge)
  I_minus_R <- diag(d_theta) - R

  # Check if (I - R) is positive definite
  eig_I_minus_R <- eigen(I_minus_R, symmetric = FALSE, only.values = TRUE)$values
  if (any(Re(eig_I_minus_R) <= 0)) {
    warning(
      "(I - R) has non-positive eigenvalues. ",
      "This suggests EM is not contractive at convergence. ",
      "Falling back to I_obs = I_com (complete-data information)."
    )
    I_obs <- I_com
  } else {
    I_obs <- I_com_sqrt %*% I_minus_R %*% I_com_sqrt
  }

  # Symmetrize to correct numerical errors
  I_obs <- (I_obs + t(I_obs)) / 2

  # Check for positive definiteness
  eig_I_obs <- eigen(I_obs, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_I_obs)

  if (min_eig <= 0) {
    warning(
      "Observed information matrix is not positive definite ",
      "(min eigenvalue = ", sprintf("%.2e", min_eig), "). ",
      "This may indicate weak identifiability or insufficient overlap. ",
      "SEs may be unreliable; consider using bootstrap instead."
    )
    # Add small ridge to make it invertible
    I_obs <- I_obs + diag(abs(min_eig) + ridge, d_theta)
    eig_I_obs <- eigen(I_obs, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig_I_obs)
  }

  # Step 5: Invert to get variance of theta_hat
  Var_theta <- tryCatch({
    chol2inv(chol(I_obs))
  }, error = function(e) {
    warning("Failed to invert I_obs via Cholesky; using eigendecomposition.")
    eig_decomp <- eigen(I_obs, symmetric = TRUE)
    eig_vals_inv <- 1 / pmax(eig_decomp$values, ridge)
    eig_decomp$vectors %*% diag(eig_vals_inv) %*% t(eig_decomp$vectors)
  })

  # Step 6: Delta method to get variance of Sigma_hat
  J <- jacobian_logchol_to_sigma(theta_hat)
  Var_Sigma_vec <- J %*% Var_theta %*% t(J)

  # Extract variances for each element of Sigma
  Sigma_se <- matrix(NA_real_, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      idx <- (j - 1) * p + i  # Column-major vectorization
      Sigma_se[i, j] <- sqrt(max(Var_Sigma_vec[idx, idx], 0))
    }
  }

  dimnames(Sigma_se) <- dim_names

  # Compute condition number
  cond_num <- max(eig_I_obs) / min_eig

  list(
    Sigma_se = Sigma_se,
    I_obs = I_obs,
    I_com = I_com,
    R = R,
    condition_number = cond_num,
    min_eigenvalue = min_eig
  )
}

#' Compute Complete-Data Information Matrix
#'
#' Computes the negative Hessian of the complete-data log-likelihood
#' (Q function) with respect to log-cholesky parameters at convergence.
#' Uses numerical differentiation for robustness.
#'
#' @param theta_hat Log-cholesky parameters at convergence
#' @param internal_data Preprocessed data structure from covcomb
#' @param ridge Ridge parameter for numerical stability
#'
#' @return Complete-data information matrix (negative Hessian of Q)
#'
#' @keywords internal
compute_complete_data_info <- function(theta_hat, internal_data, ridge = 1e-8) {
  d_theta <- length(theta_hat)

  # Compute Q function at theta_hat
  Q_hat <- compute_Q_function(theta_hat, theta_hat, internal_data, ridge)

  # Numerical Hessian via central differences
  H <- matrix(0, d_theta, d_theta)
  h <- 1e-6

  for (i in 1:d_theta) {
    for (j in i:d_theta) {
      # Four-point formula for mixed partials
      theta_pp <- theta_hat
      theta_pp[i] <- theta_pp[i] + h
      theta_pp[j] <- theta_pp[j] + h
      Q_pp <- compute_Q_function(theta_pp, theta_hat, internal_data, ridge)

      theta_pm <- theta_hat
      theta_pm[i] <- theta_pm[i] + h
      theta_pm[j] <- theta_pm[j] - h
      Q_pm <- compute_Q_function(theta_pm, theta_hat, internal_data, ridge)

      theta_mp <- theta_hat
      theta_mp[i] <- theta_mp[i] - h
      theta_mp[j] <- theta_mp[j] + h
      Q_mp <- compute_Q_function(theta_mp, theta_hat, internal_data, ridge)

      theta_mm <- theta_hat
      theta_mm[i] <- theta_mm[i] - h
      theta_mm[j] <- theta_mm[j] - h
      Q_mm <- compute_Q_function(theta_mm, theta_hat, internal_data, ridge)

      H[i, j] <- (Q_pp - Q_pm - Q_mp + Q_mm) / (4 * h * h)
      if (i != j) H[j, i] <- H[i, j]
    }
  }

  # Return negative Hessian (information)
  I_com <- -H
  I_com <- (I_com + t(I_com)) / 2  # Symmetrize

  # Ensure positive definiteness
  eig_I_com <- eigen(I_com, symmetric = TRUE)
  if (min(eig_I_com$values) < ridge) {
    warning(
      "Complete-data information I_com has small/negative eigenvalues. ",
      "Adding ridge for numerical stability."
    )
    I_com <- I_com + diag(ridge, nrow(I_com))
  }

  I_com
}

#' Compute Q Function
#'
#' Computes the Q function Q(theta | theta_ref) = E[complete log-lik | obs, theta_ref]
#' This is the expected complete-data log-likelihood.
#'
#' @param theta Parameters at which to evaluate Q
#' @param theta_ref Reference parameters for taking expectation
#' @param internal_data Preprocessed data structure
#' @param ridge Ridge for numerical stability
#'
#' @return Scalar value of Q function
#'
#' @keywords internal
compute_Q_function <- function(theta, theta_ref, internal_data, ridge = 1e-8) {
  # Convert theta to Sigma
  Sigma <- logchol_to_sigma(theta)
  Sigma_ref <- logchol_to_sigma(theta_ref)
  p <- nrow(Sigma)

  # Compute E-step expectations using Sigma_ref
  Q_val <- 0

  for (s in internal_data$samples) {
    O_k <- s$O_k
    M_k <- s$M_k
    W_k <- s$W_k
    nu_k <- s$nu
    alpha_k <- s$alpha_k
    pk <- length(O_k)

    # Compute conditional expectation of full W given observed W_OO
    # using Sigma_ref for the E-step
    W_tilde <- .e_step_k(Sigma_ref, s, ridge)

    # Evaluate complete-data log-likelihood at theta (Sigma)
    # using E[W | obs, Sigma_ref] = W_tilde
    Sigma_inv <- .safe_chol_inverse(Sigma, ridge = ridge)

    log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus
    tr_term <- sum(diag(Sigma_inv %*% W_tilde))

    # Complete-data log-likelihood (ignoring constants not depending on Sigma)
    # l_c = (nu - p - 1)/2 * log|W| - nu*p/2 * log(alpha) - nu/2 * log|Sigma|
    #       - 1/(2*alpha) * tr(Sigma^{-1} W)
    # Taking expectation over W given obs:
    # Q = (nu - p - 1)/2 * E[log|W|] - nu*p/2 * log(alpha) - nu/2 * log|Sigma|
    #     - 1/(2*alpha) * tr(Sigma^{-1} E[W])

    log_det_W_tilde <- determinant(W_tilde, logarithm = TRUE)$modulus

    Q_k <- 0.5 * (nu_k - p - 1) * as.numeric(log_det_W_tilde) -
           0.5 * nu_k * p * log(alpha_k) -
           0.5 * nu_k * as.numeric(log_det_Sigma) -
           0.5 * tr_term / alpha_k

    Q_val <- Q_val + Q_k
  }

  as.numeric(Q_val)
}

#' Compute EM Rate Matrix
#'
#' Estimates the EM rate matrix R = dM/dtheta^T where M is the EM map.
#' Uses finite differences: perturb each parameter, run one EM step, and
#' compute the Jacobian.
#'
#' @param theta_hat Log-cholesky parameters at convergence
#' @param S_list Original sample covariance matrices
#' @param nu Degrees of freedom vector
#' @param scale_method Scaling method
#' @param alpha_normalization Alpha normalization method
#' @param h Step size for finite differences
#' @param ridge Ridge parameter
#'
#' @return EM rate matrix R (d_theta x d_theta)
#'
#' @keywords internal
compute_em_rate_matrix <- function(theta_hat, S_list, nu, scale_method,
                                   alpha_normalization, h = 1e-6,
                                   ridge = 1e-8) {
  d_theta <- length(theta_hat)
  R <- matrix(0, d_theta, d_theta)

  # Compute M(theta_hat) - should equal theta_hat at convergence
  M_hat <- em_map_one_step(theta_hat, S_list, nu, scale_method,
                           alpha_normalization, ridge)

  # Compute columns of R by finite differences
  for (j in 1:d_theta) {
    theta_pert <- theta_hat
    theta_pert[j] <- theta_pert[j] + h

    M_pert <- em_map_one_step(theta_pert, S_list, nu, scale_method,
                              alpha_normalization, ridge)

    R[, j] <- (M_pert - M_hat) / h
  }

  R
}

#' EM Map: One Step
#'
#' Performs one iteration of the EM algorithm starting from given parameters.
#' Returns the updated parameters after one E-step and M-step.
#'
#' @param theta Current log-cholesky parameters
#' @param S_list Sample covariance matrices
#' @param nu Degrees of freedom
#' @param scale_method Scaling method
#' @param alpha_normalization Alpha normalization
#' @param ridge Ridge parameter
#'
#' @return Updated parameters theta_new after one EM iteration
#'
#' @keywords internal
em_map_one_step <- function(theta, S_list, nu, scale_method,
                            alpha_normalization, ridge = 1e-8) {
  # Convert theta to Sigma
  Sigma <- logchol_to_sigma(theta)

  # Preprocess data
  internal_data <- .preprocess_data(S_list, nu, scale_method,
                                    alpha_normalization)

  # E-step
  w_tilde_list <- lapply(internal_data$samples, function(s) {
    .e_step_k(Sigma, s, ridge)
  })

  # M-step
  m_result <- .m_step(w_tilde_list, internal_data, Sigma, ridge = ridge)
  Sigma_new <- m_result$sigma_new

  # Convert back to log-cholesky
  theta_new <- sigma_to_logchol(Sigma_new)

  theta_new
}

#' Matrix Square Root
#'
#' Computes the symmetric square root of a positive semidefinite matrix.
#'
#' @param A Positive semidefinite matrix
#' @param ridge Ridge parameter for numerical stability
#'
#' @return Square root matrix such that A = A_sqrt %*% A_sqrt
#'
#' @keywords internal
matrix_sqrt <- function(A, ridge = 1e-8) {
  eig_decomp <- eigen(A, symmetric = TRUE)
  eig_vals_sqrt <- sqrt(pmax(eig_decomp$values, ridge))
  A_sqrt <- eig_decomp$vectors %*% diag(eig_vals_sqrt) %*%
            t(eig_decomp$vectors)
  A_sqrt
}