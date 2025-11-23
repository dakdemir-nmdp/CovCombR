#' CovCombR: EM Algorithm for Combining Incomplete Covariance Matrices
#'
#' @description
#' This package implements the Expectation-Maximization (EM) algorithm to recover
#' a combined covariance (or genetic relationship) matrix from incomplete covariance
#' inputs. The algorithm estimates the Wishart scale parameter (per-df covariance)
#' and provides an aggregated estimate suitable for downstream prediction tasks.
#'
#' @details
#' The main function is \code{\link{fit_covcomb}}, which estimates a covariance
#' matrix from a collection of incomplete Wishart matrices. The algorithm uses:
#' \itemize{
#'   \item E-step: Conditional expectations for partitioned Wishart distributions
#'   \item M-step: Maximum likelihood updates for covariance and scale parameters
#'   \item Numerical stability: Cholesky decomposition for matrix inversions
#' }
#'
#' @references
#' McLachlan, G. J., & Krishnan, T. (2008). The EM Algorithm and Extensions (2nd ed.).
#' Wiley Series in Probability and Statistics.
#'
#' @importFrom stats cov2cor setNames
#' @keywords internal
"_PACKAGE"

# -- S3 Methods --------------------------------

#' Print method for covcomb objects
#'
#' @param x An object of class \code{covcomb}
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object
#' @export
#' @method print covcomb
print.covcomb <- function(x, ...) {
  cat("Wishart EM Covariance Combination\n")
  cat("Call: ")
  print(x$call)

  status <- if (x$convergence$converged) "converged" else "did not converge"
  cat(sprintf(
    "\nStatus: %s in %d iterations\n",
    status, x$convergence$iterations
  ))
  cat(sprintf("Final relative change: %.2e\n", x$convergence$final_rel_change))
  cat(sprintf(
    "\nCombined covariance matrix (Sigma_hat): %d x %d\n",
    nrow(x$Sigma_hat), ncol(x$Sigma_hat)
  ))
  cat(sprintf("  Mean diagonal: %.4f\n", mean(diag(x$Sigma_hat))))
  cat(sprintf("  Condition number: %.2e\n", kappa(x$Sigma_hat, exact = FALSE)))
  if (!is.null(x$Sigma_se)) {
    cat(sprintf("  Mean SE: %.4f\n", mean(x$Sigma_se, na.rm = TRUE)))
  }
  if (x$S_hat_scale != 1) {
    cat(sprintf("\nS_hat (rescaled): %d x %d\n", nrow(x$S_hat), ncol(x$S_hat)))
    cat(sprintf("  S_hat_scale: %.4f\n", x$S_hat_scale))
    cat(sprintf("  Mean diagonal: %.4f\n", mean(diag(x$S_hat))))
  }
  if (!is.null(x$alpha_hat) && !all(abs(x$alpha_hat - 1) < 1e-10)) {
    cat("\nEstimated scale factors (alpha_hat):\n")
    print(round(sort(x$alpha_hat), 4))
  }

  invisible(x)
}

#' Summary method for covcomb objects
#'
#' @param object An object of class \code{covcomb}
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object
#' @export
#' @method summary covcomb
summary.covcomb <- function(object, ...) {
  print(object)

  cat("\nEigenvalue spectrum (Sigma_hat):\n")
  eig_sigma <- eigen(object$Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  print(summary(eig_sigma))

  if (!is.null(object$Sigma_se)) {
    cat(sprintf("\nStandard errors (method: %s):\n", object$se_method))
    print(summary(as.numeric(object$Sigma_se)))
  }

  invisible(object)
}

#' Extract covariance matrix from covcomb object
#'
#' @param object An object of class \code{covcomb}
#' @param ... Additional arguments (unused)
#' @return The estimated covariance matrix
#' @export
#' @method coef covcomb
coef.covcomb <- function(object, ...) {
  object$Sigma_hat
}

#' Fitted covariance matrix for covcomb object
#'
#' @param object An object of class \code{covcomb}
#' @param ... Additional arguments (unused)
#' @return The combined covariance estimate on the aggregated scale (S_hat)
#' @export
#' @method fitted covcomb
fitted.covcomb <- function(object, ...) {
  object$S_hat
}

# -- Main Function --------------------------------

#' Fit Wishart EM Model
#'
#' @description
#' Estimates a common covariance matrix from incomplete sample covariance matrices using
#' the EM algorithm. Handles missing data patterns and heterogeneous scaling.
#'
#' @param S_list Named list of sample covariance matrices (e.g., from \code{cov(X)}).
#'   Each matrix is a principal submatrix with row/column names identifying observed variables.
#' @param nu Named numeric vector of degrees of freedom (sample sizes) for each sample.
#'   Must be at least the number of observed variables for each corresponding sample.
#' @param scale_method Scaling method: \code{"none"} (default) or \code{"estimate"}
#' @param alpha_normalization Normalization method for scale factors when
#'   \code{scale_method = "estimate"}. Options:
#'   \itemize{
#'     \item \code{"geometric"} (default): Constrains the geometric mean to 1,
#'           i.e., \eqn{(\prod_{k=1}^K \alpha_k)^{1/K} = 1}. Preferred for multiplicative effects.
#'     \item \code{"arithmetic"}: Constrains the arithmetic mean to 1,
#'           i.e., \eqn{(1/K) \sum_{k=1}^K \alpha_k = 1}.
#'   }
#'   This parameter is ignored when \code{scale_method = "none"}.
#' @param init_sigma Initial covariance matrix or initialization method
#'   (\code{"avg_padded"} or \code{"identity"})
#' @param control List of control parameters:
#'   \itemize{
#'     \item \code{max_iter}: Maximum iterations (default: 500)
#'     \item \code{tol}: Convergence tolerance (default: 1e-7)
#'     \item \code{ridge}: Ridge parameter for stability (default: 1e-8)
#'     \item \code{min_eigen}: Minimum eigenvalue threshold (default: 1e-10)
#'     \item \code{bootstrap}: Optional list controlling bootstrap SE computation with elements:
#'       \itemize{
#'         \item \code{B}: Number of bootstrap replicates (default: 200)
#'         \item \code{seed}: Optional RNG seed for reproducibility
#'         \item \code{progress}: Logical; print simple progress updates if \code{TRUE}
#'         \item \code{verbose}: Logical; emit messages for failed replicates when \code{TRUE}
#'         \item \code{retain_samples}: Logical; keep the bootstrap sample array in the result when \code{TRUE}
#'         \item \code{init_sigma}: Initialization passed to bootstrap fits (defaults to the converged \code{Sigma_hat})
#'       }
#'     \item \code{sem}: Optional list controlling SEM SE computation with elements:
#'       \itemize{
#'         \item \code{h}: Finite difference step size for computing EM rate matrix (default: 1e-6)
#'         \item \code{ridge}: Ridge parameter for numerical stability (default: uses \code{control$ridge})
#'       }
#'   }
#' @param se_method Standard error method: \code{"none"}, \code{"plugin"} (default),
#'   \code{"bootstrap"}, or \code{"sem"} (experimental).
#'
#' @details
#' The EM algorithm iterates between:
#' \itemize{
#'   \item \strong{E-step}: Computes conditional expectations of missing Wishart blocks
#'   \item \strong{M-step}: Updates covariance matrix and optional scale factors
#' }
#'
#' \strong{Standard Errors:}
#' \itemize{
#'   \item \code{"plugin"}: Fast, using the closed-form Wishart variance with per-entry
#'         coverage weights. \strong{Important:} Plugin SEs assume complete data and
#'         only account for sampling variance, not EM imputation uncertainty. They may
#'         substantially underestimate true uncertainty for entries with missing data.
#'         Returns \code{NA} for entries never jointly observed. Use for exploratory
#'         analysis only.
#'   \item \code{"bootstrap"}: Parametric bootstrap that simulates new incomplete Wishart
#'         samples from the fitted model \code{B} times and refits the EM algorithm for
#'         each replicate. Accounts for both sampling variance and imputation uncertainty.
#'         \strong{Recommended for formal inference and hypothesis testing.}
#'   \item \code{"sem"} (experimental): Supplemented EM method (Meng & Rubin, 1991) that
#'         computes asymptotic standard errors using the observed information matrix.
#'         Faster than bootstrap but may be less reliable for small samples or weak overlap.
#'         Uses log-Cholesky parameterization and finite differences to estimate the EM
#'         rate matrix. Returns diagnostic information (condition number, eigenvalues).
#'         \strong{Use with caution and validate against bootstrap for critical applications.}
#' }
#'
#' Built-in bootstrap SEs can be requested with \code{se_method = "bootstrap"}. The algorithm generates
#' \code{B} synthetic datasets (default: 200) matching the original missing-pattern structure, refits the EM model for each
#' replicate, and reports the empirical standard deviation of the resulting covariance estimates.
#' Bootstrap behaviour can be customised through \code{control$bootstrap}; see the argument details.
#'
#' SEM standard errors can be requested with \code{se_method = "sem"}. This method is experimental
#' and provides fast asymptotic approximations. Control parameters can be specified via
#' \code{control$sem} with elements \code{h} (finite difference step, default 1e-6) and
#' \code{ridge} (regularization, default uses \code{control$ridge}).
#'
#' @return Object of class \code{covcomb} containing:
#' \item{Sigma_hat}{Estimated combined covariance matrix. This is the primary output,
#'   at the same scale as the input sample covariances. Use this for downstream analysis
#'   and prediction, or with \code{Sigma_se} for inference.}
#' \item{S_hat}{Rescaled presentation of \code{Sigma_hat}. When \code{scale_method = "none"}
#'   or \code{alpha_normalization = "geometric"}, this equals \code{Sigma_hat}. When
#'   \code{scale_method = "estimate"} with \code{alpha_normalization = "arithmetic"},
#'   this equals \eqn{S\_hat\_scale \times \Sigma_{hat}}. Users should primarily use
#'   \code{Sigma_hat} for analysis; \code{S_hat} is provided for compatibility.}
#' \item{S_hat_scale}{Scalar rescaling factor applied to Sigma_hat. Equals 1 when
#'   \code{scale_method = "none"} or \code{alpha_normalization = "geometric"}.}
#' \item{alpha_hat}{Scale factors for each input sample. Equals 1 when \code{scale_method != "estimate"}.}
#' \item{Sigma_se}{Standard error matrix on the per-df scale (if \code{se_method != "none"}).
#'   This provides standard errors for \code{Sigma_hat}.}
#' \item{S_hat_se}{Standard error matrix on the data scale (if \code{se_method != "none"})}
#' \item{bootstrap}{Bootstrap metadata (if \code{se_method = "bootstrap"})}
#' \item{bootstrap_samples}{Bootstrap covariance estimates (if \code{retain_samples = TRUE})}
#' \item{sem}{SEM diagnostics (if \code{se_method = "sem"}), including:
#'   \itemize{
#'     \item \code{I_obs}: Observed information matrix in log-Cholesky parameterization
#'     \item \code{I_com}: Complete-data information matrix
#'     \item \code{R}: EM rate matrix
#'     \item \code{condition_number}: Condition number of observed information (for diagnostics)
#'     \item \code{min_eigenvalue}: Smallest eigenvalue of observed information (for diagnostics)
#'   }}
#' \item{convergence}{Convergence information (status, iterations, final change)}
#' \item{history}{Iteration history (relative change, log-likelihood with constant terms)}
#' \item{call}{Matched call}
#'
#' @examples
#' # Simulate incomplete sample covariance matrices
#' set.seed(2025)
#' p <- 8
#' var_names <- paste0("V", 1:p)
#' true_Sigma <- diag(p)
#' dimnames(true_Sigma) <- list(var_names, var_names)
#' true_Sigma[1:4, 1:4] <- 0.7
#' diag(true_Sigma) <- 1
#'
#' # Generate 3 sample covariances with different missing patterns
#' # Simulate from Wishart, then convert to sample covariances
#' W1 <- stats::rWishart(1, 50, true_Sigma[1:5, 1:5])[, , 1]
#' S1 <- W1 / 50
#' dimnames(S1) <- list(var_names[1:5], var_names[1:5])
#'
#' W2 <- stats::rWishart(1, 60, true_Sigma[3:8, 3:8])[, , 1]
#' S2 <- W2 / 60
#' dimnames(S2) <- list(var_names[3:8], var_names[3:8])
#'
#' W3 <- stats::rWishart(1, 55, true_Sigma[c(1, 3, 5, 7), c(1, 3, 5, 7)])[, , 1]
#' S3 <- W3 / 55
#' dimnames(S3) <- list(var_names[c(1, 3, 5, 7)], var_names[c(1, 3, 5, 7)])
#'
#' S_list <- list(sample1 = S1, sample2 = S2, sample3 = S3)
#' nu <- c(sample1 = 50, sample2 = 60, sample3 = 55)
#'
#' # Fit model
#' fit <- fit_covcomb(S_list, nu, se_method = "plugin")
#' print(fit)
#'
#' # Extract combined covariance estimate (same scale as input sample covariances)
#' Sigma_combined <- fit$Sigma_hat # or coef(fit) or fit$S_hat
#'
#' # The output is at the same scale as input sample covariances
#' # Compare to true population covariance
#' print(Sigma_combined)
#' print(true_Sigma)
#'
#' # Standard errors (if requested)
#' if (!is.null(fit$Sigma_se)) {
#'   print(fit$Sigma_se)
#' }
#'
#' # Bootstrap standard errors (more accurate, ~100x slower)
#' # fit_boot <- fit_covcomb(
#' #   S_list, nu,
#' #   se_method = "bootstrap",
#' #   control = list(bootstrap = list(B = 200, seed = 2025, progress = TRUE))
#' # )
#' # fit_boot$Sigma_se  # Bootstrap SEs on per-df scale

#'
#' @export

fit_covcomb <- function(S_list, nu,
                        scale_method = "none",
                        alpha_normalization = "geometric",
                        init_sigma = "identity",
                        control = list(),
                        se_method = "plugin") {
  call <- match.call()

  scale_method <- match.arg(scale_method, c("none", "estimate"))
  alpha_normalization <- match.arg(alpha_normalization, c("geometric", "arithmetic"))
  se_method <- match.arg(se_method, c("none", "plugin", "bootstrap", "sem"))

  # Validate inputs
  stopifnot(is.list(S_list), !is.null(names(S_list)))
  stopifnot(is.numeric(nu), !is.null(names(nu)))
  stopifnot(all(names(S_list) %in% names(nu)))
  nu <- nu[names(S_list)]
  if (any(is.na(nu))) {
    stop("Degrees of freedom vector must include entries for all samples in S_list.", call. = FALSE)
  }

  # Set controls
  ctrl <- list(max_iter = 500, tol = 1e-7, ridge = 1e-8, min_eigen = 1e-10)
  ctrl[names(control)] <- control

  # Pre-process data
  # Convert sample covariances to Wishart matrices internally
  internal_data <- .preprocess_data(S_list, nu, scale_method, alpha_normalization)

  # Check graph connectivity
  observed_sets <- lapply(internal_data$samples, function(s) s$O_k)
  connectivity <- .check_graph_connectivity(internal_data$p, observed_sets)

  if (!connectivity$is_connected) {
    warning(
      sprintf(
        "The observation pattern creates %d disconnected component%s. ",
        connectivity$num_components,
        if (connectivity$num_components > 1) "s" else ""
      ),
      "Variables in different components cannot be jointly estimated. ",
      "Consider using only samples that observe overlapping variable sets.",
      call. = FALSE
    )
  }

  # Initialize
  sigma_current <- .initialize_sigma(internal_data, init_sigma)
  sigma_current <- .project_to_pd(sigma_current, ctrl$min_eigen)

  # EM loop
  history <- data.frame(
    iteration = integer(), rel_change = numeric(),
    log_likelihood = numeric()
  )
  converged <- FALSE
  rel_change <- NA_real_
  iterations_used <- 0L

  for (i in 1:ctrl$max_iter) {
    # E-step
    w_tilde_list <- lapply(internal_data$samples, function(s) {
      .e_step_k(sigma_current, s, ctrl$ridge)
    })

    # M-step
    m_result <- .m_step(w_tilde_list, internal_data, sigma_current, ridge = ctrl$ridge)
    sigma_new <- m_result$sigma_new
    internal_data <- m_result$internal_data

    # Ensure PD
    sigma_new <- .project_to_pd(sigma_new, ctrl$min_eigen)

    denom <- max(norm(sigma_current, "F"), .Machine$double.eps)
    rel_change <- norm(sigma_new - sigma_current, "F") / denom
    loglik <- .compute_loglik(sigma_new, internal_data)
    history[i, ] <- list(
      iteration = i,
      rel_change = rel_change,
      log_likelihood = loglik
    )

    sigma_current <- sigma_new
    iterations_used <- i

    if (rel_change < ctrl$tol) {
      converged <- TRUE
      break
    }
  }

  # Format output
  Sigma_hat <- sigma_current
  dimnames(Sigma_hat) <- list(internal_data$all_ids, internal_data$all_ids)

  # Compute alpha_hat
  sample_names <- names(S_list)
  if (is.null(sample_names)) sample_names <- as.character(seq_along(S_list))
  alpha_hat <- if (scale_method == "estimate") {
    setNames(sapply(internal_data$samples, `[[`, "alpha_k"), sample_names)
  } else {
    setNames(rep(1, length(S_list)), sample_names)
  }

  # Compute S_hat on data scale
  if (scale_method == "estimate") {
    # With geometric normalization, mean(alpha_hat) ≈ 1 by construction
    # so S_hat_scale adds no information. Set to 1 for clarity.
    if (alpha_normalization == "geometric") {
      S_hat_scale <- 1
    } else {
      # With arithmetic normalization, use the mean for rescaling
      S_hat_scale <- mean(alpha_hat)
    }
    S_hat <- S_hat_scale * Sigma_hat
  } else {
    S_hat_scale <- 1
    S_hat <- Sigma_hat
  }
  dimnames(S_hat) <- dimnames(Sigma_hat)

  coverage_mat <- .compute_coverage(internal_data)

  # Check for identifiability issues
  K <- length(S_list)
  min_eig <- min(eigen(Sigma_hat, symmetric = TRUE, only.values = TRUE)$values)
  cond_num <- kappa(Sigma_hat, exact = FALSE)

  if (K == 1) {
    warning("Only one sample provided. The model is not identifiable - ",
      "unobserved variables will have arbitrary values. ",
      "Consider providing multiple samples with overlapping variables.",
      call. = FALSE
    )
  }

  if (cond_num > 1e6) {
    warning("Estimated covariance matrix is very ill-conditioned (condition number = ",
      sprintf("%.2e", cond_num), "). ",
      "This often indicates insufficient overlap between samples. ",
      "Check that variable pairs are jointly observed in at least one sample.",
      call. = FALSE
    )
  } else if (cond_num > 1e4) {
    message(
      "Note: Estimated covariance has high condition number (",
      sprintf("%.2e", cond_num), "). ",
      "Results may be sensitive to missing data patterns."
    )
  }

  if (min_eig < ctrl$min_eigen * 10) {
    warning("Smallest eigenvalue (", sprintf("%.2e", min_eig), ") is very close to the ",
      "minimum threshold. This may indicate rank deficiency or poor identifiability.",
      call. = FALSE
    )
  }

  result <- list(
    Sigma_hat = Sigma_hat,
    S_hat = S_hat,
    S_hat_scale = S_hat_scale,
    alpha_hat = alpha_hat,
    se_method = se_method,
    Sigma_se = NULL,
    S_hat_se = NULL,
    bootstrap = NULL,
    bootstrap_samples = NULL,
    convergence = list(
      converged = converged,
      iterations = iterations_used,
      final_rel_change = rel_change
    ),
    history = history,
    call = call
  )

  # Compute standard errors
  if (se_method == "plugin") {
    result$Sigma_se <- compute_se_plugin(result$Sigma_hat, coverage_mat)
    # Warn about plugin SE limitations
    message(
      "Note: Plugin standard errors assume complete data and may underestimate ",
      "uncertainty for entries with substantial missing data. ",
      "They account for sampling variance but not EM imputation uncertainty. ",
      "For accurate inference, use se_method = 'bootstrap'."
    )
  } else if (se_method == "bootstrap") {
    boot_res <- compute_se_bootstrap(
      S_list = S_list,
      nu = nu,
      scale_method = scale_method,
      alpha_normalization = alpha_normalization,
      init_sigma = init_sigma,
      control = control,
      bootstrap_ctrl = ctrl$bootstrap,
      Sigma_hat = result$Sigma_hat,
      alpha_hat = result$alpha_hat
    )
    result$Sigma_se <- boot_res$Sigma_se
    result$bootstrap <- boot_res$meta
    if (!is.null(boot_res$samples)) {
      result$bootstrap_samples <- boot_res$samples
    }
  } else if (se_method == "sem") {
    message(
      "Note: SEM standard errors are experimental. ",
      "They provide fast asymptotic approximations but may be unreliable ",
      "for small samples or weak overlap. Bootstrap is recommended for formal inference."
    )

    sem_ctrl <- if (!is.null(ctrl$sem)) ctrl$sem else list()
    sem_h <- if (!is.null(sem_ctrl$h)) sem_ctrl$h else 1e-6
    sem_ridge <- if (!is.null(sem_ctrl$ridge)) sem_ctrl$ridge else ctrl$ridge

    sem_res <- compute_se_sem(
      fit_result = result,
      S_list = S_list,
      nu = nu,
      scale_method = scale_method,
      alpha_normalization = alpha_normalization,
      h = sem_h,
      ridge = sem_ridge
    )

    result$Sigma_se <- sem_res$Sigma_se
    result$sem <- list(
      I_obs = sem_res$I_obs,
      I_com = sem_res$I_com,
      R = sem_res$R,
      condition_number = sem_res$condition_number,
      min_eigenvalue = sem_res$min_eigenvalue
    )
  }

  if (!is.null(result$Sigma_se)) {
    result$S_hat_se <- result$S_hat_scale * result$Sigma_se
  }

  structure(result, class = "covcomb")
}

# -- Internal Functions --------------------------------

# -- Helper Functions for Numerical Stability --------------------------------

#' Safe Cholesky Decomposition and Inversion
#'
#' @param A Matrix to invert via Cholesky
#' @param ridge Initial ridge value for regularization
#' @param max_attempts Maximum number of ridge increase attempts
#' @param return_chol Logical; if TRUE, return a list with both inverse and Cholesky factor
#' @return If return_chol=FALSE (default), returns the inverse of A.
#'   If return_chol=TRUE, returns a list with components:
#'   \itemize{
#'     \item inverse: Inverse of A computed as chol2inv(chol(A + ridge*I))
#'     \item chol: Upper triangular Cholesky factor
#'     \item ridge_used: Final ridge value applied
#'   }
#' @keywords internal
.safe_chol_inverse <- function(A, ridge = 1e-8, max_attempts = 5, return_chol = FALSE) {
  p <- nrow(A)
  ridge_current <- ridge

  for (attempt in seq_len(max_attempts)) {
    chol_try <- tryCatch(
      chol(A + diag(ridge_current, p)),
      error = function(e) NULL
    )
    if (!is.null(chol_try)) {
      inv <- chol2inv(chol_try)
      if (return_chol) {
        return(list(inverse = inv, chol = chol_try, ridge_used = ridge_current))
      } else {
        return(inv)
      }
    }
    if (attempt < max_attempts) {
      ridge_current <- ridge_current * 10
    }
  }

  stop(sprintf(
    "Failed Cholesky decomposition after %d attempts; final ridge value tried was %.2e",
    max_attempts, ridge_current
  ))
}

# -- Data Processing Functions --------------------------------

.preprocess_data <- function(S_list, nu, scale_method, alpha_normalization = "geometric") {
  sample_ids <- names(S_list)
  all_ids <- sort(unique(unlist(lapply(S_list, rownames))))
  p <- length(all_ids)
  id_map <- setNames(1:p, all_ids)

  samples <- lapply(sample_ids, function(id) {
    S_k_samp <- S_list[[id]]
    stopifnot(all(rownames(S_k_samp) == colnames(S_k_samp)))
    nu_k <- nu[[id]]
    if (is.null(nu_k) || !is.finite(nu_k) || nu_k <= 0) {
      stop(sprintf("Degrees of freedom for sample '%s' must be positive.", id), call. = FALSE)
    }

    # Convert sample covariance to Wishart matrix: W = nu * S_samp
    # This allows the algorithm to work with E[W] = nu * Sigma
    W_k <- nu_k * S_k_samp

    list(
      id = id,
      W_k = W_k, # Store Wishart matrix for internal EM algorithm
      O_k = id_map[rownames(S_k_samp)],
      M_k = setdiff(1:p, id_map[rownames(S_k_samp)]),
      nu = nu_k,
      alpha_k = 1.0
    )
  })

  list(
    p = p, K = length(S_list), all_ids = all_ids, id_map = id_map,
    samples = setNames(samples, sample_ids), scale_method = scale_method,
    alpha_normalization = alpha_normalization
  )
}

.initialize_sigma <- function(internal_data, method) {
  p <- internal_data$p

  if (is.matrix(method)) {
    stopifnot(nrow(method) == p, ncol(method) == p)
    return(method)
  }

  if (method == "identity") {
    return(diag(p))
  }

  if (method == "avg_padded") {
    sum_mat <- count_mat <- matrix(0, p, p)

    for (s in internal_data$samples) {
      idx <- s$O_k
      # s$W_k is the Wishart matrix (nu * sample_cov)
      # Divide by nu to get back to sample covariance scale for initialization
      S_k_unname <- unname(s$W_k) / s$nu
      sum_mat[idx, idx] <- sum_mat[idx, idx] + S_k_unname
      count_mat[idx, idx] <- count_mat[idx, idx] + 1
    }

    avg_mat <- sum_mat / pmax(count_mat, 1)

    # Find a reasonable scale from the observed variances
    obs_vars <- diag(avg_mat)[diag(count_mat) > 0]
    mean_obs_var <- if (length(obs_vars) > 0) mean(obs_vars, na.rm = TRUE) else 1.0

    # Create a full-rank "filler" matrix (scaled identity)
    filler_mat <- diag(mean_obs_var, p)

    # Find all entries (diag and off-diag) that were never observed
    unobs_idx <- which(count_mat == 0)

    if (length(unobs_idx) > 0) {
      # Replace all zero-count entries with the filler values
      avg_mat[unobs_idx] <- filler_mat[unobs_idx]
    }

    avg_mat <- .symmetrize(avg_mat)

    # Check for unobserved variables and warn
    unobs <- which(diag(count_mat) == 0)
    if (length(unobs) > 0) {
      warning(
        sprintf(
          "Variables %s are never observed in any sample. Their covariance parameters will not be identifiable. Consider removing these variables or providing additional samples.",
          paste(internal_data$all_ids[unobs], collapse = ", ")
        ),
        call. = FALSE
      )
    }

    dimnames(avg_mat) <- list(internal_data$all_ids, internal_data$all_ids)
    return(avg_mat)
  }

  stop("Invalid initialization method")
}

.e_step_k <- function(sigma, s, ridge) {
  p <- nrow(sigma)
  O_k <- s$O_k
  M_k <- s$M_k
  W_k <- s$W_k # Wishart matrix for observed block
  nu_k <- s$nu
  if (is.null(nu_k) || !is.finite(nu_k)) {
    stop("Sample information must include a finite 'nu' value.", call. = FALSE)
  }
  alpha_k <- if (is.null(s$alpha_k)) 1 else s$alpha_k

  mk <- length(M_k)

  if (mk == 0L) {
    W_tilde <- matrix(0, p, p)
    W_tilde[O_k, O_k] <- W_k
    return(.symmetrize(W_tilde))
  }

  Sigma_OO <- sigma[O_k, O_k, drop = FALSE]
  Sigma_MO <- sigma[M_k, O_k, drop = FALSE]
  Sigma_MM <- sigma[M_k, M_k, drop = FALSE]

  Sigma_OO_inv <- .safe_chol_inverse(Sigma_OO, ridge = ridge)

  B_k <- Sigma_MO %*% Sigma_OO_inv
  Delta_k <- Sigma_MM - B_k %*% t(Sigma_MO)

  # Conditional Wishart distribution has degrees of freedom (nu_k - p_O)
  # where p_O = |O_k| is the number of observed variables
  # Reference: Anderson (2003), Theorem 7.3.4
  p_O <- length(O_k)

  # Validity check: conditional Wishart requires nu_k >= p_O
  if (nu_k < p_O) {
    stop(sprintf(
      "Sample degrees of freedom nu_k = %g must be >= |O_k| = %d for conditional Wishart distribution to be defined. Sample has insufficient degrees of freedom.",
      nu_k, p_O
    ), call. = FALSE)
  }

  E_W_OO <- W_k
  E_W_MO <- B_k %*% W_k
  E_W_MM <- (nu_k - p_O) * alpha_k * Delta_k + B_k %*% W_k %*% t(B_k)

  W_tilde <- matrix(0, p, p)
  W_tilde[O_k, O_k] <- E_W_OO
  W_tilde[M_k, O_k] <- E_W_MO
  W_tilde[O_k, M_k] <- t(E_W_MO)
  W_tilde[M_k, M_k] <- E_W_MM

  .symmetrize(W_tilde)
}

.m_step <- function(w_tilde_list, internal_data, sigma_old, ridge = 1e-8) {
  sample_ids <- names(internal_data$samples)
  nu_total <- sum(sapply(internal_data$samples, `[[`, "nu"))
  p <- internal_data$p

  alpha_vec <- vapply(sample_ids, function(id) internal_data$samples[[id]]$alpha_k, numeric(1))

  if (internal_data$scale_method == "estimate") {
    tol_alpha <- 1e-6
    tol_sigma <- 1e-6
    max_inner <- 25L
    sigma_prev <- sigma_old

    for (inner in seq_len(max_inner)) {
      w_sum_inner <- Reduce(`+`, lapply(sample_ids, function(id) {
        w_tilde_list[[id]] / alpha_vec[[id]]
      }))
      sigma_candidate <- w_sum_inner / nu_total

      # Symmetrize to correct numerical errors from floating-point arithmetic
      # (w_tilde matrices are already symmetric from E-step, but accumulation may introduce asymmetry)
      sigma_candidate <- .symmetrize(sigma_candidate)
      # Use safe inversion with progressive ridge regularization
      sigma_inv <- .safe_chol_inverse(sigma_candidate, ridge = ridge)

      # Alpha MLE formula (with gamma=1 fixed)
      # Formula: alpha_k = tr(Sigma^{-1} W) / (nu_k * p)
      # Reference: Complete-data log-likelihood derivative
      alpha_update <- vapply(sample_ids, function(id) {
        w_tilde <- w_tilde_list[[id]]
        nu_k <- internal_data$samples[[id]]$nu
        val <- sum(diag(sigma_inv %*% w_tilde)) / (nu_k * p)
        max(val, 1e-8)
      }, numeric(1))

      # Apply normalization constraint
      if (internal_data$alpha_normalization == "geometric") {
        # Geometric mean normalization: (∏ α_k)^(1/K) = 1
        # Implementation: α_k^norm = α_k / exp(mean(log(α_k)))
        K <- length(alpha_update)
        geom_mean <- exp(mean(log(alpha_update)))
        if (is.na(geom_mean) || geom_mean <= 0) {
          alpha_update[] <- 1
        } else {
          alpha_update <- alpha_update / geom_mean
        }
      } else if (internal_data$alpha_normalization == "arithmetic") {
        # Arithmetic mean normalization: (1/K) ∑ α_k = 1
        alpha_mean <- mean(alpha_update)
        if (is.na(alpha_mean) || alpha_mean <= 0) {
          alpha_update[] <- 1
        } else {
          alpha_update <- alpha_update / alpha_mean
        }
      }

      change_alpha <- max(abs(alpha_update - alpha_vec) / pmax(abs(alpha_vec), 1e-8))
      change_sigma <- norm(sigma_candidate - sigma_prev, "F") /
        pmax(norm(sigma_prev, "F"), .Machine$double.eps)

      alpha_vec <- alpha_update
      sigma_prev <- sigma_candidate

      if (change_alpha < tol_alpha && change_sigma < tol_sigma) {
        break
      }
    }

    for (id in sample_ids) {
      internal_data$samples[[id]]$alpha_k <- alpha_vec[[id]]
    }
  }

  w_sum <- Reduce(`+`, lapply(sample_ids, function(id) {
    w_tilde_list[[id]] / internal_data$samples[[id]]$alpha_k
  }))

  sigma_new <- w_sum / nu_total

  list(sigma_new = sigma_new, internal_data = internal_data)
}

.compute_loglik <- function(sigma, internal_data) {
  loglik <- sum(sapply(internal_data$samples, function(s) {
    W_k <- s$W_k # Wishart matrix for observed block
    O_k <- s$O_k
    nu_k <- s$nu
    pk <- length(O_k)
    alpha_k <- s$alpha_k

    Sigma_OO <- sigma[O_k, O_k, drop = FALSE]
    Sigma_OO_sym <- .symmetrize(Sigma_OO)
    # Use safe inversion with progressive ridge regularization
    # Get both inverse and Cholesky factor to avoid redundant decomposition
    chol_result <- .safe_chol_inverse(Sigma_OO_sym, ridge = 1e-10, return_chol = TRUE)
    Sigma_OO_inv <- chol_result$inverse
    chol_Sigma <- chol_result$chol
    log_det_Sigma <- 2 * sum(log(diag(chol_Sigma)))

    log_det_W <- determinant(W_k, logarithm = TRUE)$modulus
    tr_term <- sum(diag(Sigma_OO_inv %*% W_k))

    # Wishart log-likelihood with constant terms
    # Reference: Anderson (2003), Theorem 7.2.2
    # log p(W | Sigma, nu) = const + (nu - p - 1)/2 * log|W| - nu/2 * log|Sigma| - 1/2 * tr(Sigma^{-1} W)
    # where const = -nu*p/2 * log(2) - p(p-1)/4 * log(pi) - sum_{j=1}^p log Gamma((nu + 1 - j)/2)

    # Compute constant terms
    log_const <- -nu_k * pk / 2 * log(2) - pk * (pk - 1) / 4 * log(pi)
    for (j in 1:pk) {
      log_const <- log_const - lgamma((nu_k + 1 - j) / 2)
    }

    log_const +
      0.5 * (nu_k - pk - 1) * as.numeric(log_det_W) -
      0.5 * nu_k * (pk * log(alpha_k) + log_det_Sigma) -
      0.5 * (tr_term / alpha_k)
  }))
  as.numeric(loglik)
}

.compute_coverage <- function(internal_data) {
  p <- internal_data$p
  coverage <- matrix(0, p, p)
  for (s in internal_data$samples) {
    idx <- s$O_k
    if (length(idx) == 0L) next
    nu_k <- s$nu
    coverage[idx, idx] <- coverage[idx, idx] + nu_k
  }
  coverage
}

#' Check Graph Connectivity of Observation Pattern
#'
#' @description
#' Performs a connectivity analysis on the observation pattern to determine
#' whether all variables can be estimated jointly. Variables are connected if
#' they appear together in at least one observation set, either directly or
#' through a chain of observations.
#'
#' @param p Integer. Total number of variables.
#' @param observed_sets List of integer vectors. Each element contains the
#'   indices of variables observed together in one sample.
#'
#' @return A list with components:
#'   \item{is_connected}{Logical. TRUE if all variables form a single connected component.}
#'   \item{num_components}{Integer. Number of connected components.}
#'   \item{component_map}{Integer vector of length p. Maps each variable to its component number.}
#'
#' @details
#' This function constructs an undirected graph where:
#' \itemize{
#'   \item Vertices represent variables (1 to p)
#'   \item An edge exists between variables i and j if they appear together in any observation set
#' }
#'
#' The function uses breadth-first search (BFS) to identify connected components.
#' A connected graph (num_components == 1) is required for identifiable parameter estimation.
#'
#' Variables that never appear in any observation set form singleton components.
#'
#' @keywords internal
.check_graph_connectivity <- function(p, observed_sets) {
  # Input validation
  if (!is.numeric(p) || length(p) != 1L || p < 1L) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(p)

  if (!is.list(observed_sets)) {
    stop("observed_sets must be a list of integer vectors.", call. = FALSE)
  }

  # Initialize adjacency matrix (p x p)
  adjacency_matrix <- matrix(0L, nrow = p, ncol = p)

  # Build adjacency matrix from observed sets
  for (obs_set in observed_sets) {
    # Skip empty observation sets
    if (length(obs_set) == 0L) next

    # Add edges for all pairs in this observation set
    if (length(obs_set) >= 2L) {
      for (i in seq_along(obs_set)) {
        for (j in seq_along(obs_set)) {
          if (i != j) {
            adjacency_matrix[obs_set[i], obs_set[j]] <- 1L
            adjacency_matrix[obs_set[j], obs_set[i]] <- 1L
          }
        }
      }
    }
  }

  # Initialize BFS data structures
  num_components <- 0L
  visited <- logical(p) # All FALSE initially
  component_map <- integer(p)

  # Run BFS from each unvisited node
  for (start_node in seq_len(p)) {
    if (visited[start_node]) next

    # Start a new component
    num_components <- num_components + 1L

    # BFS queue (using a simple vector, front pointer for efficiency)
    queue <- integer(p)
    queue_front <- 1L
    queue_back <- 1L

    # Initialize with start node
    queue[queue_back] <- start_node
    queue_back <- queue_back + 1L
    visited[start_node] <- TRUE
    component_map[start_node] <- num_components

    # BFS traversal
    while (queue_front < queue_back) {
      # Dequeue
      current_node <- queue[queue_front]
      queue_front <- queue_front + 1L

      # Visit all unvisited neighbors
      for (neighbor in seq_len(p)) {
        if (adjacency_matrix[current_node, neighbor] == 1L && !visited[neighbor]) {
          visited[neighbor] <- TRUE
          component_map[neighbor] <- num_components

          # Enqueue
          queue[queue_back] <- neighbor
          queue_back <- queue_back + 1L
        }
      }
    }
  }

  # Return results
  list(
    is_connected = (num_components == 1L),
    num_components = num_components,
    component_map = component_map
  )
}

.symmetrize <- function(A) {
  (A + t(A)) / 2
}

.project_to_pd <- function(A, min_eigen = 1e-10) {
  if (nrow(A) == 0 || ncol(A) == 0) {
    return(A)
  }
  A_sym <- .symmetrize(A)
  dn <- dimnames(A_sym) # Save dimnames before eigen decomposition
  eig <- eigen(A_sym, symmetric = TRUE)
  eig$values <- pmax(eig$values, min_eigen)
  # Fix: Use diag(x, nrow=length(x)) to handle 1x1 matrices correctly
  # R's diag(scalar) interprets scalar as dimension, not diagonal value
  A_pd <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  if (!is.null(dn)) dimnames(A_pd) <- dn # Restore dimnames
  A_pd
}
