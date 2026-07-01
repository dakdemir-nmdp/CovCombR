#' Fit Factor-Model Wishart EM
#'
#' @description
#' Estimates a k-factor covariance model Sigma = Lambda Lambda^T + Psi
#' (Psi diagonal) from incomplete sample covariance matrices using a
#' Wishart EM algorithm with Rubin-Thayer factor-analytic M-step.
#'
#' @param S_list Named list of sample covariance matrices.
#' @param nu Named numeric vector of degrees of freedom.
#' @param k Integer number of latent factors.
#' @param scale_method Scaling method: \code{"none"} (default) or \code{"estimate"}.
#' @param alpha_normalization Normalization for scale factors: \code{"geometric"}
#'   (default) or \code{"arithmetic"}.
#' @param init_sigma Initial covariance matrix or method (\code{"identity"} or
#'   \code{"avg_padded"}).
#' @param control List of control parameters (max_iter, tol, ridge, min_eigen).
#'
#' @return An object of class \code{"covcomb"} with the same fields as
#'   \code{fit_covcomb()} plus:
#'   \itemize{
#'     \item \code{Lambda_hat}: p x k factor loadings matrix (varimax-rotated when k > 1)
#'     \item \code{Psi_hat}: length-p vector of unique variances
#'     \item \code{n_factors}: k
#'     \item \code{model}: \code{"factor"}
#'   }
#'
#' @export
fit_fa_em <- function(S_list, nu, k,
                      scale_method = "none",
                      alpha_normalization = "geometric",
                      init_sigma = "identity",
                      control = list()) {
  call <- match.call()

  scale_method       <- match.arg(scale_method, c("none", "estimate"))
  alpha_normalization <- match.arg(alpha_normalization, c("geometric", "arithmetic"))

  # Validate inputs
  stopifnot(is.list(S_list), !is.null(names(S_list)))
  stopifnot(is.numeric(nu), !is.null(names(nu)))
  stopifnot(all(names(S_list) %in% names(nu)))
  nu <- nu[names(S_list)]
  if (any(is.na(nu))) {
    stop("Degrees of freedom vector must include entries for all samples in S_list.",
         call. = FALSE)
  }

  k <- as.integer(k)
  if (k < 1L) stop("k must be >= 1", call. = FALSE)

  # Set controls
  ctrl <- list(max_iter = 500L, tol = 1e-7, ridge = 1e-8, min_eigen = 1e-10)
  ctrl[names(control)] <- control

  # Pre-process data
  internal_data <- .preprocess_data(S_list, nu, scale_method, alpha_normalization)
  p <- internal_data$p

  if (k > p) {
    stop(sprintf("k (%d) cannot exceed p (%d).", k, p), call. = FALSE)
  }

  # Check graph connectivity (warning only, same as fit_covcomb)
  observed_sets <- lapply(internal_data$samples, function(s) s$O_k)
  connectivity  <- .check_graph_connectivity(p, observed_sets)
  if (!connectivity$is_connected) {
    warning(
      sprintf(
        "The observation pattern creates %d disconnected component%s. ",
        connectivity$num_components,
        if (connectivity$num_components > 1L) "s" else ""
      ),
      "Variables in different components cannot be jointly estimated. ",
      "Consider using only samples that observe overlapping variable sets.",
      call. = FALSE
    )
  }

  # Initialize full covariance Sigma from which we derive Lambda, Psi.
  # For "identity" init we use avg_padded to get sensible starting loadings,
  # then run the outer EM from the factor reconstruction of that avg matrix.
  # This avoids the degenerate case where identity's eigenvectors point along
  # coordinate axes, creating Psi=0 for one variable on the first iteration.
  init_for_loadings <- if (identical(init_sigma, "identity")) "avg_padded" else init_sigma
  sigma_init <- .initialize_sigma(internal_data, init_for_loadings)
  sigma_init  <- .project_to_pd(sigma_init, ctrl$min_eigen)

  # Derive initial Lambda and Psi from sigma_init via eigen-decomposition
  eig_init    <- eigen(sigma_init, symmetric = TRUE)
  # Use top-k eigenvectors scaled by sqrt(eigenvalue) as initial loadings
  Lambda_cur  <- eig_init$vectors[, seq_len(k), drop = FALSE] %*%
                   diag(sqrt(pmax(eig_init$values[seq_len(k)], 0)),
                        nrow = k)
  # Clamp Psi: each diagonal element must be at least 0.01 * diagonal of sigma_init
  # to avoid degenerate near-zero unique variances at initialization.
  min_psi_init <- pmax(0.01 * diag(sigma_init), 1e-8)
  Psi_cur      <- pmax(diag(sigma_init) - rowSums(Lambda_cur^2), min_psi_init)

  # Reconstruct starting Sigma from factor structure
  sigma_current <- Lambda_cur %*% t(Lambda_cur) + diag(Psi_cur, nrow = p)
  sigma_current <- .project_to_pd(sigma_current, ctrl$min_eigen)

  # EM loop
  history <- data.frame(
    iteration            = integer(),
    rel_change           = numeric(),
    log_likelihood       = numeric(),
    log_likelihood_per_df = numeric()
  )
  converged      <- FALSE
  rel_change     <- NA_real_
  iterations_used <- 0L
  total_nu       <- sum(sapply(internal_data$samples, function(s) s$nu))

  for (i in seq_len(ctrl$max_iter)) {

    # ---- E-step (unchanged: uses full p x p Sigma) ----
    w_tilde_list <- lapply(internal_data$samples, function(s) {
      .e_step_k(sigma_current, s, ctrl$ridge)
    })

    # ---- Build sufficient statistic T ----
    # Optionally update alpha_k values first (mirrors .m_step for scale_method="estimate")
    sample_ids <- names(internal_data$samples)
    alpha_vec  <- vapply(sample_ids,
                         function(id) internal_data$samples[[id]]$alpha_k,
                         numeric(1L))

    if (internal_data$scale_method == "estimate") {
      tol_alpha  <- 1e-6
      tol_sigma_inner <- 1e-6
      max_inner  <- 25L
      sigma_prev_inner <- sigma_current

      for (inner in seq_len(max_inner)) {
        w_sum_inner <- Reduce(`+`, lapply(sample_ids, function(id) {
          w_tilde_list[[id]] / alpha_vec[[id]]
        }))
        sigma_candidate <- .symmetrize(w_sum_inner / total_nu)
        sigma_inv_inner  <- .safe_chol_inverse(sigma_candidate, ridge = ctrl$ridge)

        alpha_update <- vapply(sample_ids, function(id) {
          w_tilde <- w_tilde_list[[id]]
          nu_k    <- internal_data$samples[[id]]$nu
          val     <- sum(diag(sigma_inv_inner %*% w_tilde)) / (nu_k * p)
          max(val, 1e-8)
        }, numeric(1L))

        if (internal_data$alpha_normalization == "geometric") {
          geom_mean <- exp(mean(log(alpha_update)))
          if (is.na(geom_mean) || geom_mean <= 0) {
            alpha_update[] <- 1
          } else {
            alpha_update <- alpha_update / geom_mean
          }
        } else {
          alpha_mean <- mean(alpha_update)
          if (is.na(alpha_mean) || alpha_mean <= 0) {
            alpha_update[] <- 1
          } else {
            alpha_update <- alpha_update / alpha_mean
          }
        }

        change_alpha <- max(abs(alpha_update - alpha_vec) /
                              pmax(abs(alpha_vec), 1e-8))
        change_sigma_inner <- norm(sigma_candidate - sigma_prev_inner, "F") /
          max(norm(sigma_prev_inner, "F"), .Machine$double.eps)

        alpha_vec         <- alpha_update
        sigma_prev_inner  <- sigma_candidate

        if (change_alpha < tol_alpha && change_sigma_inner < tol_sigma_inner) break
      }

      for (id in sample_ids) {
        internal_data$samples[[id]]$alpha_k <- alpha_vec[[id]]
      }
    }

    # Weighted sum: T = (sum_k W_tilde_k / alpha_k) / sum_k nu_k
    w_sum <- Reduce(`+`, lapply(sample_ids, function(id) {
      w_tilde_list[[id]] / internal_data$samples[[id]]$alpha_k
    }))
    T_mat <- w_sum / total_nu  # p x p sufficient statistic

    # ---- FA M-step: Rubin-Thayer 1982 inner iterations ----
    max_fa_inner <- 50L
    tol_fa       <- ctrl$tol

    for (fa_iter in seq_len(max_fa_inner)) {
      Lambda_old <- Lambda_cur
      Psi_old    <- Psi_cur

      Psi_mat    <- diag(Psi_cur, nrow = p)
      Psi_inv    <- diag(1 / Psi_cur, nrow = p)

      # Phi (k x p): factor regression matrix
      # Phi = solve(t(Lambda) %*% Psi^{-1} %*% Lambda + I_k) %*% t(Lambda) %*% Psi^{-1}
      LtPsiInv   <- t(Lambda_cur) %*% Psi_inv        # k x p
      LtPsiInvL  <- LtPsiInv %*% Lambda_cur          # k x k
      Phi        <- .safe_chol_inverse(LtPsiInvL + diag(k),
                                        ridge = ctrl$ridge) %*% LtPsiInv  # k x p

      # Lambda_new = T %*% t(Phi) %*% solve(Phi %*% T %*% t(Phi) + I_k - Phi %*% Lambda)
      PhiT       <- Phi %*% T_mat                     # k x p
      PhiTPhi_t  <- PhiT %*% t(Phi)                   # k x k
      rhs_mat    <- PhiTPhi_t + diag(k) - Phi %*% Lambda_cur   # k x k
      Lambda_new <- T_mat %*% t(Phi) %*% .safe_chol_inverse(rhs_mat,
                                                              ridge = ctrl$ridge)

      # Psi_new = diag(T - Lambda_new %*% Phi %*% T)
      Psi_new    <- pmax(diag(T_mat - Lambda_new %*% PhiT), 1e-8)

      Lambda_cur <- Lambda_new
      Psi_cur    <- Psi_new

      # Check inner convergence
      fa_change <- max(
        max(abs(Lambda_cur - Lambda_old)) / max(max(abs(Lambda_old)), 1e-10),
        max(abs(Psi_cur    - Psi_old))    / max(max(abs(Psi_old)),    1e-10)
      )
      if (fa_change < tol_fa) break
    }

    # Reconstruct Sigma from factor model
    sigma_new <- Lambda_cur %*% t(Lambda_cur) + diag(Psi_cur, nrow = p)
    sigma_new <- .project_to_pd(sigma_new, ctrl$min_eigen)

    # Monitor convergence on Sigma (same criterion as fit_covcomb)
    denom      <- max(norm(sigma_current, "F"), .Machine$double.eps)
    rel_change <- norm(sigma_new - sigma_current, "F") / denom

    loglik         <- .compute_loglik(sigma_new, internal_data)
    loglik_per_df  <- loglik / total_nu
    history[i, ]   <- list(
      iteration            = i,
      rel_change           = rel_change,
      log_likelihood       = loglik,
      log_likelihood_per_df = loglik_per_df
    )

    sigma_current <- sigma_new
    iterations_used <- i

    if (rel_change < ctrl$tol) {
      converged <- TRUE
      break
    }
  }

  # Apply varimax rotation to Lambda for identifiability.
  # Guard against all-zero rows (variables never observed in any study), which
  # cause varimax to produce NaN via 0/0 in its internal row-normalization step.
  Lambda_hat <- Lambda_cur
  if (k > 1L) {
    row_norms  <- sqrt(rowSums(Lambda_hat^2))
    active     <- row_norms > .Machine$double.eps
    if (any(active)) {
      rot              <- stats::varimax(Lambda_hat[active, , drop = FALSE])
      Lambda_hat[active, ] <- rot$loadings[]
    }
  }
  Psi_hat <- Psi_cur

  # Reconstruct final Sigma (after possible varimax rotation; rotation preserves LL'=LL')
  Sigma_hat <- Lambda_hat %*% t(Lambda_hat) + diag(Psi_hat, nrow = p)
  dimnames(Sigma_hat) <- list(internal_data$all_ids, internal_data$all_ids)
  rownames(Lambda_hat) <- internal_data$all_ids

  # Compute alpha_hat
  sample_names <- names(S_list)
  if (is.null(sample_names)) sample_names <- as.character(seq_along(S_list))
  alpha_hat <- if (scale_method == "estimate") {
    setNames(sapply(internal_data$samples, `[[`, "alpha_k"), sample_names)
  } else {
    setNames(rep(1, length(S_list)), sample_names)
  }

  # S_hat and S_hat_scale
  if (scale_method == "estimate") {
    if (alpha_normalization == "geometric") {
      S_hat_scale <- 1
    } else {
      S_hat_scale <- mean(alpha_hat)
    }
    S_hat <- S_hat_scale * Sigma_hat
  } else {
    S_hat_scale <- 1
    S_hat       <- Sigma_hat
  }
  dimnames(S_hat) <- dimnames(Sigma_hat)

  result <- list(
    Sigma_hat  = Sigma_hat,
    S_hat      = S_hat,
    S_hat_scale = S_hat_scale,
    alpha_hat  = alpha_hat,
    scale_method = scale_method,
    se_method  = "none",
    Sigma_se   = NULL,
    S_hat_se   = NULL,
    bootstrap  = NULL,
    bootstrap_samples = NULL,
    convergence = list(
      converged    = converged,
      iterations   = iterations_used,
      final_rel_change = rel_change
    ),
    history    = history,
    # Factor model extras
    Lambda_hat = Lambda_hat,
    Psi_hat    = Psi_hat,
    n_factors  = k,
    model      = "factor",
    call       = call
  )

  structure(result, class = "covcomb")
}
