#' Report Which Covariance Entries a Design Can Identify
#'
#' @description
#' Determines, from the observation design alone and without fitting anything,
#' which entries of the combined covariance matrix are determined by the data
#' and which are not.
#'
#' A pair of variables \eqn{(i,j)} that is never jointly observed in any input
#' sample carries no information about \eqn{\Sigma_{ij}}. Under the
#' free-\eqn{\Sigma} model the observed-data log-likelihood is exactly flat in
#' that coordinate: the Fisher information is zero, no unbiased estimator
#' exists, and the value the EM algorithm returns there is a deterministic
#' function of the starting value \code{init_sigma} rather than of the data.
#' Different initializations agree to machine precision on every jointly
#' observed entry and disagree on these.
#'
#' Structured models supply an identifying assumption that can pin such entries
#' down. The k-factor model does exactly this, provided the variables share
#' factors through overlapping variables. That is a real resolution, not a
#' cosmetic one, but the resulting value rests on the factor assumption rather
#' than on direct observation, and should be reported that way.
#'
#' @param S_list Named list of covariance submatrices, as passed to
#'   \code{\link{fit_covcomb}}. Only the \code{rownames} are used.
#' @param n_factors The model that will be fitted, using the same encoding as
#'   \code{\link{fit_covcomb}}: \code{"auto"} or an integer \eqn{k \ge 1} for a
#'   factor model, \code{NULL} or \code{"free"} for the unconstrained model.
#'   Affects only the advice printed, never the design facts computed.
#'
#' @return An object of class \code{covcomb_identifiability}, a list with:
#'   \item{variables}{Character vector of all variables, in the order used by
#'     \code{fit_covcomb} (the sorted union of \code{rownames}).}
#'   \item{p}{Number of variables.}
#'   \item{n_samples}{Number of input samples.}
#'   \item{co_observation}{Integer matrix; entry \eqn{(i,j)} counts the samples
#'     observing both \eqn{i} and \eqn{j}.}
#'   \item{never_mask}{Logical matrix, \code{TRUE} where an off-diagonal pair is
#'     never jointly observed.}
#'   \item{unidentified_pairs}{Data frame of those pairs.}
#'   \item{n_unidentified_pairs, n_total_pairs, fraction_unidentified}{Counts and
#'     the proportion of off-diagonal entries not determined by the data.}
#'   \item{unobserved_variables}{Variables absent from every sample.}
#'   \item{is_connected, n_components}{Connectivity of the co-observation graph.
#'     Note that connectivity is strictly weaker than pairwise co-observation:
#'     a fully connected design can still leave most pairs unidentified.}
#'   \item{k_max}{Saturation point \eqn{\lfloor(2p+1-\sqrt{8p+1})/2\rfloor},
#'     beyond which a factor model adds no parameters relative to free-\eqn{\Sigma}.}
#'   \item{n_factors}{The model encoding that was passed in.}
#'
#' @details
#' One boundary case is worth knowing. If either conditional variance entering
#' the feasible interval for \eqn{\Sigma_{ij}} is exactly zero, the entry is
#' identified after all, forced by the positive-definiteness constraint alone.
#' That situation is not generic and this function does not attempt to detect
#' it; a pair reported here as unidentified is unidentified except in that
#' degenerate case.
#'
#' Connectivity of the co-observation graph, which
#' \code{fit_covcomb} already checks, is a much weaker requirement. A chain
#' design in which sample 1 observes variables 1-8, sample 2 observes 6-15 and
#' sample 3 observes 13-20 is fully connected, yet every pair spanning the first
#' and third blocks is never jointly observed: half the off-diagonal entries.
#'
#' @examples
#' set.seed(1)
#' mk <- function(ids) {
#'   m <- diag(length(ids))
#'   dimnames(m) <- list(ids, ids)
#'   m
#' }
#' # A connected chain design that nonetheless leaves pairs unidentified.
#' S_list <- list(
#'   a = mk(c("v1", "v2", "v3")),
#'   b = mk(c("v3", "v4", "v5"))
#' )
#' identifiability_report(S_list, n_factors = NULL)
#'
#' @seealso \code{\link{fit_covcomb}}
#' @export
identifiability_report <- function(S_list, n_factors = "auto") {
  if (!is.list(S_list) || length(S_list) == 0L) {
    stop("S_list must be a non-empty list of covariance matrices.", call. = FALSE)
  }

  id_sets <- lapply(S_list, function(S_k) {
    ids <- rownames(S_k)
    if (is.null(ids)) {
      stop(
        "Every element of S_list must have rownames identifying its variables.",
        call. = FALSE
      )
    }
    unique(ids)
  })

  # Same variable universe and ordering that .preprocess_data() builds, so that
  # indices here line up with the rows/columns of a fitted Sigma_hat.
  all_ids <- sort(unique(unlist(id_sets)))
  p <- length(all_ids)
  if (p < 2L) {
    stop("At least two distinct variables are required.", call. = FALSE)
  }
  id_map <- stats::setNames(seq_len(p), all_ids)

  co_obs <- matrix(0L, p, p, dimnames = list(all_ids, all_ids))
  for (ids in id_sets) {
    idx <- id_map[ids]
    co_obs[idx, idx] <- co_obs[idx, idx] + 1L
  }

  never_mask <- co_obs == 0L
  diag(never_mask) <- FALSE

  # Variables absent from every sample: their diagonal was never incremented.
  unobserved <- all_ids[diag(co_obs) == 0L]

  upper <- upper.tri(never_mask)
  hits <- which(never_mask & upper, arr.ind = TRUE)
  unidentified_pairs <- data.frame(
    var1 = all_ids[hits[, "row"]],
    var2 = all_ids[hits[, "col"]],
    stringsAsFactors = FALSE
  )
  if (nrow(unidentified_pairs) > 1L) {
    ord <- order(unidentified_pairs$var1, unidentified_pairs$var2)
    unidentified_pairs <- unidentified_pairs[ord, , drop = FALSE]
    rownames(unidentified_pairs) <- NULL
  }

  n_total_pairs <- p * (p - 1L) / 2L
  n_unid <- nrow(unidentified_pairs)

  observed_sets <- lapply(id_sets, function(ids) unname(id_map[ids]))
  connectivity <- .check_graph_connectivity(p, observed_sets)

  structure(
    list(
      variables = all_ids,
      p = p,
      n_samples = length(S_list),
      co_observation = co_obs,
      never_mask = never_mask,
      unidentified_pairs = unidentified_pairs,
      n_unidentified_pairs = n_unid,
      n_total_pairs = n_total_pairs,
      fraction_unidentified = if (n_total_pairs > 0L) n_unid / n_total_pairs else 0,
      unobserved_variables = unobserved,
      is_connected = connectivity$is_connected,
      n_components = connectivity$num_components,
      k_max = .k_max(p),
      n_factors = n_factors
    ),
    class = "covcomb_identifiability"
  )
}


#' Print an Identifiability Report
#'
#' @param x An object of class \code{covcomb_identifiability}.
#' @param max_pairs Maximum number of unidentified pairs to list individually.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print covcomb_identifiability
print.covcomb_identifiability <- function(x, max_pairs = 10L, ...) {
  cat("CovCombR identifiability report\n")
  cat(sprintf(
    "Variables: %d   Samples: %d   Off-diagonal entries: %d\n",
    x$p, x$n_samples, x$n_total_pairs
  ))

  cat(sprintf(
    "Co-observation graph: %s (%d component%s)\n",
    if (x$is_connected) "connected" else "DISCONNECTED",
    x$n_components,
    if (x$n_components == 1L) "" else "s"
  ))

  if (length(x$unobserved_variables) > 0L) {
    cat(sprintf(
      "\nVariables never observed in any sample (%d): %s\n",
      length(x$unobserved_variables),
      paste(x$unobserved_variables, collapse = ", ")
    ))
  }

  if (x$n_unidentified_pairs == 0L) {
    cat("\nEvery variable pair is jointly observed in at least one sample.\n")
    cat("All entries are determined by the data under either model.\n")
    return(invisible(x))
  }

  cat(sprintf(
    "\nPairs NEVER jointly observed: %d of %d (%.1f%%)\n",
    x$n_unidentified_pairs, x$n_total_pairs, 100 * x$fraction_unidentified
  ))

  shown <- utils::head(x$unidentified_pairs, max_pairs)
  cat(paste0(
    "  ", shown$var1, " ~ ", shown$var2,
    collapse = "\n"
  ), "\n", sep = "")
  if (x$n_unidentified_pairs > nrow(shown)) {
    cat(sprintf(
      "  ... and %d more (see $unidentified_pairs)\n",
      x$n_unidentified_pairs - nrow(shown)
    ))
  }

  if (.is_free_model(x$n_factors)) {
    cat(
      "\nModel: free Sigma. The log-likelihood is exactly flat in these ",
      "entries.\nTheir fitted values are determined by init_sigma, not by the ",
      "data, and\nwill change if you change the starting value. Do not ",
      "interpret them, and\ndo not read conditional independence off the ",
      "corresponding precision entries.\n",
      sep = ""
    )
    cat(sprintf(
      "\nTo identify them, fit a factor model: fit_covcomb(..., n_factors = \"auto\")\n(k_max here is %d), or add a sample observing the pairs above.\n",
      x$k_max
    ))
  } else {
    cat(sprintf(
      "\nModel: %s factor. These entries are identified through the factor\nstructure, not by direct observation, so they rest on that assumption.\nReport them as such, and where the design has redundancy, check the\nassumption rather than trusting it.\n",
      if (identical(x$n_factors, "auto")) "auto-selected" else paste0("k = ", x$n_factors)
    ))
  }

  invisible(x)
}


#' Saturation Point for the Number of Factors
#'
#' @param p Number of variables.
#' @return Largest k beyond which a k-factor model adds no parameters relative
#'   to the free-\eqn{\Sigma} model.
#' @keywords internal
.k_max <- function(p) {
  as.integer(floor((2 * p + 1 - sqrt(8 * p + 1)) / 2))
}


#' Is This n_factors Encoding the Free-Sigma Model?
#'
#' @param n_factors Value using \code{fit_covcomb}'s encoding.
#' @return \code{TRUE} for \code{NULL} and \code{"free"}.
#' @keywords internal
.is_free_model <- function(n_factors) {
  is.null(n_factors) || identical(n_factors, "free")
}


#' Summarize Unidentified Pairs From a Coverage Matrix
#'
#' @description
#' Internal counterpart to \code{\link{identifiability_report}} for use inside
#' \code{fit_covcomb}, where the nu-weighted coverage matrix has already been
#' computed. Since every \code{nu_k} is required to be positive, a zero
#' off-diagonal coverage entry is exactly a never-jointly-observed pair.
#'
#' @param coverage_mat Coverage matrix from \code{.compute_coverage}.
#' @param all_ids Character vector of variable names, in matrix order.
#' @return List with \code{n}, \code{total}, \code{fraction} and \code{pairs}.
#' @keywords internal
.unidentified_from_coverage <- function(coverage_mat, all_ids) {
  p <- nrow(coverage_mat)
  never <- coverage_mat == 0
  diag(never) <- FALSE
  hits <- which(never & upper.tri(never), arr.ind = TRUE)

  pairs <- data.frame(
    var1 = all_ids[hits[, "row"]],
    var2 = all_ids[hits[, "col"]],
    stringsAsFactors = FALSE
  )
  total <- p * (p - 1L) / 2L

  list(
    n = nrow(pairs),
    total = total,
    fraction = if (total > 0L) nrow(pairs) / total else 0,
    pairs = pairs,
    mask = never
  )
}


#' Format a Few Example Pairs for a Warning Message
#'
#' @param pairs Data frame with \code{var1} and \code{var2}.
#' @param n_show Number of pairs to name.
#' @return A single string.
#' @keywords internal
.format_pairs <- function(pairs, n_show = 3L) {
  shown <- utils::head(pairs, n_show)
  txt <- paste0(shown$var1, "~", shown$var2, collapse = ", ")
  if (nrow(pairs) > nrow(shown)) {
    txt <- paste0(txt, ", ...")
  }
  txt
}
