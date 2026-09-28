test_that("identifiability_report finds no unidentified pairs in a saturated design", {
  ids <- c("v1", "v2", "v3")
  S <- diag(3)
  dimnames(S) <- list(ids, ids)

  rep <- identifiability_report(list(a = S), n_factors = NULL)

  expect_s3_class(rep, "covcomb_identifiability")
  expect_equal(rep$n_unidentified_pairs, 0L)
  expect_equal(rep$n_total_pairs, 3)
  expect_equal(rep$fraction_unidentified, 0)
  expect_true(rep$is_connected)
  expect_length(rep$unobserved_variables, 0L)
})


test_that("a connected chain design still leaves pairs unidentified", {
  # The exact trap: connectivity holds, pairwise co-observation does not.
  mk <- function(ids) {
    m <- diag(length(ids))
    dimnames(m) <- list(ids, ids)
    m
  }
  S_list <- list(
    a = mk(c("v1", "v2", "v3")),
    b = mk(c("v3", "v4", "v5"))
  )

  rep <- identifiability_report(S_list, n_factors = NULL)

  # Connected through v3 ...
  expect_true(rep$is_connected)
  expect_equal(rep$n_components, 1L)

  # ... yet v1,v2 x v4,v5 are never co-observed: 4 pairs out of C(5,2)=10.
  expect_equal(rep$n_unidentified_pairs, 4L)
  expect_equal(rep$n_total_pairs, 10)
  expect_equal(rep$fraction_unidentified, 0.4)

  got <- paste(rep$unidentified_pairs$var1, rep$unidentified_pairs$var2)
  expect_setequal(got, c("v1 v4", "v1 v5", "v2 v4", "v2 v5"))
})


test_that("co_observation counts samples and orders variables like the fitter", {
  mk <- function(ids) {
    m <- diag(length(ids))
    dimnames(m) <- list(ids, ids)
    m
  }
  # Deliberately unsorted names, to pin the sorted-union ordering contract.
  S_list <- list(
    a = mk(c("b", "a")),
    b = mk(c("b", "c"))
  )

  rep <- identifiability_report(S_list)

  expect_equal(rep$variables, c("a", "b", "c"))
  expect_equal(rep$co_observation["a", "b"], 1L)
  expect_equal(rep$co_observation["b", "c"], 1L)
  expect_equal(rep$co_observation["a", "c"], 0L)
  expect_equal(rep$co_observation["b", "b"], 2L)
  expect_equal(rep$n_unidentified_pairs, 1L)
})


test_that("variables absent from every sample are reported", {
  mk <- function(ids) {
    m <- diag(length(ids))
    dimnames(m) <- list(ids, ids)
    m
  }
  rep <- identifiability_report(list(a = mk(c("v1", "v2"))))
  expect_length(rep$unobserved_variables, 0L)

  # A variable can only be "unobserved" relative to the union, so build a
  # coverage matrix directly to exercise the branch.
  cov_mat <- matrix(0, 3, 3)
  cov_mat[1:2, 1:2] <- 5
  u <- .unidentified_from_coverage(cov_mat, c("v1", "v2", "v3"))
  expect_equal(u$n, 2L)
  expect_equal(u$total, 3)
})


test_that("identifiability_report rejects malformed input", {
  expect_error(identifiability_report(list()), "non-empty")
  expect_error(identifiability_report(matrix(1)), "non-empty|list")

  S <- diag(2) # no dimnames
  expect_error(identifiability_report(list(a = S)), "rownames")

  one <- matrix(1, 1, 1, dimnames = list("v1", "v1"))
  expect_error(identifiability_report(list(a = one)), "two distinct variables")
})


test_that(".k_max matches the documented saturation formula", {
  for (p in c(2, 5, 10, 37, 200)) {
    expect_equal(.k_max(p), as.integer(floor((2 * p + 1 - sqrt(8 * p + 1)) / 2)))
  }
})


test_that(".is_free_model recognises both free-Sigma encodings", {
  expect_true(.is_free_model(NULL))
  expect_true(.is_free_model("free"))
  expect_false(.is_free_model("auto"))
  expect_false(.is_free_model(2L))
})


test_that("fit_covcomb warns on unidentified pairs under the free model", {
  set.seed(42)
  mk <- function(ids) {
    m <- diag(length(ids)) + 0.3
    dimnames(m) <- list(ids, ids)
    m
  }
  S_list <- list(
    a = mk(c("v1", "v2", "v3")),
    b = mk(c("v3", "v4", "v5"))
  )
  nu <- c(a = 50, b = 50)

  expect_warning(
    fit <- fit_covcomb(S_list, nu, n_factors = NULL, se_method = "none"),
    "never jointly observed"
  )

  expect_equal(fit$identifiability$n_unidentified_pairs, 4L)
  expect_equal(fit$identifiability$model, "free")
  expect_true(is.matrix(fit$identifiability$never_mask))
  # The mask must align with the fitted matrix it describes.
  expect_equal(dim(fit$identifiability$never_mask), dim(fit$Sigma_hat))
})


test_that("fit_covcomb stays quiet when every pair is jointly observed", {
  set.seed(42)
  ids <- c("v1", "v2", "v3")
  mk <- function() {
    m <- diag(3) + 0.3
    dimnames(m) <- list(ids, ids)
    m
  }
  S_list <- list(a = mk(), b = mk())
  nu <- c(a = 50, b = 50)

  expect_no_warning(
    fit <- fit_covcomb(S_list, nu, n_factors = NULL, se_method = "none")
  )
  expect_equal(fit$identifiability$n_unidentified_pairs, 0L)
})


test_that("the warning is the flat-likelihood claim, verified by initialization", {
  # The substantive point: on unidentified entries the fit tracks init_sigma,
  # while identified entries agree across initializations. If this ever stops
  # holding, the warning text is wrong and must change.
  set.seed(7)
  mk <- function(ids, r) {
    m <- diag(length(ids)) + r
    dimnames(m) <- list(ids, ids)
    m
  }
  S_list <- list(
    a = mk(c("v1", "v2", "v3"), 0.4),
    b = mk(c("v3", "v4", "v5"), 0.4)
  )
  nu <- c(a = 100, b = 100)

  fit_from <- function(seed) {
    set.seed(seed)
    p <- 5
    init <- diag(p) + 0.5 * seed / 10
    dimnames(init) <- list(
      c("v1", "v2", "v3", "v4", "v5"),
      c("v1", "v2", "v3", "v4", "v5")
    )
    suppressWarnings(
      fit_covcomb(S_list, nu,
        n_factors = NULL, se_method = "none",
        init_sigma = init
      )
    )
  }

  f1 <- fit_from(1)
  f2 <- fit_from(9)

  mask <- f1$identifiability$never_mask
  expect_true(any(mask))

  # Identified off-diagonal entries: agree closely across initializations.
  ident <- !mask
  diag(ident) <- FALSE
  expect_lt(
    max(abs(f1$Sigma_hat[ident] - f2$Sigma_hat[ident])),
    1e-4
  )

  # Unidentified entries: the whole point is that they need not agree.
  # Assert only that they are recorded, not that they differ by any fixed
  # amount, since that depends on how far apart the initializations are.
  expect_equal(sum(mask) / 2, f1$identifiability$n_unidentified_pairs)
})


test_that("print method reports both regimes without error", {
  mk <- function(ids) {
    m <- diag(length(ids))
    dimnames(m) <- list(ids, ids)
    m
  }
  S_list <- list(
    a = mk(c("v1", "v2", "v3")),
    b = mk(c("v3", "v4", "v5"))
  )

  free_out <- capture.output(print(identifiability_report(S_list, n_factors = NULL)))
  expect_true(any(grepl("free Sigma", free_out)))
  expect_true(any(grepl("init_sigma", free_out)))
  expect_true(any(grepl("4 of 10", free_out)))

  fa_out <- capture.output(print(identifiability_report(S_list, n_factors = 2L)))
  expect_true(any(grepl("factor", fa_out)))
  expect_true(any(grepl("rest on that assumption", fa_out)))

  clean <- capture.output(print(identifiability_report(list(a = mk(c("v1", "v2"))))))
  expect_true(any(grepl("Every variable pair", clean)))
})


test_that("print truncates long pair lists", {
  mk <- function(ids) {
    m <- diag(length(ids))
    dimnames(m) <- list(ids, ids)
    m
  }
  S_list <- list(
    a = mk(paste0("v", 1:6)),
    b = mk(paste0("v", 6:12))
  )
  rep <- identifiability_report(S_list, n_factors = NULL)
  expect_gt(rep$n_unidentified_pairs, 10L)

  out <- capture.output(print(rep, max_pairs = 3L))
  expect_true(any(grepl("and \\d+ more", out)))
})
