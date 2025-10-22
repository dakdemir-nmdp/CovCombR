test_that("fully connected graph with overlapping observations", {
  # Test case: 4 variables, all connected through observation patterns
  # Sample 1: variables 1,2,3
  # Sample 2: variables 2,3,4
  # This creates a connected graph: 1-2-3-4

  observed_sets <- list(
    c(1, 2, 3),
    c(2, 3, 4)
  )

  result <- .check_graph_connectivity(p = 4, observed_sets = observed_sets)

  expect_true(result$is_connected)
  expect_equal(result$num_components, 1)
  expect_length(result$component_map, 4)
  expect_true(all(result$component_map == 1))
})

test_that("disconnected graph with two components", {
  # Test case: 6 variables in 2 disconnected components
  # Component 1: variables 1,2,3 (connected through samples)
  # Component 2: variables 4,5,6 (connected through samples)
  # No sample observes variables from both components

  observed_sets <- list(
    c(1, 2, 3),    # Component 1
    c(1, 2),       # Component 1
    c(4, 5, 6),    # Component 2
    c(5, 6)        # Component 2
  )

  result <- .check_graph_connectivity(p = 6, observed_sets = observed_sets)

  expect_false(result$is_connected)
  expect_equal(result$num_components, 2)
  expect_length(result$component_map, 6)

  # Check that variables 1,2,3 are in one component
  # and variables 4,5,6 are in another component
  component_1 <- result$component_map[1]
  component_2 <- result$component_map[4]

  expect_equal(result$component_map[1], component_1)
  expect_equal(result$component_map[2], component_1)
  expect_equal(result$component_map[3], component_1)

  expect_equal(result$component_map[4], component_2)
  expect_equal(result$component_map[5], component_2)
  expect_equal(result$component_map[6], component_2)

  expect_false(component_1 == component_2)
})

test_that("graph with isolated variables", {
  # Test case: 5 variables
  # Variables 1,2,3 are connected
  # Variable 4 is isolated (observed but never with other variables)
  # Variable 5 is unobserved

  observed_sets <- list(
    c(1, 2, 3),
    c(1, 2),
    c(4)        # isolated variable
  )

  result <- .check_graph_connectivity(p = 5, observed_sets = observed_sets)

  expect_false(result$is_connected)
  expect_equal(result$num_components, 3)  # {1,2,3}, {4}, {5}
  expect_length(result$component_map, 5)

  # Variable 5 should be in its own component (unobserved)
  # Variable 4 should be in its own component (isolated)
  # Variables 1,2,3 should be in the same component

  component_main <- result$component_map[1]
  expect_equal(result$component_map[2], component_main)
  expect_equal(result$component_map[3], component_main)
  expect_false(result$component_map[4] == component_main)
  expect_false(result$component_map[5] == component_main)
  expect_false(result$component_map[4] == result$component_map[5])
})

test_that("all variables unobserved", {
  # Edge case: no observations at all

  observed_sets <- list()

  result <- .check_graph_connectivity(p = 3, observed_sets = observed_sets)

  expect_false(result$is_connected)
  expect_equal(result$num_components, 3)
  expect_length(result$component_map, 3)

  # Each variable should be in its own component
  expect_equal(result$component_map, c(1, 2, 3))
})

test_that("single variable observed", {
  # Edge case: only one variable in the system

  observed_sets <- list(c(1))

  result <- .check_graph_connectivity(p = 1, observed_sets = observed_sets)

  expect_true(result$is_connected)
  expect_equal(result$num_components, 1)
  expect_equal(result$component_map, 1)
})

test_that("empty observation sets are ignored", {
  # Edge case: some samples have no observations

  observed_sets <- list(
    c(1, 2),
    integer(0),  # empty observation
    c(2, 3),
    integer(0)   # empty observation
  )

  result <- .check_graph_connectivity(p = 3, observed_sets = observed_sets)

  expect_true(result$is_connected)
  expect_equal(result$num_components, 1)
})

test_that("three components with complex pattern", {
  # More complex disconnected graph
  # Component 1: variables 1,2,3,4
  # Component 2: variables 5,6
  # Component 3: variables 7,8,9

  observed_sets <- list(
    c(1, 2),
    c(2, 3, 4),
    c(1, 4),
    c(5, 6),
    c(5, 6),
    c(7, 8),
    c(8, 9),
    c(7, 9)
  )

  result <- .check_graph_connectivity(p = 9, observed_sets = observed_sets)

  expect_false(result$is_connected)
  expect_equal(result$num_components, 3)

  # Verify component assignments
  comp1 <- result$component_map[1]
  comp2 <- result$component_map[5]
  comp3 <- result$component_map[7]

  # Component 1: {1,2,3,4}
  expect_true(all(result$component_map[1:4] == comp1))

  # Component 2: {5,6}
  expect_true(all(result$component_map[5:6] == comp2))

  # Component 3: {7,8,9}
  expect_true(all(result$component_map[7:9] == comp3))

  # All components are distinct
  expect_true(length(unique(c(comp1, comp2, comp3))) == 3)
})

test_that("linear chain of variables", {
  # Test case: variables connected in a linear chain
  # 1-2, 2-3, 3-4, 4-5
  # Should still be fully connected

  observed_sets <- list(
    c(1, 2),
    c(2, 3),
    c(3, 4),
    c(4, 5)
  )

  result <- .check_graph_connectivity(p = 5, observed_sets = observed_sets)

  expect_true(result$is_connected)
  expect_equal(result$num_components, 1)
  expect_true(all(result$component_map == 1))
})

test_that("star topology (all connected through center)", {
  # Test case: variable 1 observed with all others, but others not with each other
  # Still forms a connected graph

  observed_sets <- list(
    c(1, 2),
    c(1, 3),
    c(1, 4),
    c(1, 5)
  )

  result <- .check_graph_connectivity(p = 5, observed_sets = observed_sets)

  expect_true(result$is_connected)
  expect_equal(result$num_components, 1)
  expect_true(all(result$component_map == 1))
})

test_that("duplicate observation sets don't affect connectivity", {
  # Same observation pattern repeated multiple times

  observed_sets <- list(
    c(1, 2, 3),
    c(1, 2, 3),
    c(1, 2, 3),
    c(3, 4),
    c(3, 4)
  )

  result <- .check_graph_connectivity(p = 4, observed_sets = observed_sets)

  expect_true(result$is_connected)
  expect_equal(result$num_components, 1)
})

test_that("large p with sparse observations", {
  # Test case: 100 variables, only a few observed
  # Variables 1-5 connected, rest unobserved

  observed_sets <- list(
    c(1, 2, 3),
    c(3, 4, 5)
  )

  result <- .check_graph_connectivity(p = 100, observed_sets = observed_sets)

  expect_false(result$is_connected)
  expect_equal(result$num_components, 96)  # {1,2,3,4,5} plus 95 singletons

  # Variables 1-5 should be in same component
  comp_main <- result$component_map[1]
  expect_true(all(result$component_map[1:5] == comp_main))

  # Variables 6-100 should each be in their own component
  expect_equal(length(unique(result$component_map[6:100])), 95)
})

# Integration tests with fit_covcomb()
test_that("fit_covcomb warns on disconnected graph", {
  # Create two disconnected components
  # Component 1: variables x1, x2
  # Component 2: variables x3, x4

  S1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  dimnames(S1) <- list(c("x1", "x2"), c("x1", "x2"))

  S2 <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
  dimnames(S2) <- list(c("x3", "x4"), c("x3", "x4"))

  S_list <- list(sample1 = S1, sample2 = S2)
  nu <- c(sample1 = 10, sample2 = 10)

  # Should produce a warning about disconnected components
  expect_warning(
    fit_covcomb(S_list, nu, scale_method = "none"),
    "The observation pattern creates 2 disconnected components"
  )
})

test_that("fit_covcomb does not warn on connected graph", {
  # Create a connected graph
  # Sample 1: x1, x2
  # Sample 2: x2, x3  (x2 connects the two samples)

  S1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  dimnames(S1) <- list(c("x1", "x2"), c("x1", "x2"))

  S2 <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
  dimnames(S2) <- list(c("x2", "x3"), c("x2", "x3"))

  S_list <- list(sample1 = S1, sample2 = S2)
  nu <- c(sample1 = 10, sample2 = 10)

  # Should NOT produce a connectivity warning
  result <- fit_covcomb(S_list, nu, scale_method = "none")
  expect_s3_class(result, "covcomb")
})

test_that("fit_covcomb connectivity check with three components", {
  # Create three disconnected components
  S1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  dimnames(S1) <- list(c("v1", "v2"), c("v1", "v2"))

  S2 <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
  dimnames(S2) <- list(c("v3", "v4"), c("v3", "v4"))

  S3 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  dimnames(S3) <- list(c("v5", "v6"), c("v5", "v6"))

  S_list <- list(s1 = S1, s2 = S2, s3 = S3)
  nu <- c(s1 = 20, s2 = 20, s3 = 20)

  expect_warning(
    fit_covcomb(S_list, nu),
    "3 disconnected components"
  )
})
