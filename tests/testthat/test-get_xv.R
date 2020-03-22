library(copula)
library(ncopula)
library(caret)

test_that("get_xv works, no children", {
  tau <- c(100)
  th <- iTau(archmCopula("Clayton"), tau = tau)
  set.seed(9999)
  X <- matrix(rnorm(n = 200, mean = 0, sd = 1), ncol = 2, nrow = 100)
  nac_node_no_children <- new_nac_node("Clayton", th, 1:2, list())

  expect_equal(get_xv(X, nac_node_no_children, k = 10), -0.0147, tolerance = 0.001)
})

test_that("get_xv works, with NULL", {
  tau <- c(0.1, 0.2, 0.3)
  th <- iTau(archmCopula("Joe"), tau = tau)
  set.seed(1234)
  X <- matrix(rnorm(n = 400, mean = 0, sd = 1), ncol = 4, nrow = 100)
  nac_node_child1 <- new_nac_node("Joe", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Joe", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Joe", th[1], NULL, list(nac_node_child1, nac_node_child2))

  expect_equal(get_xv(X, nac_node_null, k = 10), 2.22189, tolerance = 0.001)
})

test_that("get_xv works, normal nesting structure with Gumbel", {
  tau <- c(2, 10)
  th <- iTau(archmCopula("Gumbel"), tau = tau)
  set.seed(1234)
  X <- matrix(rnorm(n = 400, mean = 0, sd = 1), ncol = 4, nrow = 100)
  nac_node_child <- new_nac_node("Gumbel", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Gumbel", th[1], 1:2, list(nac_node_child))

  expect_equal(get_xv(X, nac_node_normal, k = 10), -0.0702, tolerance = 0.001)
})

test_that("get_xv works, normal nesting structure with amh", {
  th <- c(-0.5, 0.9)
  set.seed(9999)
  X <- matrix(rnorm(n = 300, mean = 0, sd = 1), ncol = 3, nrow = 100)
  nac_node_child <- new_nac_node("Amh", th[2], 2:3, list())
  nac_node_normal_amh <- new_nac_node("Amh", th[1], 1, list(nac_node_child))

  expect_equal(get_xv(X, nac_node_normal_amh, k = 5), 0.2940288, tolerance = 0.001)
})

test_that("get_xv works, normal nesting structure with frank", {
  tau <- c(0.1, 0.2, 0.3)
  th <- iTau(archmCopula("Frank"), tau = tau)

  set.seed(9999)
  X <- matrix(rnorm(n = 600, mean = 0, sd = 1), ncol = 6, nrow = 100)

  nac_node_child1 <- new_nac_node("Frank", th[2], 3:4, list())
  nac_node_child2 <- new_nac_node("Frank", th[3], 5:6, list())
  nac_node_normal <- new_nac_node("Frank", th[1], 1:2, list(nac_node_child1, nac_node_child2))

  expect_equal(get_xv(X, nac_node_normal, k = 7), 1.2, tolerance = 0.001)
})

test_that("get_xv works, normal nesting structure with clayton", {
  tau <- c(0.6, 0.8)
  th <- iTau(archmCopula("Clayton"), tau = tau)

  set.seed(9999)
  X <- matrix(rnorm(n = 400, mean = 0, sd = 1), ncol = 4, nrow = 100)

  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  expect_equal(get_xv(X, nac_node_normal, k = 10), 0.085, tolerance = 0.001)
})

test_that("get_xv works, NAC with complex nesting structure", {
  th <- c(-2.251, 0.201, 0.405, 0.847, 1.386, 2.197, 2.890)

  set.seed(9999)
  X <- matrix(rnorm(n = 1100, mean = 0, sd = 1), ncol = 11, nrow = 100)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1,
                                                             nac_node_child2,
                                                             nac_node_child3))
  expect_equal(get_xv(X, nac_node_full, k = 15), -Inf, tolerance = 0.001)
})
