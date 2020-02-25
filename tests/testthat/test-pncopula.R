library(copula)
library(ncopula)
library(lattice)

test_that("pncopula works, no children", {
  tau <- c(0.05)
  th <- iTau(archmCopula("Frank"), tau = tau)

  nac_node_no_children <- new_nac_node("Frank", th[1], 1:2, list())
  U <- matrix(runif(n=200), nrow = 100, ncol = 2)

  cc_no_children <- frankCopula(th[1], dim = 2)
  p_test_no_children <- pCopula(U[,1:2], copula = cc_no_children)

  expect_equal(pncopula(nac_node_no_children, U), p_test_no_children)
})

test_that("pncopula works, with NULL", {
  tau <- c(0.05, 0.1, 0.2)
  th <- iTau(archmCopula("Frank"), tau = tau)

  nac_node_child1 <- new_nac_node("Frank", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Frank", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Frank", th[1], NULL, list(nac_node_child1, nac_node_child2))

  U <- matrix(runif(n=400), nrow = 100, ncol = 4)

  nested_cc2 <- frankCopula(th[3], dim = 2)
  nested_p_test2 <- pCopula(U[,3:4], copula = nested_cc2)
  nested_cc1 <- frankCopula(th[2], dim = 2)
  nested_p_test1 <- pCopula(U[,1:2], copula = nested_cc1)
  cc_null <- frankCopula(th[1], dim = 2)
  p_test_null <- pCopula(cbind(nested_p_test1, nested_p_test2), copula = cc_null)

  expect_equal(pncopula(nac_node_null, U), p_test_null)
})

test_that("pncopula works, normal nesting structure", {
  tau <- c(0.05, 0.1)
  th <- iTau(archmCopula("Clayton"), tau = tau)

  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))
  U <- matrix(runif(n=400), nrow = 100, ncol = 4)

  nested_cc <- claytonCopula(th[2], dim = 2)
  nested_p_test <- pCopula(U[,3:4], copula = nested_cc)
  cc <- claytonCopula(th[1], dim = 3)
  p_test <- pCopula(cbind(U[,1:2], nested_p_test), copula = cc)

  expect_equal(pncopula(nac_node_normal, U), p_test)
})

test_that("pncopula works, complicated nesting structure", {
  tau <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9)
  th <- iTau(archmCopula("Clayton"), tau = tau)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))
  U <- matrix(runif(n=1100), nrow = 100, ncol = 11)

  nested_complicated_cc3 <- claytonCopula(th[7], dim = 3)
  nested_complicated_p_test3 <- pCopula(U[,8:11], copula = nested_complicated_cc3)
  complicated_cc3 <- claytonCopula(th[6], dim = 2)
  p_test3 <- pCopula(cbind(U[,7], nested_complicated_p_test3), copula = complicated_cc3)

  nested_complicated_cc2 <- claytonCopula(th[5], dim = 2)
  nested_complicated_p_test2 <- pCopula(U[,5:6], copula = nested_complicated_cc2)
  complicated_cc2 <- claytonCopula(th[4], dim = 2)
  p_test2 <- pCopula(cbind(U[,4], nested_complicated_p_test2), copula = complicated_cc2)

  nested_complicated_cc1 <- claytonCopula(th[3], dim = 2)
  nested_complicated_p_test1 <- pCopula(U[,2:3], copula = nested_complicated_cc1)
  complicated_cc1 <- claytonCopula(th[2], dim = 2)
  p_test1 <- pCopula(cbind(U[,1], nested_complicated_p_test1), copula = complicated_cc1)

  parent_copula <- claytonCopula(th[1], dim = 3)
  p_test <- pCopula(cbind(p_test1, p_test2, p_test3), copula = parent_copula)

  expect_equal(pncopula(nac_node_full, U), p_test)
})
