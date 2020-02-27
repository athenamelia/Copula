library(copula)
library(lattice)

test_that("get_loglik works, no children", {
  set.seed(1234)
  tau <- c(0.05) # Kendall's tau
  th <- iTau(archmCopula("Gumbel"), tau = tau) # corresponding parameters

  nac_node_no_children <- new_nac_node("Gumbel", th[1], 1:2, list())
  U <- matrix(runif(n=200), nrow = 100, ncol = 2)

  expect_equal(get_loglik(nac_node_no_children, U, th), -267.128, tolerance = 0.001)
})

test_that("get_loglik works, with NULL", {
  set.seed(1234)
  tau <- c(0.05, 0.1, 0.2) # Kendall's tau
  th <- iTau(archmCopula("Frank"), tau = tau) # corresponding parameters

  nac_node_child1 <- new_nac_node("Frank", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Frank", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Frank", th[1], NULL, list(nac_node_child1, nac_node_child2))
  U <- matrix(runif(n=400), nrow = 100, ncol = 4)

  expect_equal(get_loglik(nac_node_null, U, th), -3.892, tolerance = 0.001)
})

test_that("get_loglik works, normal nesting structure", {
  set.seed(1234)
  tau <- c(0.05, 0.1) # Kendall's tau
  th <- iTau(archmCopula("Clayton"), tau = tau) # corresponding parameters

  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  U <- matrix(runif(n=400), nrow = 100, ncol = 4)
  expect_equal(get_loglik(nac_node_normal, U, th), -77.610, tolerance = 0.001)
})

test_that("get_loglik works, normal nesting structure with new set of parameters", {
  set.seed(1234)
  tau <- c(0.05, 0.1) # Kendall's tau
  th <- iTau(archmCopula("Clayton"), tau = tau) # corresponding parameters

  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  U <- matrix(runif(n=400), nrow = 100, ncol = 4)
  new_par <- c(0.123456, 2.22446688)

  expect_equal(get_loglik(nac_node_normal, U, new_par), -665.605, tolerance = 0.001)
})


