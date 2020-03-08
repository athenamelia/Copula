library(copula)
library(ncopula)

test_that("get_cor_theta works, Clayton with NULL and subcopulas", {
  set.seed(1234)
  tau <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9)
  th <- iTau(archmCopula("Clayton"), tau = tau)

  U <- matrix(runif(n = 1100), nrow = 100, ncol = 11)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))

  expect_equal(get_cor_theta(nac_node_full, U), 0.0758, tolerance = 0.001)
})

test_that("get_cor_theta works, Gumbel with NULL and subcopulas", {
  set.seed(4321)
  tau <- c(0.4, 0.6, 0.8, 0.9, 0.95)
  th <- iTau(archmCopula("Gumbel"), tau = tau)

  U <- matrix(runif(n = 600), nrow = 100, ncol = 6)

  nac_node_child11 <- new_nac_node("Gumbel", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Gumbel", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Gumbel", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Gumbel", th[4], 4, list(nac_node_child21))

  nac_node_full <- new_nac_node("Gumbel", th[1], NULL, list(nac_node_child1, nac_node_child2))
  expect_equal(get_cor_theta(nac_node_full, U), 1.0356, tolerance = 0.001)
})

test_that("get_cor_theta works, Frank with no subcopulas", {
  set.seed(9999)
  tau <- c(0.7)
  th <- iTau(archmCopula("Frank"), tau = tau)

  U <- matrix(runif(n = 300), nrow = 100, ncol = 3)
  nac_node <- new_nac_node("Frank", th[1], 1:3, list())

  expect_equal(get_cor_theta(nac_node, U), 0.05, tolerance = 0.001)
})
