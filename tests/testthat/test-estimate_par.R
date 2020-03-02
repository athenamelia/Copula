library(copula)
library(ncopula)
library(lattice)
test_that("estimate_par works, with subcopulas, independence copula", {
  set.seed(9999)
  U <- matrix(runif(n=1100), nrow = 100, ncol = 11)

  new_child <- new_nac_node("Clayton", 70, 8:11, list())
  new_node <- new_nac_node("Clayton", 60, 7, list(new_child))

  expect_equal(estimate_par(new_node, U), c(0.03723496, 0.07648524), tolerance = 0.001)
})

test_that("estimate_par works, with subcopulas, simulate data from Clayton copula", {
  set.seed(9999)
  cc <- claytonCopula(13, dim = 3)
  U <- rCopula(n = 300, copula = cc)

  new_child <- new_nac_node("Clayton", 70, 2:3, list())
  new_node <- new_nac_node("Clayton", 1000, 1, list(new_child))
  expect_equal(estimate_par(new_node, U), c(13.010, 13.647), tolerance = 0.001)
})


test_that("estimate_par works, with subcopulas, simulate data from Gumbel copula", {
  set.seed(9999)
  gc <- gumbelCopula(10, dim = 3)
  U <- rCopula(n = 300, copula = gc)

  new_child <- new_nac_node("Gumbel", 700, 2:3, list())
  new_node <- new_nac_node("Gumbel", 1000, 1, list(new_child))

  expect_equal(estimate_par(new_node, U), c(11.211, 11.071), tolerance = 0.001)
})

test_that("estimate_par works, with NULL and subcopulas", {
  tau <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9)
  th <- iTau(archmCopula("Clayton"), tau = tau)

  #set.seed(5432) # Error in copClayton@dacopula: theta is outside its defined interval
  set.seed(9999)
  cc <- claytonCopula(10, dim = 11)
  U <- rCopula(n = 1100, copula = cc)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))

  expect_equal(estimate_par(nac_node_full, U), , tolerance = 0.001)
})
