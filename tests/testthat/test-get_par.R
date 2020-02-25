library(copula)
library(ncopula)
library(lattice)

test_that("get_par works", {
  tau <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9) # Kendall's tau
  th <- iTau(archmCopula("Clayton"), tau = tau) # corresponding parameters

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))

  U <- matrix(runif(n=1100), nrow = 100, ncol = 11)

  expect_equal(get_par(nac_node_full), c(0.105, 0.222, 0.5, 1.333, 3, 8, 18), tolerance = 0.001)
})
