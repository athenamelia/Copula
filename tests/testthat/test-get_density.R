library(copula)
library(lattice)

test_that("get_density works, no children", {
  tau <- c(0.05, 0.1) # Kendall's tau
  th <- iTau(archmCopula("Clayton"), tau = tau) # corresponding parameters

  nlist_no_children <- list("Clayton", th[1], 1:2)
  U <- matrix(runif(n=200), nrow = 100, ncol = 2)

  density_no_children <- dCopula(U[,1:2], claytonCopula(th[1], dim = 2), log = TRUE)
  expect_equal(get_density(nlist_no_children, U), density_no_children)
})

test_that("get_density works, with NULL", {
  tau <- c(0.05, 0.1, 0.2) # Kendall's tau
  th <- iTau(archmCopula("Frank"), tau = tau) # corresponding parameters

  nlist_null <- list("Frank", th[1], NULL, list(list("Frank", th[2], 1:2),
                                                list("Frank", th[3], 3:4)))
  U <- matrix(runif(n=200), nrow = 100, ncol = 4)

  density_child1 <- dCopula(U[,1:2], frankCopula(th[2], dim = 2), log = TRUE)
  density_child2 <- dCopula(U[,3:4], frankCopula(th[3], dim = 2), log = TRUE)

  cc1 <- frankCopula(th[2], dim = 2)
  cdf_child1 <- pCopula(U[,1:2], copula = cc1)
  cc2 <- frankCopula(th[3], dim = 2)
  cdf_child2 <- pCopula(U[,3:4], copula = cc2)

  density_parent <- dCopula(cbind(cdf_child1, cdf_child2), frankCopula(th[1], dim = 2), log = TRUE)
  density_null <- density_parent + density_child1 + density_child2

  expect_equal(get_density(nlist_null, U), density_null)
})

test_that("get_density works, normal listing structure", {
  tau <- c(0.05, 0.1) # Kendall's tau
  th <- iTau(archmCopula("Clayton"), tau = tau) # corresponding parameters

  nlist <- list("Clayton", th[1], 1:2, list(list("Clayton", th[2], 3:4)))
  U <- matrix(runif(n=200), nrow = 100, ncol = 4)

  density_child <- dCopula(U[,3:4], claytonCopula(th[2], dim = 2), log = TRUE)
  cc <- claytonCopula(th[2], dim = 2)
  cdf_child <- pCopula(U[,3:4], copula = cc)

  density_parent <- dCopula(cbind(U[,1:2], cdf_child), claytonCopula(th[1], dim = 3), log = TRUE)
  density_test <- density_parent + density_child

  expect_equal(get_density(nlist, U), density_test)
})


test_that("get_density works, complicated nesting structure", {
  tau <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9) # Kendall's tau
  th <- iTau(archmCopula("Clayton"), tau = tau) # corresponding parameters

  nlist_complicated <- list("Clayton", th[1], NULL, list(list("Clayton", th[2], 1, list(list("Clayton", th[3], 2:3))), # NAC structure
                                              list("Clayton", th[4], 4, list(list("Clayton", th[5], 5:6))),
                                              list("Clayton", th[6], 7, list(list("Clayton", th[7], 8:11)))))
  U <- matrix(runif(n=220), nrow = 100, ncol = 11)

  # density of nested children copulas
  d_test_nested_cplex1 <- dCopula(U[,2:3], claytonCopula(th[3], dim = 2), log = TRUE)
  d_test_nested_cplex2 <- dCopula(U[,5:6], claytonCopula(th[5], dim = 2), log = TRUE)
  d_test_nested_cplex3 <- dCopula(U[,8:11], claytonCopula(th[7], dim = 4), log = TRUE)

  # cdf of nested children copulas
  cc_nested_cplex1 <- claytonCopula(th[3], dim = 2)
  cdf_nested_child1 <- pCopula(U[,2:3], copula = cc_nested_cplex1)
  cc_nested_cplex2 <- claytonCopula(th[5], dim = 2)
  cdf_nested_child2 <- pCopula(U[,5:6], copula = cc_nested_cplex2)
  cc_nested_cplex3 <- claytonCopula(th[7], dim = 4)
  cdf_nested_child3 <- pCopula(U[,8:11], copula = cc_nested_cplex3)

  # density of children copulas
  d_test_cplex1 <- dCopula(cbind(U[,1], cdf_nested_child1), claytonCopula(th[2], dim = 2), log = TRUE)
  d_test_cplex2 <- dCopula(cbind(U[,4], cdf_nested_child2), claytonCopula(th[4], dim = 2), log = TRUE)
  d_test_cplex3 <- dCopula(cbind(U[,7], cdf_nested_child3), claytonCopula(th[6], dim = 2), log = TRUE)

  # cdf of children copulas
  cc_cplex1 <- claytonCopula(th[2], dim = 2)
  cdf_child1 <- pCopula(cbind(U[,1], cdf_nested_child1), copula = cc_cplex1)
  cc_cplex2 <- claytonCopula(th[4], dim = 2)
  cdf_child2 <- pCopula(cbind(U[,4], cdf_nested_child2), copula = cc_cplex2)
  cc_cplex3 <- claytonCopula(th[6], dim = 2)
  cdf_child3 <- pCopula(cbind(U[,7], cdf_nested_child3), copula = cc_cplex3)

  # density of big tree
  d_test_parent <- dCopula(cbind(cdf_child1, cdf_child2, cdf_child3), claytonCopula(th[1], dim = 3), log = TRUE)
  d_test <-  (d_test_nested_cplex1 + d_test_cplex1 + (d_test_nested_cplex2 + d_test_cplex2) + (d_test_nested_cplex3 + d_test_cplex3) + d_test_parent)

  # check
  expect_equal(get_density(nlist_complicated, U), d_test)
})


