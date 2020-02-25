library(copula)
library(ncopula)

test_that("transform_theta_unbounded_to_bounded works, no children", {
  tau <- c(0.1)
  th <- iTau(archmCopula("Clayton"), tau = tau) # 0.2222222
  nac_node_no_children <- new_nac_node("Clayton", th[1], 1:2, list())

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_no_children), 0.249, tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, with NULL", {
  tau <- c(0.1, 0.2, 0.3)
  th <- iTau(archmCopula("Joe"), tau = tau) # 1.194410 1.443813 1.772105
  nac_node_child1 <- new_nac_node("Joe", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Joe", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Joe", th[1], NULL, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_null), c(4.302, 5.237, 6.883), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, normal nesting structure with Gumbel", {
  tau <- c(0.05, 0.1)
  th <- iTau(archmCopula("Gumbel"), tau = tau) # 1.052632 1.111111
  nac_node_child <- new_nac_node("Gumbel", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Gumbel", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal), c(3.865, 4.038), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, normal nesting structure with amh", {
  th <- c(-0.5, 0.9)
  nac_node_child <- new_nac_node("Amh", th[2], 3:4, list())
  nac_node_normal_amh <- new_nac_node("Amh", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal_amh), c(0.378, 0.422), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, normal nesting structure with frank", {
  tau <- c(0.1, 0.2, 0.3)
  th <- iTau(archmCopula("Frank"), tau = tau) # 0.9073676 1.8608838 2.9174344

  nac_node_child1 <- new_nac_node("Frank", th[2], 3:4, list())
  nac_node_child2 <- new_nac_node("Frank", th[3], 5:6, list())
  nac_node_normal <- new_nac_node("Frank", th[1], 1:2, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal), c(2.478, 1.861, 2.917), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, normal nesting structure with clayton", {
  tau <- c(0.6, 0.8)
  th <- iTau(archmCopula("Clayton"), tau = tau) # 3 8

  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal), c(20.086, 2979.958), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, NAC with complex nesting structure", {
  th <- c(-2.251, 0.201, 0.405, 0.847, 1.386, 2.197, 2.890)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))
  expect_equal(transform_theta_unbounded_to_bounded(nac_node_full), c(0.105, 0.222, 0.5, 1.333, 3, 8, 18), tolerance = 0.001)
})
