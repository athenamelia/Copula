library(copula)
library(ncopula)

test_that("transform_theta_unbounded_to_bounded works, no children", {
  tau <- c(0.1)
  th <- iTau(archmCopula("Clayton"), tau = tau) # 0.2222222
  nac_node_no_children <- new_nac_node("Clayton", th[1], 1:2, list())

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_no_children), -0.1896, tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, with NULL", {
  tau <- c(0.1, 0.2, 0.3)
  th <- iTau(archmCopula("Joe"), tau = tau) # 1.194410 1.443813 1.772105
  nac_node_child1 <- new_nac_node("Joe", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Joe", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Joe", th[1], NULL, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_null), c(2.459, 2.656, 2.929), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, normal nesting structure with Gumbel", {
  tau <- c(0.05, 0.1)
  th <- iTau(archmCopula("Gumbel"), tau = tau) # 1.052632 1.111111
  nac_node_child <- new_nac_node("Gumbel", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Gumbel", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal), c(2.352, 2.396), tolerance = 0.001)
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

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal), c(1.246, 1.861, 2.917), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, normal nesting structure with clayton", {
  tau <- c(0.6, 0.8)
  th <- iTau(archmCopula("Clayton"), tau = tau) # 3 8

  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_normal), c(3.049, 7), tolerance = 0.001)
})

test_that("transform_theta_unbounded_to_bounded works, NAC with complex nesting structure", {
  tau <- c(-2.251, 0.201, 0.405, 0.847, 1.386, 2.197, 2.890)
  th <- iTau(archmCopula("Clayton"), tau = tau)
  # -1.3848047  0.5031289  1.3613445 11.0718954 -7.1813472 -3.6708438 -3.0582011

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))

  expect_equal(transform_theta_unbounded_to_bounded(nac_node_full), c(0.22, -0.024, 0.589, 10.072, -0.999, -0.975, 0.046), tolerance = 0.001)
})
