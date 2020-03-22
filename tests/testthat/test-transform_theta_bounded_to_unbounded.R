library(copula)
library(ncopula)

test_that("transform_theta_bounded_to_unbounded works, no children", {
  th <- -0.1896
  nac_node_no_children <- new_nac_node("Clayton", th, 1:2, list())

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_no_children), 0.222189, tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, with NULL", {
  th <- c(2.459, 2.656, 2.929)
  nac_node_child1 <- new_nac_node("Joe", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Joe", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Joe", th[1], NULL, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_null), c(1.194424, 1.444166, 1.772003), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with Gumbel", {
  th <- c(2.352009, 2.395683)
  nac_node_child <- new_nac_node("Gumbel", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Gumbel", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal), c(1.052632, 1.111111), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with amh", {
  th <- c(0.3775407, 0.4218990)
  nac_node_child <- new_nac_node("Amh", th[2], 3:5, list())
  nac_node_normal_amh <- new_nac_node("Amh", th[1], 1, list(nac_node_child))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal_amh), c(0.7943768, -0.3149826), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with frank", {
  th <- c(1.246397, 1.860884, 2.917434)
  nac_node_child1 <- new_nac_node("Frank", th[2], 3:4, list())
  nac_node_child2 <- new_nac_node("Frank", th[3], 5:6, list())
  nac_node_normal <- new_nac_node("Frank", th[1], 1:2, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal), c(0.9073676, 1.8608838, 2.9174344), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with clayton", {
  th <- c(3.048587, 7.000335)
  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal), c(3, 8), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, NAC with complex nesting structure", {
  th <- c(0.1, -0.201, -0.084, 0.204, 0.609, 1.302, 2.944)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))
  expect_equal(transform_theta_bounded_to_unbounded(nac_node_full), c(-2.251, 0.201, 0.405, 0.847, 1.386, 2.197, 2.890), tolerance = 0.001)
})
