library(copula)
library(ncopula)

test_that("transform_theta_bounded_to_unbounded works, no children", {
  th <- 0.2488489
  nac_node_no_children <- new_nac_node("Clayton", th, 1:2, list())

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_no_children), 0.222222, tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, with NULL", {
  th <- c(4.301608, 5.236820, 6.883223)
  nac_node_child1 <- new_nac_node("Joe", th[2], 1:2, list())
  nac_node_child2 <- new_nac_node("Joe", th[3], 3:4, list())
  nac_node_null <- new_nac_node("Joe", th[1], NULL, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_null), c(1.194410, 1.443813, 1.772105), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with Gumbel", {
  th <- c(3.865181, 4.037732)
  nac_node_child <- new_nac_node("Gumbel", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Gumbel", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal), c(1.052632, 1.111111), tolerance = 0.001)
})

# NaNs produced
# test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with amh", {
#   th <- c(0.3775407, 0.4218990)
#   nac_node_child <- new_nac_node("Amh", th[2], 3:5, list())
#   nac_node_normal_amh <- new_nac_node("Amh", th[1], 1, list(nac_node_child))
#
#   expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal_amh), c(-0.5, 0.9), tolerance = 0.001)
# })

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with frank", {
  th <- c(2.477791, 1.860884, 2.917434)
  nac_node_child1 <- new_nac_node("Frank", th[2], 3:4, list())
  nac_node_child2 <- new_nac_node("Frank", th[3], 5:6, list())
  nac_node_normal <- new_nac_node("Frank", th[1], 1:2, list(nac_node_child1, nac_node_child2))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal), c(0.9073674, 1.8608840, 2.9174340), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, normal nesting structure with clayton", {
  th <- c(20.08554, 2979.95799)
  nac_node_child <- new_nac_node("Clayton", th[2], 3:4, list())
  nac_node_normal <- new_nac_node("Clayton", th[1], 1:2, list(nac_node_child))

  expect_equal(transform_theta_bounded_to_unbounded(nac_node_normal), c(3, 8), tolerance = 0.001)
})

test_that("transform_theta_bounded_to_unbounded works, NAC with complex nesting structure", {
  th <- c(0.105, 0.222, 0.5, 1.333, 3, 8, 18)

  nac_node_child11 <- new_nac_node("Clayton", th[3], 2:3, list())
  nac_node_child1 <- new_nac_node("Clayton", th[2], 1, list(nac_node_child11))

  nac_node_child21 <- new_nac_node("Clayton", th[5], 5:6, list())
  nac_node_child2 <- new_nac_node("Clayton", th[4], 4, list(nac_node_child21))

  nac_node_child31 <- new_nac_node("Clayton", th[7], 8:11, list())
  nac_node_child3 <- new_nac_node("Clayton", th[6], 7, list(nac_node_child31))

  nac_node_full <- new_nac_node("Clayton", th[1], NULL, list(nac_node_child1, nac_node_child2, nac_node_child3))
  expect_equal(transform_theta_bounded_to_unbounded(nac_node_full), c(-2.251, 0.201, 0.405, 0.847, 1.386, 2.197, 2.890), tolerance = 0.001)
})
