library(ncopula)

test_that("is.nac_Node works, no children", {
  nac_node_no_children <- new_nac_node("Frank", 1, 1:2, list())
  U <- matrix(runif(n=200), nrow = 100, ncol = 2)

  expect_equal(is.nac_Node(nac_node_no_children), TRUE)
})

test_that("is.nac_Node works, with children", {
  child_node <- new_nac_node("Clayton", 2, 3:6, list())
  nac_node <- new_nac_node("Frank", 1, 1:2, list(child_node))
  U <- matrix(runif(n=200), nrow = 100, ncol = 2)

  expect_equal(is.nac_Node(nac_node), TRUE)
})

