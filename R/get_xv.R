#' Cross validation for nested Archimedean copulas
#'
#' @param nac_Node an NAC node
#' @param x a data matrix to be transformed to pseudo-observations
#' @param k the number of data blocks; if k = NULL, nrow(x) blocks are considered
#' @return A real number equal to the cross-validation criterion
#' @export
get_xv <- function(x, nac_Node, k = NULL) {
  if (!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }

  d <- ncol(x)
  n <- nrow(x)
  stopifnot(is.numeric(x), d > 1, n > 0)

  if (is.null(k)) {k <- n}
  else {as.integer(k)}

  density <- 0
  folds <- caret::createFolds(1:n, k, list = TRUE, returnTrain = FALSE)

  for (i in 1:k) {
    test_ind <- folds[[i]]
    train_set <- pobs(x[-test_ind, , drop = FALSE], ties.method = 'average')
    new_nac_Node <- estimate_par(nac_Node, train_set)

    full_set <- pobs(x, ties.method = 'average')
    test_set <- full_set[test_ind, ]
    density <- density + sum(dncopula(new_nac_Node, test_set, log = TRUE))
  }
  return(density/k)
}
