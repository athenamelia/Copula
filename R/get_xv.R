#' Cross Validation for nested A. copulas
#'
#' @param nac_Node an NAC node
#' @param x a data matrix that will be transformed to pseudo-observations
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

  if (is.null(k)) { k <- n}
  else {as.integer(k)}

  density <- 0
  folds <- createFolds(1:n, k, list = TRUE, returnTrain = FALSE)

  for (i in 1:k) {
    test_ind <- c(folds[[i]])
    train_set <- pobs(x[-test_ind, , drop = FALSE], ties.method = 'average')
    optim_results <- optim(par = param,
                           fn = get_loglik,
                           U = train_set,
                           nac_Node = nac_Node,
                           method ="BFGS",
                           control = list(fnscale = -1))

    full_set <- pobs(x, ties.method = 'average')
    test_set <- full_set[test_ind, ]

    new_nac_Node <- update_par(nac_Node, exp(optim_results$par), current_index = 1)[[1]]
    density <- density + sum(dncopula(new_nac_Node, test_set, log = TRUE))
  }
  density/k
}
