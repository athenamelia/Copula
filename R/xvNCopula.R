# Cross Validation for nested A. copulas
library(caret)

set.seed(1234)
xvNCopula <- function(x, naclist, k = NULL) {
  if (!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  stopifnot(is.numeric(x), (d <- ncol(x)) > 1, (n <- nrow(x)) > 0)
  
  if (is.null(k)) { k <- n}
  else {as.integer(k)}
  
  density <- 0
  folds <- createFolds(1:n, k, list = TRUE, returnTrain = FALSE)

  for (i in 1:k) {
    test_ind <- c(folds[[i]])
    train_set <- pobs(x[-test_ind, , drop = FALSE], ties.method = 'average')
    optim_results <- optim(par = param,
                           fn = loglikNC,
                           U = train_set,
                           naclist = naclist,
                           method ="BFGS",
                           control = list(fnscale = -1))
    
    full_set <- pobs(x, ties.method = 'average')
    test_set <- full_set[test_ind, ]
    
    update_naclist <- updatePars(exp(optim_results$par), naclist, current_index = 1)
    density <- density + sum(dNC(test_set, update_naclist[[1]], log = TRUE))
  }
  density
}


