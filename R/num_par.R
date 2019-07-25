# Cross Validation for nested A. copulas
library(caret)
set.seed(1234)

xvNCopula1 <- function(x, naclist, k = NULL) {
  if (!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  stopifnot(is.numeric(x), (d <- ncol(x)) > 1, (n <- nrow(x)) > 0)
  
  x <- x[sample(nrow(x)),]
  n <- nrow(x)
  
  if (is.null(k)) { k <- n}
  else {as.integer(k)}
  
  density <- 0
  
  p <- n%/%k
  m <- rep(p, k)
  r <- n - k * p
  if (r > 0) 
    m[seq_len(r)] <- p + 1
  
  b <- c(0, cumsum(m))
  
  for (i in seq_len(k)) {
    sel <- (b[i] + 1):b[i + 1]
    u <- pobs(x[-sel, , drop = FALSE], ties.method = 'average')
    optim_results <- optim(par = param,
                           fn = loglikNC,
                           U = u,
                           family = "Clayton",
                           naclist = naclist,
                           method ="BFGS",
                           control = list(fnscale = -1))
    
    u_full <- pobs(x, ties.method = 'average')
    v <- u_full[sel, ]
    
    update_naclist <- updatePars(exp(optim_results$par), naclist, current_index = 1)
    
    density <- density + dNC(v, update_naclist[[1]], log = TRUE)
  }
  density/n
}

test <- function(n, k) {
  p <- n%/%k
  m <- rep(p, k)
  r <- n - k * p
  if (r > 0) 
    m[seq_len(r)] <- p + 1
  
  b <- c(0, cumsum(m))
  for (i in seq_len(k)) {
    sel <- (b[i] + 1):b[i + 1]
    print(sel)
  }
}

# folds not equal sizes
test1 <- function(n, k) {
  fold <- createFolds(rnorm(n), k, list = TRUE, returnTrain = FALSE)
  for (i in 1:k) {
    test_ind <- c(fold[[i]])
    print(length(test_ind))
  }
}

# number of parameters
num_param <- function(naclist) {
  num_param <- 0
  if (length(naclist) == 3) {
    num_param <- num_param + 1
    sub_list <- naclist[[3]]
    
    for (i in 1:length(sub_list)) {
      num_param <- num_param + num_param(sub_list[[i]])
    }
  }
  
  else if (length(naclist) == 2) {
    num_param <- num_param + 1
  }
  
  return(num_param)
}


