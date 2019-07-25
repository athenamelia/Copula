# estimation of nested A. copula parameters

loglikNC <- function(par, U, naclist) {
  # transform provided par (vector of real numbers) to a vector of parameter values
  # that are within the bounds for the copula we are estimating
  
  # need par[[i]] for each copula in the tree
  # if (length(naclist[[3]]) + length(naclist[[2]]) == 2) {
  #   par <- exp(par)
  # } else if (length(naclist[[3]]) == 3) {
  #   par <- -1 + exp(par)
  # }
  
  par <- exp(par)
  print(par)
  # evaluate log-likelihood at provided parameter values
  naclist <- updatePars(par, naclist, current_index = 1)
  updated_naclist <- naclist[[1]]
  result <- sum(dNC(U,updated_naclist, log = TRUE))
  print(result)
  return(result)
               
}

param_test <- c(.6, .3)
nlist_test <- list(th[1], 1:2, list(list(th[2], 3:4)))
U_test <- matrix(runif(n= 300), nrow = 100, ncol = 4)

optim_results <- optim(par = param_test,
                       fn = loglikNC,
                       U = U_test,
                       naclist = nlist_test,
                       method ="BFGS",
                       control = list(fnscale = -1))
optim_results


