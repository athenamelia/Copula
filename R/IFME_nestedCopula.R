# estimation of nested A. copula parameters

loglikNC <- function(par, U, family, naclist) {
  # transform provided par (vector of real numbers) to a vector of parameter values
  # that are within the bounds for the copula we are estimating
  par <- exp(par)
  
  # evaluate log-likelihood at provided parameter values
  naclist <- updatePars(par, naclist, current_index = 1)
  updated_naclist <- naclist[[1]]
  sum(dNC(U,updated_naclist, log = TRUE))
}

param_test <- c(.6, .3)
nlist_test <- list(th[1], 1:2, list(list(th[2], 3:4)))
U_test <- matrix(runif(n= 300), nrow = 100, ncol = 4)

optim_results <- optim(par = param_test,
  fn = loglikNC,
  U = U_test,
  family = "Clayton",
  naclist = nlist_test,
  method ="BFGS",
  control = list(fnscale = -1))
optim_results


