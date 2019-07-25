# loop optim function many many times 

fitNC <- function(num_of_times, U, naclist) {
  n <- num_param(naclist)
  loglik_value <- -Inf
  results <- 0
  
  for (i in 1: num_of_times) {
    param <- iTau(archmCopula('Clayton'), tau = runif(n = n, min = 0, max = 1))
    print(param)
    optim_results <- optim(par = log(param+1),
                           fn = loglikNC,
                           U = U,
                           naclist = naclist,
                           method ="BFGS",
                           control = list(fnscale = -1))
    
    if (loglik_value < optim_results$value) {
      loglik_value <- optim_results$value
      results <- optim_results
    }
  }
  return(results)
}

