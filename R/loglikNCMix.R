# estimation of mixture copula parameters

loglikNCMix <- function(parlist, U, mixlist) {
  num_of_copulas <- length(mixlist)
  sum_params <- 0
  W <- vector()
  
  # softmax transformation 
  # transform par into weights
  for (i in 1:num_of_copulas) {
    sum_params <- sum_params + exp(parlist[[i]])
  }
  
  for (i in 1:num_of_copulas) {
    W[[i]] <- exp(parlist[[i]]) / sum_params
  }
  
  # evaluate log-likelihood 
  tryCatch(sum(dNCMix(U, W, mixlist, log = TRUE)), 
           error = function(e) {
             return(-1 * .Machine$double.xmax)
           }) 
  
}

