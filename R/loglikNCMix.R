#' @param parlist vector of parameters of copulas in the mixture to be transformed into weights
#' @param mixlist list of nesting structures of all mixture copulas
#' @param U pseudo observations
#' @return estimation of mixture copula parameters

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

