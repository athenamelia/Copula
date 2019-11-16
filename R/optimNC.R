
#' Repeat optim function many many times - for a nested copula that belongs to 1 specific family
#'
#' @param num_of_times number of times to run optim function
#' @param U matix of pseudo observation
#' @param naclist list specifying nesting structure
#' @param num_of_times number of times to call optim function on NC
#' @param naclist list specifying nesting structure
#' @param U pseudo observations
#' @return A largest result obtained from running optim function many times

optimNC <- function(num_of_times, U, naclist) {
  n <- num_par(naclist)
  loglik_value <- -Inf
  results <- 0
  family <- naclist[[1]]

  for (i in 1: num_of_times) {

      if (family == 'Clayton') {
          param <- iTau(archmCopula('Clayton'), tau = runif(n = n, min = 0, max = 1))
      } else if (family == 'Frank') {
          param <- iTau(archmCopula('Frank'), tau = runif(n = n, min = 0, max = 1))
      } else if (family == 'Gumbel') {
          param <- iTau(archmCopula('Gumbel'), tau = runif(n = n, min = 0, max = 1))
      } else if (family == 'Joe') {
          param <- iTau(archmCopula('Joe'), tau = runif(n = n, min = 0, max = 1))
      } else if (family == 'Amh') {
          param <- iTau(archmCopula('Amh'), tau = runif(n = n, min = 0, max = 1))
      }

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

