#'estimation of nested A. copula parameters
#'
#' @param nac_Node an NAC node 
#' @param U pseudo-observations
#' @param par a vector of parameters 
#' @return result obtained from optimizing provided parameters
get_loglik <- function(nac_Node, U, par) {
  # transform provided par (vector of real numbers) to a vector of parameter values
  # that are within the bounds for the copula we are estimating
  orig_par <- par
  par <- set_theta_from_unbounded(nac_Node, par, current_index = 1)[[1]]  
  
  # if(is.na(par)) {
  #   browser()
  # }
  
  # evaluate log-likelihood at provided parameter values
  new_nac_Node <- update_par(nac_Node, par, current_index = 1)[[1]]
  result <- sum(get_density(new_nac_Node, U, log = TRUE))
  return(result)
}
