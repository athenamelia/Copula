#' Log likelihood of nested A. copula parameters
#'
#' @param nac_Node object of class
#' @param U pseudo-observations
#' @return Log-likelihood
#' @export
get_loglik <- function(nac_Node, U, par) {
  # transform provided par (vector of real numbers) to bounded
  new_theta <- transform_theta_unbounded_to_bounded(nac_Node)

  # if(is.na(par)) {
  #   browser()
  # }

  # evaluate log-likelihood at provided parameter values
  new_nac_Node <- set_par(nac_Node, new_theta)[[1]]
  result <- sum(dncopula(new_nac_Node, U, log = TRUE))
  return(result)
}
