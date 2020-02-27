#' Log likelihood of nested A. copula parameters
#'
#' @param nac_Node object of class
#' @param U pseudo-observations
#' @param par a vector of copula parameters on unbounded scale (real numbers)
#' @return Log-likelihood
#' @export
get_loglik <- function(nac_Node, U, par) {

  nac_Node <- set_par(nac_Node, par)[[1]]
  new_par <- transform_theta_unbounded_to_bounded(nac_Node)

  # if(is.na(par)) {
  #   browser()
  # }

  nac_Node <- set_par(nac_Node, new_par)[[1]]
  # evaluate log-likelihood at provided parameter values
  result <- sum(dncopula(nac_Node, U, log = TRUE))
  return(result)
}

