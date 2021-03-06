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
  print(paste0("in get_loglik, bounded param are ", new_par))

  nac_Node <- set_par(nac_Node, new_par)[[1]]
  print(paste0("after set_par, par = ", get_par(nac_Node)))
  result <- sum(dncopula(nac_Node, U, log = TRUE))
  print(paste0("result = ", result))

  if (result == -Inf || is.nan(result)) {
    result <- -.Machine$double.xmax
  }
  print(paste0("result = ", result))
  return(result)
}

