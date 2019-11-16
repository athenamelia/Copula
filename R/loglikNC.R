#'estimation of nested A. copula parameters
#'
#' @param par vector of copula parameters
#' @param naclist list specifying nesting structure
#' @param U observations
#' @return result obtained from optimizing provided parameters
#' @export
loglikNC <- function(par, U, naclist) {
  # transform provided par (vector of real numbers) to a vector of parameter values
  # that are within the bounds for the copula we are estimating
  orig_par <- par
  par <- transform_par(par, naclist)[[1]]
  # if(is.na(par)) {
  #   browser()
  # }

  # evaluate log-likelihood at provided parameter values
  updated_naclist <- updatePars(par, naclist, current_index = 1)[[1]]
  result <- sum(dNC(U, updated_naclist, log = TRUE))
  return(result)
}
