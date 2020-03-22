#' estimate copula parameters inside out
#'
#' @param U pseudo-observations
#' @param nac_Node an nac node
#' @return a list of nac node, estimated param, and result of optim function
#' @export
estimate_par <- function(nac_Node, U) {

  if (has_subcopula(nac_Node) == FALSE) {
    new_theta <- get_cor_theta(nac_Node, U)
    nac_Node <- set_theta(nac_Node, new_theta)
    nac_Node <- estimate_nac(nac_Node, U)
  }

  else if (has_subcopula(nac_Node)) {
    num_subcopulas <- count_subcopula(nac_Node)

    for (i in seq_len(num_subcopulas)) {
      child_copula <- get_subcopula(nac_Node, i)
      new_child_copula <- estimate_par(child_copula, U)
      nac_Node <- set_subcopula(nac_Node, i, new_child_copula)
    }

    new_theta <- get_cor_theta(nac_Node, U)
    nac_Node <- set_theta(nac_Node, new_theta)
    nac_Node <- estimate_nac(nac_Node, U)
  }
  return(nac_Node)
}

#' A helper method to estimate the parameters of one Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @param U pseudo-observations
#' @return nac_Node a nested Archimedean copula with new estimates
#' @export
estimate_nac <- function(nac_Node, U) {
  par <- transform_theta_bounded_to_unbounded(nac_Node)
  print(paste0("in estimate function, unbounded estimate before calling optim is ", par))

  result <- optim(par = par,
                  fn = get_loglik,
                  U = U,
                  nac_Node = nac_Node,
                  method ="BFGS",
                  control = list(fnscale = -1))

  estimate <- result$par
  print(paste0("in estimate function, unbounded estimate is ", estimate))
  bounded_estimate <- transform_theta_unbounded_to_bounded(set_par(nac_Node, estimate)[[1]])
  print(paste0("in estimate function, bounded estimate is ", estimate))
  nac_Node <- set_par(nac_Node, bounded_estimate)[[1]]

  return(nac_Node)
}
