#' estimate copula parameters inside out
#'
#' @param U pseudo-observations
#' @param nac_Node an nac node
#' @return a list of nac node, estimated param, and result of optim function
#' @export
estimate_param <- function(nac_Node, U) {
  param_list <- c()

  if (has_subcopulas(nac_Node)) {
    subcopulas_list <- get_subcopula(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

    # get initial parameter estimates for each descendant copula
    for (i in 1:length(subcopulas_list)) {
      child_copula <- subcopulas_list[[i]]
      temp <- estimate_param(child_copula, U)
      nac_Node <- set_subcopula(nac_Node, i, temp[[1]]) # update the naclist of the parent copula with new child-list
      child_copula <- temp[[1]]
      param_list <- c(param_list, temp[[2]]) # save new param to be used for optim fcn later
      V[,i] <- get_cdf(child_copula, U)  # estimate cdf using child copula with new param
    }

    # get initial parameter estimates for the current copula
    # get theta based on U and child copula's cdf,
    # then call optim to update current copula and all child copulas
    list_pNC <- nac_Node[-4]
    UV_concat <- cbind(U[, list_pNC[[3]]], V)
    list_pNC[[3]] <- seq_len(ncol(UV_concat))

    # get theta for list
    theta <- get_iniCor(list_pNC, UV_concat)
    nac_Node[[2]] <- theta
    param_list <- c(
      set_theta_unbounded(nac_Node, get_iniTheta(nac_Node))[[1]],
      param_list
    )

    result <- optim(par = param_list,
                    fn = get_loglik,
                    U = U,
                    nac_Node = nac_Node,
                    method ="BFGS",
                    control = list(fnscale = -1))
    nac_Node <- update_par(nac_Node, set_theta_from_unbounded(nac_Node, result$par)[[1]])[[1]]
    param_list <- result$par

  } else if (has_subcopulas(nac_Node) == FALSE) {
    # check theta
    theta <- get_iniCor(nac_Node, U)
    nac_Node[[2]] <- theta

    result <- optim(par = set_theta_unbounded(nac_Node, get_iniTheta(nac_Node))[[1]],
                    fn = get_loglik,
                    U = U,
                    nac_Node = nac_Node,
                    method ="BFGS",
                    control = list(fnscale = -1))
    nac_Node[[2]] <- set_theta_from_unbounded(nac_Node, result$par)[[1]]
    param_list <- result$par
    print(param_list)
  }
  return(list(nac_Node, param_list, result))
}
