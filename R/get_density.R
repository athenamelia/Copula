#' Calculate density of a nested archimedean copula
#'
#' @param U matrix of pseudo-observations
#' @param nac_Node a nested A. copula
#' @return vector of length nrow(U) with density values
get_density <- function(nac_Node, U, log = TRUE) {
  density <- 0
  family <- get_family(nac_Node)
  theta <- get_iniTheta(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  U_indices <- get_U_indices(nac_Node)
  
  if (has_subcopulas(nac_Node) == FALSE) {
    if (family == 'Clayton') {
      density <- density + dCopula(U[,U_indices], claytonCopula(theta, dim = ncol_U), log = TRUE)
    } else if (family == 'Frank') {
      density <- density + dCopula(U[,U_indices], frankCopula(theta, dim = ncol_U), log = TRUE)
    } else if (family == 'Gumbel') {
      density <- density + dCopula(U[,U_indices], gumbelCopula(theta, dim = ncol_U), log = TRUE)
    } else if (family == 'Independence') {
      density <- density + dCopula(U[,U_indices], indepCopula(dim = ncol_U), log = TRUE)
    } else if (family == 'Joe') {
      density <- density + dCopula(U[,U_indices], joeCopula(theta, dim = ncol_U), log = TRUE)
    } else if (family == 'Ali') {
      density <- density + dCopula(U[,U_indices], amhCopula(theta, dim = ncol_U), log = TRUE)
    }
  }
  
  else if (has_subcopulas(nac_Node)) {
    subcopulas_list <- get_subcopulas(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))
    
    for (v_ind in 1:length(subcopulas_list)) {
      density <- density + get_density(subcopulas_list[[v_ind]], U, log = TRUE)
      V[,v_ind] <- get_cdf(subcopulas_list[[v_ind]], U)
    }
    
    X <- cbind(U[,U_indices], V)
    # check family
    if (family == 'Clayton') {
      density <- density + dCopula(X, claytonCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Frank') {
      density <- density + dCopula(X, frankCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Gumbel') {
      density <- density + dCopula(X, gumbelCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Independence') {
      density <- density + dCopula(X, indepCopula(dim = ncol(X)), log = TRUE)
    } else if (family == 'Joe') {
      density <- density + dCopula(X, joeCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Ali') {
      density <- density + dCopula(X, amhCopula(theta, dim = ncol(X)), log = TRUE)
    }
  }
  
  if (log) {
    return(density)
  } else {
    return(exp(density))
  }
}