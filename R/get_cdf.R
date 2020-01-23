#' Calculate cdf of a nested archimedean copula
#'
#' @param nac_Node a nested A. copula
#' @param U matrix of pseudo-observations
#' @return vector of length nrow(U) with cdf values
get_cdf <- function(nac_Node, U) {
  nestedCopula <- 0
  family <- get_family(nac_Node)
  theta <- get_iniTheta(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  U_indices <- get_U_indices(nac_Node)
  
  if (has_subcopulas(nac_Node)) {
    subcopulas_list <- get_subcopulas(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))
    
    for (v_ind in 1:length(subcopulas_list)) {
      V[,v_ind] <- get_cdf(subcopulas_list[[v_ind]], U)
    }
    
    # check family
    if (family == 'Clayton') {
      copula <- claytonCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Frank') {
      copula <- frankCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Gumbel') {
      copula <- gumbelCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Independence') {
      copula <- indepCopula(dim = ncol_U + ncol(V))
    } else if (family == 'Joe') {
      copula <- joeCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Ali') {
      copula <- amhCopula(theta, dim = ncol_U + ncol(V))
    }
    
    nestedCopula <- nestedCopula + pCopula(cbind(U[,U_indices],V), copula = copula)
  }
  
  else if (has_subcopulas(nac_Node) == FALSE) {
    if (family == 'Clayton') {
      copula <- claytonCopula(theta, dim = ncol_U)
    } else if (family == 'Frank') {
      copula <- frankCopula(theta, dim = ncol_U)
    } else if (family == 'Gumbel') {
      copula <- gumbelCopula(theta, dim = ncol_U)
    } else if (family == 'Independence') {
      copula <- indepCopula(dim = ncol_U)
    } else if (family == 'Joe') {
      copula <- joeCopula(theta, dim = ncol_U)
    } else if (family == 'Ali') {
      copula <- amhCopula(theta, dim = ncol_U)
    }
    
    nestedCopula <- nestedCopula + pCopula(U[,U_indices], copula = copula)
  }
  return(nestedCopula)
}


