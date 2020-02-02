#' Calculate cdf of a nested archimedean copula
#'
#' @param nac_Node a nested A. copula
#' @param U matrix of pseudo-observations
#' @return vector of length nrow(U) with cdf values
#' @export
pncopula <- function(nac_Node, U) {
  nestedCopula <- 0
  family <- get_family(nac_Node)
  theta <- get_theta(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  U_indices <- get_U_indices(nac_Node)

  if (has_subcopula(nac_Node)) {
    subcopulas <- count_subcopula(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = subcopulas)

    for (v_ind in 1:subcopulas) {
      V[,v_ind] <- pncopula(get_subcopula(nac_Node, v_ind), U)
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

  else if (has_subcopula(nac_Node) == FALSE) {
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


