#' Calculate density of a nested archimedean copula
#'
#' @param U matrix of pseudo-observations
#' @param nac_Node a nested A. copula
#' @return vector of length nrow(U) with density values
#' @export
dncopula <- function(nac_Node, U, log = TRUE) {
  density <- 0
  family <- get_family(nac_Node)
  theta <- get_theta(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  U_indices <- get_U_indices(nac_Node)

  if (has_subcopula(nac_Node) == FALSE) {
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

  else if (has_subcopula(nac_Node)) {
    subcopulas <- count_subcopula(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = subcopulas)

    for (v_ind in 1:subcopulas) {
      child_copula <- get_subcopula(nac_Node, v_ind)
      density <- density + dncopula(child_copula, U, log = TRUE)
      V[,v_ind] <- pncopula(child_copula, U)
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
