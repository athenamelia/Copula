#' get initial theta for nested copula with U
#'
#' @param nac_Node a nested A. copula
#' @param U pseudo-observation
#' @return initial theta for nested copula based on U correlation
#' @export
get_iniCor <- function(nac_Node, U) {
  theta <- 0
  family <- get_family(nac_Node)
  U_indices <- get_U_indices(nac_Node)
  empirical_taus <- cor(U[,U_indices], method = "kendall")

  n <- length(U_indices)
  corr_coef <- c() # vector of correlation coefficients
  for (index in 1:(n-1)) {
    corr_coef <- c(corr_coef, empirical_taus[index, (index+1):n])
  }

  tau_init <- mean(corr_coef)
  theta <- iTau(archmCopula(family = family), tau = tau_init)

  if (theta < -1) {
    if (family == "Ali") {
      theta <- -1
    } else if (family == "Clayton") {
      if (has_subcopulas(nac_Node) == FALSE && length(get_U_indices(nac_Node)) == 2) {
        theta <- -1
      } else if (has_subcopulas(nac_Node) && length(get_U_indices(nac_Node)) + length(get_subcopulas(nac_Node)) == 2) {
        theta <- -1
      }
    }
  }
  else if (theta < 0) {
    if (family == "Clayton") {
      if (has_subcopulas(nac_Node) == FALSE && length(get_U_indices(nac_Node)) >= 3) {
        theta <- 0
      } else if (has_subcopulas(nac_Node) && length(get_U_indices(nac_Node)) + length(get_subcopulas(nac_Node)) >= 3) {
        theta <- 0
      }
    }
  }
  else if (theta <= 1) {
    if (family == "Joe" | family == "Gumbel") {
      theta <- 1.000001
    }
  }
  else if (theta == 0) {
    if (family == "Frank") {
      theta <- 0.000001
    }
  }
  else if (theta > 1) {
    if (family == "Joe") {
      theta <- 1
    }
  }

  return(theta)
}

