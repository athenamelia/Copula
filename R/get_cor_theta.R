#' Get theta for nested copula via U correlation
#'
#' @param nac_Node a nested A. copula
#' @param U pseudo-observation
#' @return theta for nested copula based on U correlation
#' @export
get_cor_theta <- function(nac_Node, U) {
  theta <- 0
  corr_coef <- c()
  family <- get_family(nac_Node)
  U_indices <- get_U_indices(nac_Node)

  if (is.null(U_indices) && count_subcopula(nac_Node) < 2) {
    return(NA)
  }

  if (has_subcopula(nac_Node) == FALSE) {
    empirical_taus <- cor(U[,U_indices], method = "kendall")
    n <- length(U_indices)
  }

  else if (has_subcopula(nac_Node)) {
    num_subcopulas <- count_subcopula(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = num_subcopulas)

    for (i in seq_len(num_subcopulas)) {
      child_copula <- get_subcopula(nac_Node, i)
      V[,i] <- pncopula(child_copula, U)
    }

    if (is.null(U_indices)) {
      empirical_taus <- cor(V, method = "kendall")
      n <- ncol(V)
    }

    else {
      UV <- cbind(U[,U_indices], V)
      empirical_taus <- cor(UV, method = "kendall")
      n <- length(U_indices) + ncol(V)
    }
  }

  for (index in 1:(n-1)) {
    corr_coef <- c(corr_coef, empirical_taus[index, (index+1):n])
  }

  tau <- mean(corr_coef)
  theta <- iTau(archmCopula(family = family), tau = tau)
  theta <- bound_theta(nac_Node, theta)
  return(theta)
}


#' A helper method to ensure the bound of NAC parameters
#'
#' @param nac_Node a nested Archimedean copula
#' @param theta a parameter to check
#' @return a valid theta
#' @export
bound_theta <- function(nac_Node, theta) {
  family <- get_family(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  nsubcopula <- count_subcopula(nac_Node)
  upper_bound <- Inf
  lower_bound <- -Inf
  dimension <- 3
  disallow_0 <- TRUE

  if (has_subcopula(nac_Node) == FALSE) {
    if (ncol_U == 2) {
      dimension <- 2
    }
  } else if (has_subcopula(nac_Node)) {
    if (nsubcopula + ncol_U == 2) {
      dimension <- 2
    }
  }

  if (family == "Clayton") {
    if (dimension == 2) {
      lower_bound <- -1
    } else {
      lower_bound <- 0.05
    }
  }

  else if (family == "Frank") {
    lower_bound <- 0.05
  }

  else if (family == "Gumbel" | family == "Joe") {
    lower_bound <- 1
  }

  else if (family == "Amh") {
    upper_bound <- 1
    disallow_0 <- FALSE
    if (dimension == 2) {
      lower_bound <- -1
    } else {
      lower_bound <- 0
    }
  }

  if (theta < lower_bound && lower_bound != -Inf) {
    theta <- lower_bound
  } else if (theta > upper_bound && upper_bound != Inf) {
    theta <- upper_bound
  }

  if (disallow_0 && (abs(theta) < sqrt(.Machine$double.eps))) {
    theta <- 2 * sqrt(.Machine$double.eps)
  }

  return(theta)
}

