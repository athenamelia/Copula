#' Transform vector of unbounded real numbers (as optim might use)
#' to a vector of parameters that are used for each copula in a NAC.
#'
#' @param nac_Node a nested Archimedean copula
#' @return a vector of transformed parameters
#' @export

transform_theta_unbounded_to_bounded <- function(nac_Node) {
  param <- c()
  param <- c(param, transform_nac(nac_Node))

  if (has_subcopula(nac_Node)) {
    for (i in seq_len(count_subcopula(nac_Node))) {
      child_copula <- get_subcopula(nac_Node, i)
      param <- c(param, transform_theta_unbounded_to_bounded(child_copula))
    }
  }
  return(param)
}


#' A helper method to transform vector of unbounded real numbers (as optim might use)
#' to a vector of parameters that are used for each copula in a NAC.
#'
#' @param nac_Node a nested Archimedean copula
#' @return a bounded theta
#' @export
transform_nac <- function(nac_Node) {
  family <- get_family(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  nsubcopula <- count_subcopula(nac_Node)
  theta <- get_theta(nac_Node)
  upper_bound <- 0
  lower_bound <- 0
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
    upper_bound <- Inf
    if (dimension == 2) {
      lower_bound <- -1 # cannot be 0
    } else {
      lower_bound <- 0
    }
  }

  else if (family == "Frank") {
    upper_bound <- Inf
    if (dimension == 2) {
      lower_bound <- -Inf # cannot be 0
    } else {
      lower_bound <- 0
    }
  }

  else if (family == "Gumbel" | family == "Joe") {
    upper_bound <- Inf
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

  if (lower_bound > -Inf & upper_bound == Inf) {
    theta <- lower_bound + exp(theta)
  } else if (lower_bound > -Inf & upper_bound < Inf) {
    theta <- lower_bound + exp(theta)/(1 + exp(theta)) * (upper_bound - lower_bound)
  }

  if (disallow_0 && (abs(theta) < sqrt(.Machine$double.eps))) {
    # Frank and Clayton families can’t have theta = 0
    theta <- 2 * sqrt(.Machine$double.eps)
  }
  return(theta)
}
