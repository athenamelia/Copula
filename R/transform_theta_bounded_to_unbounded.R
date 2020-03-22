#' Transform vector of parameters that are used for each copula in a NAC
#' to a vector of unbounded real numbers (as optim might use).
#'
#' @param nac_Node a nested Archimedean copula
#' @return a vector of unbounded parameters
#' @export

transform_theta_bounded_to_unbounded <- function(nac_Node) {
  param <- c()
  param <- c(param, untransform_nac(nac_Node))

  if (has_subcopula(nac_Node)) {
    for (i in seq_len(count_subcopula(nac_Node))) {
      child_copula <- get_subcopula(nac_Node, i)
      param <- c(param, transform_theta_bounded_to_unbounded(child_copula))
    }
  }
  return(param)
}

#' A helper method to transform vector of NAC parameters
#' to a vector of unbounded real numbers (as optim might use).
#'
#' @param nac_Node a nested Archimedean copula
#' @return an unbounded theta
#' @export
untransform_nac <- function(nac_Node) {
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
      lower_bound <- -1
    } else {
      lower_bound <- 0
    }
  }

  else if (family == "Frank") {
    upper_bound <- Inf
    if (dimension == 2) {
      lower_bound <- -Inf
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
  print(paste0("in untransform function, original theta is ", theta))
  print(paste0("in untransform function, the bounds are ", lower_bound, " and ", upper_bound))

  if (lower_bound > -Inf & upper_bound == Inf) {
    theta <- logspace_sub(theta - lower_bound, 0)

  } else if (lower_bound > -Inf & upper_bound < Inf) {
    theta <- log(1 / (1 - (theta - lower_bound) / (upper_bound - lower_bound)) - 1)
  }

  if (disallow_0 && (abs(theta) < sqrt(.Machine$double.eps))) {
    theta <- 2 * sqrt(.Machine$double.eps)
  }

  if (theta == -Inf || is.nan(theta)) {
    theta <- -.Machine$double.xmax
  }

  print(paste0("in untransform function, theta is ", theta))
  return(theta)
}
