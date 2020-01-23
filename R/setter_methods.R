#' set family for nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @param family new family
set_family <- function(nac_Node, family) {
  nac_Node[[1]] <- family
  # return(nac_Node)
}

#' set new initial theta for nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @param theta new theta
set_iniTheta <- function(nac_Node, theta) {
  nac_Node[[2]] <- theta
}

#' set new pseudo-obs for nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @param index new indices of U
set_U_indices <- function(nac_Node, index) {
  nac_Node[[3]] <- index
}

#' set new subcopula for nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @param subcopula new subcopula 
set_subcopula <- function(nac_Node, subcopula) {
  nac_Node[[4]] <- subcopula
}

#' boolean function to check if NAC has a subcopula
#'
#' @param nac_Node a nested Archimedean copula 
has_subcopulas <- function(nac_Node) {
  if (length(nac_Node) == 4) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' add subcopula to a nac_Node
#'
#' @param nac_Node a nested Archimedean copula 
#' @param subcopula new subcopula to be added
append_subcopula <- function(nac_Node, subcopula) {
  nac_Node <- append(nac_Node, subcopula)
  return(nac_Node)
}
