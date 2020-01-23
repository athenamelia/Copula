#' get family of nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @return family
get_family <- function(nac_Node) {
  return(nac_Node[[1]])
}

#' get initial theta for nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @return initial theta
get_iniTheta <- function(nac_Node) {
  return(nac_Node[[2]])
}

#' get pseudo-obs for nested Archimedean copula 
#'
#' @param nac_Node a nested Archimedean copula 
#' @return pseudo-obs
get_U_indices <- function(nac_Node) {
  tryCatch(return(nac_Node[[3]]),
           error = function(e) return(NULL))
}

#' get subcopulas for nested Archimedean copula, NULL if no child 
#'
#' @param nac_Node a nested Archimedean copula 
#' @return subcopulas 
get_subcopulas <- function(nac_Node) {
  tryCatch(return(nac_Node[[4]]),
           error = function(e) return(NULL))
}

