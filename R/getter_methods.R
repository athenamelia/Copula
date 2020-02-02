#' get family of nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @return family
#' @export
get_family <- function(nac_Node) {
  return(nac_Node$family)
}

#' get initial theta for nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @return initial theta
#' @export
get_theta <- function(nac_Node) {
  return(nac_Node$theta)
}

#' get pseudo-obs for nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @return pseudo-obs
#' @export
get_U_indices <- function(nac_Node) {
  tryCatch(return(nac_Node$U_indices),
           error = function(e) return(NULL))
}

#' get subcopulas for nested Archimedean copula, NULL if no child
#'
#' @param nac_Node a nested Archimedean copula
#' @param index index of subcopula
#' @return subcopulas
#' @export
get_subcopula <- function(nac_Node, index) {
  tryCatch(return(nac_Node$subcopula[[index]]),
           error = function(e) return(NULL))
}

