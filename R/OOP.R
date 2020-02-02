nac_Node <- structure(list(), class = "nac_node")
nac_Node <- list()
class(nac_Node) <- "nac_node"

#' boolean function to check if it's an NAC node
#'
#' @param nac_Node a nested Archimedean copula
#' @return true false
#' @export
is.nac_Node <- function(nac_Node) {
  if (inherits(nac_Node, "nac_node")) { return(TRUE) }
  else {return(FALSE)}
}

#' create an NAC node
#'
#' @param family a nested Archimedean copula family
#' @param theta initial copula theta
#' @param U_indices indices of pseudo-observation
#' @param subcopula child copulas
#' @return NAC node
#' @export
new_nac_node <- function(family = character(), theta = double(), U_indices = double(), subcopula = NULL) {
  if(length(subcopula)  != 0) {
    for (i in 1:length(subcopula)) {
      stopifnot(is.nac_Node(subcopula[[i]]))
    }
  }
  stopifnot(is.character(family))
  stopifnot(family == "Clayton" || family == "Frank" || family == "Joe" ||
              family == "Ali" || family == "Independence" || family == "Gumbel")
  stopifnot(is.double(theta))
  stopifnot(is.double(U_indices) || is.integer(U_indices) || is.null(U_indices))
  stopifnot(is.list(subcopula))

  nac_Node <- list(family = family, theta = theta, U_indices = U_indices, subcopula = subcopula)
  structure(nac_Node, class = "nac_node")
}

