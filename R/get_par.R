#' get all parameters of a nested Archimedean copula
#'
#' @param nac_Node an NAC node
#' @return par: a vector of NAC parameters
#' @export
get_par <- function(nac_Node) {
  par <- c()
  par <- c(par, get_theta(nac_Node))
  if (has_subcopula(nac_Node)) {
    subcopulas <- count_subcopula(nac_Node)
    for (child_index in seq_len(subcopulas)) {
      child_copula <- get_subcopula(nac_Node, child_index)
      par <- c(par, get_par(child_copula))
    }
  }

  return(par)
}
