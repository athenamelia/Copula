#' Update nested copula parameters
#'
#' @param nac_Node an NAC node
#' @param par a vector of new parameters
#' @param current_index index of copula
#' @return an updated nac_Node with new parameters and current index
#' @export
update_par <- function(nac_Node, par, current_index = 1) {
  if (has_subcopula(nac_Node)) {
    # if(is.na(nac_Node[[2]])) {
    #     cat("in updatePars")
    #     browser()
    # }

    nac_Node <- set_theta(nac_Node, par[current_index])
    current_index <- current_index + 1
    subcopulas <- count_subcopula(nac_Node)
    for (child_index in 1:subcopulas) {
      temp <- update_par(get_subcopula(nac_Node, child_index), par, current_index)
      nac_Node <- set_subcopula(nac_Node, child_index, temp[[1]])
      current_index <- temp[[2]]
    }
  }

  else if (has_subcopula(nac_Node) == FALSE) {
    nac_Node <- set_theta(nac_Node, par[current_index])
    # if(is.na(nac_Node[[2]])) {
    #     cat("in updatePars")
    #     browser()
    # }
    current_index <- current_index + 1
  }

  return(list(nac_Node, current_index))
}
