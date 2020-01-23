#' Update nested copula parameters
#'
#' @param nac_Node an NAC node
#' @param par a vector of new parameters
#' @param current_index index of copula
#' @return an updated nac_Node with new parameters and current index
update_par <- function(nac_Node, par, current_index = 1) {
  if (has_subcopulas(nac_Node)) {
    # set_iniTheta(nac_Node, par[current_index])
    # if(is.na(nac_Node[[2]])) {
    #     cat("in updatePars")
    #     browser()
    # }
    
    nac_Node[[2]] <- par[current_index]
    current_index <- current_index + 1
    subcopulas <- get_subcopulas(nac_Node)
    for (child_index in 1:length(subcopulas)) {
      temp <- update_par(subcopulas[[child_index]], par, current_index)
      # set_subcopula(nac_Node[[4]][[child_index]], temp[[1]])
      nac_Node[[4]][[child_index]] <- temp[[1]]
      current_index <- temp[[2]]
    }
  }
  
  else if (has_subcopulas(nac_Node) == FALSE) {
    # set_iniTheta(nac_Node, par[current_index])
    nac_Node[[2]] <- par[current_index]
    # if(is.na(nac_Node[[2]])) {
    #     cat("in updatePars")
    #     browser()
    # }
    current_index <- current_index + 1
  }
  
  return(list(nac_Node, current_index))
}
