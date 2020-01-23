#' calculate number of parameters

#' @param nac_Node an NAC node
#' @return number of parameters
get_num_par <- function(nac_Node) {
  num_par <- 0
  if (has_subcopulas(nac_Node)) {
    num_par <- num_par + 1
    sub_list <- get_subcopulas(nac_Node)
    
    for (i in 1:length(sub_list)) {
      num_par <- num_par + num_par(sub_list[[i]])
    }
  }
  
  else if (has_subcopulas(nac_Node) == FALSE) {
    num_par <- num_par + 1
  }
  
  return(num_par)
}

