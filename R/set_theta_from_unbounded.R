#' transform parameters of NAC so that they are within the bounds
#'
#' @param nac_Node a nested A. copula
#' @param par vector of copula parameters
#' @param current_index index of copula
#' @return a list of updated parameters within the bounds and current index
set_theta_from_unbounded <- function(nac_Node, par, current_index = 1) {
  family <- get_family(nac_Node)
  
  if (family == "Clayton") {
    if (has_subcopulas(nac_Node) == FALSE) {
      if (length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- -1 + exp(par[current_index])
      } else {
        par[current_index] <- exp(par[current_index])
      }
      current_index <- current_index + 1
    }
    
    else if (has_subcopulas(nac_Node)) {
      if (length(get_subcopulas(nac_Node)) + length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- -1 + exp(par[current_index])
      } else {
        par[current_index] <- exp(par[current_index])
      }
      current_index <- current_index + 1
      
      subcopulas_list <- get_subcopulas(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  else if (family == "Frank") {
    if (has_subcopulas(nac_Node) == FALSE) {
      if (par[current_index] == 0) {
        par[current_index] = 0.000001
      }
      current_index <- current_index + 1
    }
    
    else if (has_subcopulas(nac_Node)) {
      if (par[current_index] == 0) {
        par[current_index] = 0.000001
      }
      
      current_index <- current_index + 1
      subcopulas_list <- get_subcopulas(nac_Node)
      
      for (i in 1:length(subcopulas_list)) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  else if (family == "Gumbel" | family == "Joe"){
    if (has_subcopulas(nac_Node) == FALSE) {
      par[current_index] <- exp(par[current_index]) + 1
      current_index <- current_index + 1
    }
    else if (has_subcopulas(nac_Node)) {
      par[current_index] <- exp(par[current_index]) + 1
      current_index <- current_index + 1
      subcopulas_list <- get_subcopulas(nac_Node)
      
      for (i in 1:length(subcopulas_list)) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  else if (family == "Ali") {
    if (has_subcopulas(nac_Node) == FALSE) {
      par[current_index] <- 2 * (exp(par[current_index])/1+exp(par[current_index])) -1
      current_index <- current_index + 1
    }
    else if (has_subcopulas(nac_Node)) {
      par[current_index] <- 2 * (exp(par[current_index])/1+exp(par[current_index])) -1
      current_index <- current_index + 1
      subcopulas_list <- get_subcopulas(nac_Node)
      
      for (i in 1:length(subcopulas_list)) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  return(list(par, current_index))
}

