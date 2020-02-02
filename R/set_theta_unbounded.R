#' transform parameters so that they are unbounded
#'
#' @param nac_Node a nested A. copula
#' @param par vector of copula parameters
#' @param current_index index of copula
#' @return a list of unbounded parameters and current index
#' @export
set_theta_unbounded <- function(nac_Node, par, current_index = 1) {
  family <- get_family(nac_Node)

  if (family == "Clayton") {
    if (has_subcopula(nac_Node) == FALSE) {
      if (length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- log(par[current_index] + 1)
      } else {
        par[current_index] <- log(par[current_index])
      }
      current_index <- current_index + 1
    }

    else if (has_subcopula(nac_Node)) {
      if (count_subcopula(nac_Node) + length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- log(par[current_index] + 1)
      } else {
        par[current_index] <- log(par[current_index])
      }
      current_index <- current_index + 1

      subcopulas <- count_subcopula(nac_Node)
      for (i in 1:subcopulas) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }

  if (family == "Frank") {
    # do nothing for frank family
  }

  if (family == "Gumbel" | family == "Joe"){
    if (has_subcopula(nac_Node) == FALSE) {
      par[current_index] <- log(par[current_index] - 1)
      current_index <- current_index + 1
    }
    else if (has_subcopula(nac_Node)) {
      par[current_index] <- log(par[current_index] - 1)
      current_index <- current_index + 1

      subcopulas_list <- get_subcopula(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }

  if (family == "Ali") {
    if (has_subcopula(nac_Node) == FALSE) {
      par[current_index] <- 1 - 1/ (1 - (par[current_index] + 1)/2)
      current_index <- current_index + 1
    }
    else if (has_subcopula(nac_Node)) {
      par[current_index] <- 1 - 1/ (1 - (par[current_index] + 1)/2)
      current_index <- current_index + 1

      subcopulas_list <- get_subcopula(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        child_copula <- subcopulas_list[[i]]
        temp <- set_theta_unbounded(child_copula, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  return(list(par, current_index))
}
