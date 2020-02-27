#' set nested copula with new parameters
#'
#' @param nac_Node an NAC node
#' @param par a vector of new parameters
#' @param current_index index of copula
#' @return an nac_Node with new parameters and current index
#' @export
set_par <- function(nac_Node, theta, current_index = 1) {
  nac_Node <- set_theta(nac_Node, theta[current_index])
  current_index <- current_index + 1

  if (has_subcopula(nac_Node)) {
    subcopulas <- count_subcopula(nac_Node)
    for (child_index in seq_len(subcopulas)) {
      temp <- set_par(get_subcopula(nac_Node, child_index), theta, current_index)
      nac_Node <- set_subcopula(nac_Node, child_index, temp[[1]])
      current_index <- temp[[2]]
    }
  }

  return(list(nac_Node, current_index))
}
