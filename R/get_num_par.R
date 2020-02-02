#' calculate number of parameters

#' @param nac_Node an NAC node
#' @return number of parameters
#' @export
get_num_par <- function(nac_Node) {
  num_par <- 0
  if (has_subcopula(nac_Node)) {
    num_par <- num_par + 1
    num_subcopula <- count_subcopula(nac_Node)

    for (i in 1:num_subcopula) {
      num_par <- num_par + get_num_par(get_subcopula(nac_Node, i))
    }
  }

  else if (has_subcopula(nac_Node) == FALSE) {
    num_par <- num_par + 1
  }

  return(num_par)
}
