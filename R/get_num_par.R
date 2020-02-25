#' calculate number of parameters

#' @param nac_Node an NAC node
#' @return number of parameters
#' @export
get_num_par <- function(nac_Node) {
  num_par <- 0
  num_par <- num_par + 1

  if (has_subcopula(nac_Node)) {
    num_subcopula <- count_subcopula(nac_Node)

    for (i in seq_len(num_subcopula)) {
      num_par <- num_par + get_num_par(get_subcopula(nac_Node, i))
    }
  }

  return(num_par)
}
