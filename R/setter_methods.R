#' set family for nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @param new_family new family
#' @export
set_family <- function(nac_Node, new_family) {
  if (is.character(new_family)) {
    if (new_family == "Clayton" ||
        new_family == "Frank" ||
        new_family == "Gumbel" ||
        new_family == "Joe" ||
        new_family == "Independence" ||
        new_family == "Amh") {
          nac_Node$family <- new_family
          return(nac_Node)
    } else {
      print("invalid input")
      return(NA)
    }
  }

  else {
    print("invalid input")
    return(NA)
  }
}

#' set new initial theta for nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @param new_theta new theta
#' @export
set_theta <- function(nac_Node, new_theta) {
  if (is.double(new_theta)) {
      nac_Node$theta <- new_theta
      return(nac_Node)
  } else {
      print("invalid input")
      return(NA)
  }
}

#' set new pseudo-obs for nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @param new_indices new indices of U
#' @export
set_U_indices <- function(nac_Node, new_indices) {
  if (is.double(new_indices) || is.integer(new_indices)) {
      nac_Node$U_indices <- new_indices
      return(nac_Node)
  } else {
      print("invalid input")
      return(NA)
  }
}

#' set new subcopula for nested Archimedean copula
#'
#' @param nac_Node a nested Archimedean copula
#' @param index index of subcopula
#' @param new_subcopula new subcopula
#' @export
set_subcopula <- function(nac_Node, index, new_subcopula) {
  if (is.nac_Node(new_subcopula)) {
      nac_Node$subcopula[[index]] <- new_subcopula
      return(nac_Node)
  } else {
      print("invalid input")
      return(NA)
  }
}

#' boolean function to check if NAC has a subcopula
#'
#' @param nac_Node a nested Archimedean copula
#' @export
has_subcopula <- function(nac_Node) {
  if (length(nac_Node$subcopula) != 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' add subcopula to a nac_Node
#'
#' @param nac_Node a nested Archimedean copula
#' @param new_subcopula new subcopula to be added
#' @export
append_subcopula <- function(nac_Node, new_subcopula) {
  if (has_subcopula(nac_Node) == FALSE) {
    nac_Node$subcopula <- append(nac_Node$subcopula, new_subcopula)
    return(nac_Node)
  } else {
    print("nac_node has a subcopula already")
    return(NA)
  }
}

#' count number of subcopulas
#'
#' @param nac_Node a nested Archimedean copula
#' @export
count_subcopula <- function(nac_Node) {
  count <- length(nac_Node$subcopula)
  return(count)
}

