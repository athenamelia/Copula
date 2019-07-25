#' Update nested copula parameters
#'
#' @param par vector of copula parameters
#' @param naclist list specifying nesting structure
#' @param current_index index of copula
#' @return a list of updated naclist with new parameters and current index
updatePars <- function(par, naclist, current_index = 1) {
    if (length(naclist) == 3) {
        naclist[[1]] <- par[current_index]
        current_index <- current_index + 1

        for (child_index in 1:length(naclist[[3]])) {
            temp <- updatePars(par, naclist[[3]][[child_index]], current_index)
            naclist[[3]][[child_index]] <- temp[[1]]
            current_index <- temp[[2]]
        }
    }

    else if (length(naclist) == 2) {
        naclist[[1]] <- par[current_index]
        current_index <- current_index + 1
    }

    return(list(naclist, current_index))
}
