#' Transform parameters so that they are within the bounds
#'
#' @param par vector of copula parameters
#' @param naclist list specifying nesting structure
#' @param current_index index of copula
#' @return a list of updated parameters within the bounds and current index

transform_par <- function(par, naclist, current_index = 1) {
    family <- naclist[[1]]

    if (family == "Clayton") {
        if (length(naclist) == 3) {
            if (length(naclist[[3]]) == 2) {
                par[current_index] <- -1 + exp(par[current_index])
            } else {
                par[current_index] <- exp(par[current_index])
            }
            current_index <- current_index + 1
        }

        else if (length(naclist) == 4) {
            if (length(naclist[[4]]) + length(naclist[[3]]) == 2) {
                par[current_index] <- -1 + exp(par[current_index])
            } else {
                par[current_index] <- exp(par[current_index])
            }
            current_index <- current_index + 1

            subcopulas_list <- naclist[[4]]
            for (i in 1:length(subcopulas_list)) {
                list <- subcopulas_list[[i]]
                temp <- transform_par(par, list, current_index)
                par[current_index] <- temp[[1]][[current_index]]
                current_index <- temp[[2]]
            }
        }
    }

   if (family == "Frank") {
       for (index in 1:length(par)) {
           if (par[index] == 0) {
                par[index] = 0.000001
            }
        }
    }

    if (family == "Gumbel" | family == "Joe"){
        par <- exp(par) + 1
    }

    if (family == "Ali") {
        for (index in 1:length(par)) {
            par[index] <- 2 * (exp(par[index])/1+exp(par[index])) -1
        }
    }

    return(list(par, current_index))
}

