#' Calculate density of a nested archimedean copula
#'
#' @param U matrix of pseudo-observations
#' @param naclist list specifying nesting structure
#' @return vector of length nrow(U) with density values
dNC <- function(U, naclist, log = TRUE) {
    density <- 0
    family <- naclist[[1]]

    if (length(naclist) == 3) {

        if (family == 'Clayton') {
            density <- density + dCopula(U[,naclist[[3]]], claytonCopula(naclist[[2]], dim =length(naclist[[3]])), log = TRUE)
        } else if (family == 'Frank') {
            density <- density + dCopula(U[,naclist[[3]]], frankCopula(naclist[[2]], dim =length(naclist[[3]])), log = TRUE)
        } else if (family == 'Gumbel') {
            density <- density + dCopula(U[,naclist[[3]]], gumbelCopula(naclist[[2]], dim =length(naclist[[3]])), log = TRUE)
        } else if (family == 'Independence') {
            density <- density + dCopula(U[,naclist[[3]]], indepCopula(dim =length(naclist[[3]])), log = TRUE)
        } else if (family == 'Joe') {
            density <- density + dCopula(U[,naclist[[3]]], joeCopula(naclist[[2]], dim =length(naclist[[3]])), log = TRUE)
        } else if (family == 'Ali') {
            density <- density + dCopula(U[,naclist[[3]]], amhCopula(naclist[[2]], dim =length(naclist[[3]])), log = TRUE)
        }
    }

    else if (length(naclist) == 4) {
        subcopulas_list <- naclist[[4]]
        V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

        for (v_ind in 1:length(subcopulas_list)) {
            density <- density + dNC(U, subcopulas_list[[v_ind]], log = TRUE)
            V[,v_ind] <- pNC(U, subcopulas_list[[v_ind]])
        }

        X <- cbind(U[,naclist[[3]]], V)
        # check family
        if (family == 'Clayton') {
            density <- density + dCopula(X, claytonCopula(naclist[[2]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Frank') {
            density <- density + dCopula(X, frankCopula(naclist[[2]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Gumbel') {
            density <- density + dCopula(X, gumbelCopula(naclist[[2]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Independence') {
            density <- density + dCopula(X, indepCopula(dim = ncol(X)), log = TRUE)
        } else if (family == 'Joe') {
            density <- density + dCopula(X, joeCopula(naclist[[2]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Ali') {
            density <- density + dCopula(X, amhCopula(naclist[[2]], dim = ncol(X)), log = TRUE)
        }
    }

    if (log) {
      return(density)
    } else {
      return(exp(density))
    }
}
