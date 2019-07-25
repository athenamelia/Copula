#' Calculate density of a nested archimedean copula
#'
#' @param U matrix of pseudo-observations
#' @param naclist list specifying nesting structure
#' @return vector of length nrow(U) with density values
dNC <- function(U, naclist, family) {
    density <- 0
    if (length(naclist) == 2) {

        if (family == 'Clayton') {
            density <- density + dCopula(U[,naclist[[2]]], claytonCopula(naclist[[1]], dim =length(naclist[[2]])), log = TRUE)
        } else if (family == 'Frank') {
            density <- density + dCopula(U[,naclist[[2]]], frankCopula(naclist[[1]], dim =length(naclist[[2]])), log = TRUE)
        } else if (family == 'Gumbel') {
            density <- density + dCopula(U[,naclist[[2]]], gumbelCopula(naclist[[1]], dim =length(naclist[[2]])), log = TRUE)
        } else if (family == 'Independence') {
            density <- density + dCopula(U[,naclist[[2]]], indepCopula(dim =length(naclist[[2]])), log = TRUE)
        } else if (family == 'Joe') {
            density <- density + dCopula(U[,naclist[[2]]], joeCopula(naclist[[1]], dim =length(naclist[[2]])), log = TRUE)
        }
    }

    else if (length(naclist) == 3) {
        subcopulas_list <- naclist[[3]]
        V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

        for (v_ind in 1:length(subcopulas_list)) {
            density <- density + dNC(U, subcopulas_list[[v_ind]], family)
            V[,v_ind] <- pNC(U, subcopulas_list[[v_ind]], family)
        }

        X <- cbind(U[,naclist[[2]]], V)
        # check family
        if (family == 'Clayton') {
            density <- density + dCopula(X, claytonCopula(naclist[[1]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Frank') {
            density <- density + dCopula(X, frankCopula(naclist[[1]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Gumbel') {
            density <- density + dCopula(X, gumbelCopula(naclist[[1]], dim = ncol(X)), log = TRUE)
        } else if (family == 'Independence') {
            density <- density + dCopula(X, indepCopula(dim = ncol(X)), log = TRUE)
        } else if (family == 'Joe') {
            density <- density + dCopula(X, joeCopula(naclist[[1]], dim = ncol(X)), log = TRUE)
        }
    }
    return(density)
}
