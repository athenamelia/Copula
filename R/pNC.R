# Hello, world!  This is an example function named 'hello' which prints 'Hello, world!'.  You can learn more about
# package authoring with RStudio at: http://r-pkgs.had.co.nz/ Some useful keyboard shortcuts for package authoring:
# Build and Reload Package: 'Cmd + Shift + B' Check Package: 'Cmd + Shift + E' Test Package: 'Cmd + Shift + T'

#' Calculate cdf of a nested archimedean copula
#'
#' @param U matrix of pseudo-observations
#' @param naclist list specifying nesting structure
#' @return vector of length nrow(U) with cdf values
pNC <- function(U, naclist, family) {
    nestedCopula <- 0
    if (length(naclist) == 3) {
        subcopulas_list <- naclist[[3]]
        V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

        for (v_ind in 1:length(subcopulas_list)) {
            V[,v_ind] <- pNC(U, subcopulas_list[[v_ind]], family)
        }

        # check family
        if (family == 'Clayton') {
            copula <- claytonCopula(naclist[[1]], dim = length(naclist[[2]]) + ncol(V))
        } else if (family == 'Frank') {
            copula <- frankCopula(naclist[[1]], dim = length(naclist[[2]]) + ncol(V))
        } else if (family == 'Gumbel') {
            copula <- gumbelCopula(naclist[[1]], dim = length(naclist[[2]]) + ncol(V))
        } else if (family == 'Independence') {
            copula <- indepCopula(dim = length(naclist[[2]]) + ncol(V))
        } else if (family == 'Joe') {
            copula <- joeCopula(naclist[[1]], dim = length(naclist[[2]]) + ncol(V))
        }

        nestedCopula <- nestedCopula + pCopula(cbind(U[,naclist[[2]]],V), copula = copula)
    }

    else if (length(naclist) == 2) {
        # check family
        if (family == 'Clayton') {
            copula <- claytonCopula(naclist[[1]], dim = length(naclist[[2]]))
        } else if (family == 'Frank') {
            copula <- frankCopula(naclist[[1]], dim = length(naclist[[2]]))
        } else if (family == 'Gumbel') {
            copula <- gumbelCopula(naclist[[1]], dim = length(naclist[[2]]))
        } else if (family == 'Independence') {
            copula <- indepCopula(dim = length(naclist[[2]]))
        } else if (family == 'Joe') {
            copula <- joeCopula(naclist[[1]], dim = length(naclist[[2]]))
        }

        nestedCopula <- nestedCopula + pCopula(U[,naclist[[2]]], copula = copula)
    }
    return(nestedCopula)
}
