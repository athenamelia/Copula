
#' Calculate cdf of a nested archimedean copula
#'
#' @param U matrix of pseudo-observations
#' @param naclist list specifying nesting structure
#' @return vector of length nrow(U) with cdf values
pNC <- function(U, naclist) {
    nestedCopula <- 0
    family <- naclist[[1]]

    if (length(naclist) == 4) {
        subcopulas_list <- naclist[[4]]
        V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

        for (v_ind in 1:length(subcopulas_list)) {
            V[,v_ind] <- pNC(U, subcopulas_list[[v_ind]])
        }

        # check family
        if (family == 'Clayton') {
            copula <- claytonCopula(naclist[[2]], dim = length(naclist[[3]]) + ncol(V))
        } else if (family == 'Frank') {
            copula <- frankCopula(naclist[[2]], dim = length(naclist[[3]]) + ncol(V))
        } else if (family == 'Gumbel') {
            copula <- gumbelCopula(naclist[[2]], dim = length(naclist[[3]]) + ncol(V))
        } else if (family == 'Independence') {
            copula <- indepCopula(dim = length(naclist[[3]]) + ncol(V))
        } else if (family == 'Joe') {
            copula <- joeCopula(naclist[[2]], dim = length(naclist[[3]]) + ncol(V))
        } else if (family == 'Ali') {
            copula <- amhCopula(naclist[[2]], dim = length(naclist[[3]]) + ncol(V))
        }

        nestedCopula <- nestedCopula + pCopula(cbind(U[,naclist[[3]]],V), copula = copula)
    }

    else if (length(naclist) == 3) {
        # check family
        if (family == 'Clayton') {
            copula <- claytonCopula(naclist[[2]], dim = length(naclist[[3]]))
        } else if (family == 'Frank') {
            copula <- frankCopula(naclist[[2]], dim = length(naclist[[3]]))
        } else if (family == 'Gumbel') {
            copula <- gumbelCopula(naclist[[2]], dim = length(naclist[[3]]))
        } else if (family == 'Independence') {
            copula <- indepCopula(dim = length(naclist[[3]]))
        } else if (family == 'Joe') {
            copula <- joeCopula(naclist[[2]], dim = length(naclist[[3]]))
        } else if (family == 'Ali') {
            copula <- amhCopula(naclist[[2]], dim = length(naclist[[3]]))
        }

        nestedCopula <- nestedCopula + pCopula(U[,naclist[[3]]], copula = copula)
    }
    return(nestedCopula)
}
