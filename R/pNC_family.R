## Define a nested copula
library(copula)
family <- "Clayton"
tau <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9) # Kendall's tau
th <- iTau(archmCopula(family), tau = tau) # corresponding parameters

# pNestedCopula with family specification 

pNC_update <- function(U, naclist, family) {
  nestedCopula <- 0
  if (length(naclist) == 3) {
    subcopulas_list <- naclist[[3]]
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))
    
    for (v_ind in 1:length(subcopulas_list)) {
      V[,v_ind] <- pNC_update(U, subcopulas_list[[v_ind]], family)
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

# test
nlist <- list(th[1], 1:2, list(list(th[2], 3:4)))

nlist_null <- list(th[1], NULL, list(list(th[2], 1:2), list(th[3], 3:4)))

nlist_no_children <- list(th[1], 1:2)

nlist_very_long <- list(th[1], NULL, list(list(th[2], 1, list(list(th[3], 2:3))), # NAC structure
                                          list(th[4], 4, list(list(th[5], 5:6))),
                                          list(th[6], 7, list(list(th[7], 8:11)))))

pNC(U[,1:2], nlist_no_children)
pNC(U[,1:4], nlist)
pNC(U[,1:4], nlist_null)

## no children
cc_no_children <- frankCopula(nlist_no_children[[1]], dim = 2)
p_test_no_children <- pCopula(U[,1:2], copula = cc_no_children)
p_test_no_children == pNC_update(U, nlist_no_children, family = 'Frank')

## normal list
nested_cc <- claytonCopula(th[2], dim = 2)
nested_p_test <- pCopula(U[,3:4], copula = nested_cc)
cc <- claytonCopula(th[1], dim = 3)
p_test <- pCopula(cbind(U[,1:2], nested_p_test), copula = cc)
p_test == pNC_update(U[,1:4], nlist, family = 'Clayton')


## list with NULL 
nested_cc2 <- claytonCopula(th[3], dim = 2)
nested_p_test2 <- pCopula(U[,3:4], copula = nested_cc2)
nested_cc1 <- claytonCopula(th[2], dim = 2)
nested_p_test1 <- pCopula(U[,1:2], copula = nested_cc1)
cc_null <- claytonCopula(th[1], dim = 2)
p_test_null <- pCopula(cbind(nested_p_test1, nested_p_test2), copula = cc_null)
p_test_null == pNC_update(U[,1:4], nlist_null, family = 'Clayton')

## 
nested_verylong_cc3 <- claytonCopula(th[7], dim = 3)
nested_verylong_p_test3 <- pCopula(U[,8:11], copula = nested_verylong_cc3)
verylong_cc3 <- claytonCopula(th[6], dim = 2)
p_test3 <- pCopula(cbind(U[,7], nested_verylong_p_test3), copula = verylong_cc3)

nested_verylong_cc2 <- claytonCopula(th[5], dim = 2)
nested_verylong_p_test2 <- pCopula(U[,5:6], copula = nested_verylong_cc2)
verylong_cc2 <- claytonCopula(th[4], dim = 2)
p_test2 <- pCopula(cbind(U[,4], nested_verylong_p_test2), copula = verylong_cc2)

nested_verylong_cc1 <- claytonCopula(th[3], dim = 2)
nested_verylong_p_test1 <- pCopula(U[,2:3], copula = nested_verylong_cc1)
verylong_cc1 <- claytonCopula(th[2], dim = 2)
p_test1 <- pCopula(cbind(U[,1], nested_verylong_p_test1), copula = verylong_cc1)

cc_verylong <- claytonCopula(th[1], dim = 3)
p_test_verylong <- pCopula(cbind(p_test1, p_test2, p_test3), copula = cc_verylong)
p_test_verylong == pNC_update(U[,1:11], nlist_very_long, family = 'Clayton')

# clayton Copula parameter negative iff dim = 2