## Define a nested copula
library(copula)
family <- "Clayton"
tau <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9) # Kendall's tau
th <- iTau(archmCopula(family), tau = tau) # corresponding parameters

dNC <- function(U, naclist, log = TRUE) {
  density <- 0
  if (length(naclist) == 2) {
    density <- density + dCopula(U[,naclist[[2]]], claytonCopula(naclist[[1]], dim =length(naclist[[2]])), log = TRUE)
  }
  
  else if (length(naclist) == 3) {
    subcopulas_list <- naclist[[3]]
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))
    
    for (v_ind in 1:length(subcopulas_list)) {
      density <- density + dNC(U, subcopulas_list[[v_ind]])
      V[,v_ind] <- pNC(U, subcopulas_list[[v_ind]])
    }
    
    X <- cbind(U[,naclist[[2]]], V)
    density <- density + dCopula(X, claytonCopula(naclist[[1]], dim = ncol(X)), log = TRUE)
  }
  
  if(log) {
    return(density)  
  } else {
    return(exp(density))
  }
}

# no children
nlist_no_children <- list(th[1], 1:2)

dNC(U[,1:2], nlist_no_children)
density_no_children <- dCopula(U[,1:2], claytonCopula(th[1], dim = 2), log = TRUE)
density_no_children == dNC(U[,1:2], nlist_no_children)

# list with NULL
nlist_null <- list(th[1], NULL, list(list(th[2], 1:2), 
                                     list(th[3], 3:4)))
dNC(U[,1:4], nlist_null)

density_child1 <- dCopula(U[,1:2], claytonCopula(th[2], dim = 2), log = TRUE)
density_child2 <- dCopula(U[,3:4], claytonCopula(th[3], dim = 2), log = TRUE)

cc1 <- claytonCopula(th[2], dim = 2)
cdf_child1 <- pCopula(U[,1:2], copula = cc1)
cc2 <- claytonCopula(th[3], dim = 2)
cdf_child2 <- pCopula(U[,3:4], copula = cc2)

density_null <- dCopula(cbind(cdf_child1, cdf_child2), claytonCopula(th[1], dim = 2), log = TRUE)

# wrong_density dNC(U[,1:4], nlist_null) == density_null
dNC(U, nlist_null) == density_null + density_child1 + density_child2

# simple nesting structure
nlist <- list(th[1], 1:2, list(list(th[2], 3:4)))

dNC(U[,1:4], nlist)

density_child <- dCopula(U[,3:4], claytonCopula(th[2], dim = 2), log = TRUE)
cc <- claytonCopula(th[2], dim = 2)
cdf_child <- pCopula(U[,3:4], copula = cc)

density_simple <- dCopula(cbind(U[,1:2], cdf_child), claytonCopula(th[1], dim = 3), log = TRUE)
dNC(U[,1:4], nlist) == density_simple + density_child

###
# complicated nesting structure
nlist_very_long <- list(th[1], NULL, list(list(th[2], 1, list(list(th[3], 2:3))), # NAC structure
                                          list(th[4], 4, list(list(th[5], 5:6))),
                                          list(th[6], 7, list(list(th[7], 8:11)))))

# density of nested children copulas
d_test_nested_cplex1 <- dCopula(U[,2:3], claytonCopula(th[3], dim = 2), log = TRUE)
d_test_nested_cplex2 <- dCopula(U[,5:6], claytonCopula(th[5], dim = 2), log = TRUE)
d_test_nested_cplex3 <- dCopula(U[,8:11], claytonCopula(th[7], dim = 4), log = TRUE)

# cdf of nested children copulas
cc_nested_cplex1 <- claytonCopula(th[3], dim = 2)
cdf_nested_child1 <- pCopula(U[,2:3], copula = cc_nested_cplex1)
cc_nested_cplex2 <- claytonCopula(th[5], dim = 2)
cdf_nested_child2 <- pCopula(U[,5:6], copula = cc_nested_cplex2)
cc_nested_cplex3 <- claytonCopula(th[7], dim = 4)
cdf_nested_child3 <- pCopula(U[,8:11], copula = cc_nested_cplex3)

# density of children copulas
d_test_cplex1 <- dCopula(cbind(U[,1], cdf_nested_child1), claytonCopula(th[2], dim = 2), log = TRUE)
d_test_cplex2 <- dCopula(cbind(U[,4], cdf_nested_child2), claytonCopula(th[4], dim = 2), log = TRUE)
d_test_cplex3 <- dCopula(cbind(U[,7], cdf_nested_child3), claytonCopula(th[6], dim = 2), log = TRUE)

# cdf of children copulas
cc_cplex1 <- claytonCopula(th[2], dim = 2)
cdf_child1 <- pCopula(cbind(U[,1], cdf_nested_child1), copula = cc_cplex1)
cc_cplex2 <- claytonCopula(th[4], dim = 2)
cdf_child2 <- pCopula(cbind(U[,4], cdf_nested_child2), copula = cc_cplex2)
cc_cplex3 <- claytonCopula(th[6], dim = 2)
cdf_child3 <- pCopula(cbind(U[,7], cdf_nested_child3), copula = cc_cplex3)

# density of big tree
d_test_parent <- dCopula(cbind(cdf_child1, cdf_child2, cdf_child3), claytonCopula(th[1], dim = 3), log = TRUE)

d_test <-  (d_test_nested_cplex1 + d_test_cplex1 + (d_test_nested_cplex2 + d_test_cplex2) + (d_test_nested_cplex3 + d_test_cplex3) + d_test_parent)
# check 
d_test == dNC(U[,1:11], nlist_very_long) 

### another check 
nlist <- list(th[1], 1:2, list(list(th[2], 3:4), # NAC structure
                               list(th[3], 5, list(list(th[4], 6:7)))))
d1 <- dCopula(U[,6:7], claytonCopula(th[4], dim = 2), log = TRUE)

cc1 <- claytonCopula(th[4], dim = 2)
cdf1 <- pCopula(U[,6:7], copula = cc1)

d2 <- dCopula(cbind(U[,5], cdf1), claytonCopula(th[3], dim = 2), log = TRUE)
cc2 <- claytonCopula(th[3], dim = 2)
cdf2 <- pCopula(cbind(U[,5], cdf1), copula = cc2)

d3 <- dCopula(U[,3:4], claytonCopula(th[2], dim = 2), log = TRUE)
cc3 <- claytonCopula(th[2], dim = 2)
cdf3 <- pCopula(U[,3:4], copula = cc3)

d4 <- dCopula(cbind(U[,1:2], cdf2, cdf3), claytonCopula(th[1], dim = 4), log = TRUE)

d <- d1 + d2 + d3 + d4

d == dNC(U, nlist)