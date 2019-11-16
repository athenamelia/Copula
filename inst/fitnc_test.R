library(ncopula)
library(copula)

new_list <- list("Clayton", 1, 1:2, list(list("Frank", 2, 3:4), # NAC structure
                                         list("Joe", 3, 5:7, list(list("Frank", 4, 8:9))),
                                         list("Clayton", 5, 10:11)))

new_list_short <- list("Clayton", 1, 1:2, list(list("Frank", 2, 3:4)))

tau <- c(0.2, 0.4, 0.6, 0.8) # Kendall’s tau
th <- iTau(archmCopula(family), tau = tau) # corresponding parameters



### Evan worked with examples below here:
family <- "Gumbel" # copula family
nlist <- list(1, 1:2, list(list(2, 3:4), # NAC structure
                           list(3, 5:7, list(list(4, 8:9))),
                           list(5, 10:11)))

new_list <- list("Gumbel", 1, 1:2, list(list("Gumbel", 2, 3:4), # NAC structure
                                        list("Gumbel", 3, 5:7, list(list("Gumbel", 4, 8:9))),
                                        list("Gumbel", 5, 10:11)))

NAC <- onacopulaL(family, nacList = nlist) # NAC copula
## Sample
set.seed(271) # set seed (for reproducibility)
U <- rCopula(1000, copula = NAC) # sample
pairs(U)

fitNCopula(U, new_list)[[1]]



### Fitting an inner copula
U <- U[, 3:4]
inner_naclist <- new_list[[4]][[1]]
inner_naclist[[3]] <- 1:2
theta <- get_theta(U, inner_naclist)
inner_naclist[[2]] <- theta

result <- optim(par = transform_par_unbounded(inner_naclist[[2]], inner_naclist)[[1]],
                fn = loglikNC,
                U = U,
                naclist = inner_naclist,
                method ="BFGS",
                control = list(fnscale = -1))

transform_par(result$par, inner_naclist)[[1]]


fitNCopula(U, new_list_short)[[3]]
