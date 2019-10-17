# nested copula with U and inner copulas.. get_theta for list with length of 4
get_theta <- function(U, naclist) {
  theta <- 0
  family <- naclist[[1]]
  if (length(naclist) == 3) {
      empirical_taus <- cor(U[,naclist[[3]]], method = "kendall")

      n <- length(naclist[[3]])
      corr_coef <- c() # vector of correlation coefficients
      for (index in 1:(n-1)) {
        corr_coef <- cbind(corr_coef, empirical_taus[index, (index+1):n])
      }
  }

  else if (length(naclist) == 4) {
      empirical_taus <- cor(U[,naclist[[3]]], naclist[[4]]), method = "kendall")
      corr_coef <- empirical_taus
  }

  tau_init <- mean(corr_coef)
  theta <- iTau(archmCopula(family = family), tau = tau_init)

  if (theta < -1) {
    if (family == "Ali") {
      theta <- -1
    } else if (family == "Clayton") {
        if (length(naclist) == 3 && length(naclist[[3]]) == 2) {
          theta <- -1
        } else if (length(naclist) == 4 && length(naclist[[3]]) + length(naclist[[4]]) == 2) {
          theta <- -1
        }
    }
  }
  else if (theta < 0) {
    if (family == "Clayton") {
      if (length(naclist) == 3 && length(naclist[[3]]) >= 3) {
        theta <- 0
      } else if (length(naclist) == 4 && length(naclist[[3]]) + length(naclist[[4]]) >= 3) {
        theta <- 0
      }
    }
  }
  else if (theta < 1) {
    if (family == "Joe" | family == "Gumbel") {
      theta <- 1.1
    }
  }
  else if (theta == 0) {
    if (family == "Frank") {
      theta <- 0.000001
    }
  }
  else if (theta > 1) {
    if (family == "Joe") {
      theta <- 1
    }
  }

  return(theta)
}

new_list <- list("Clayton", 1, 1:2, list(list("Frank", 2, 3:4), # NAC structure
                                         list("Joe", 3, 5:7, list(list("Frank", 4, 8:9))),
                                         list("Clayton", 5, 10:11)))

new_list_short <- list("Clayton", 1, 1:2, list(list("Frank", 2, 3:4)))

# fitNCopula function
fitNCopula <- function(U, naclist) {
  param_list <- c()

  if (length(naclist) == 4) {
    subcopulas_list <- naclist[[4]]
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

    for (i in 1:length(subcopulas_list)) {
      list <- subcopulas_list[[i]]
      temp <- fitNCopula(U, list)
      naclist[[4]][[i]] <- temp[[1]] # update the naclist of the parent copula with new child-list
      list <- temp[[1]]
      param_list <- c(param_list, temp[[2]]) # save new param to be used for optim fcn later
      V[,i] <- pNC(U, list)  # estimate cdf using list with new param
    }

    # check theta
    list_pNC <- naclist # get theta based on U and child copula's cdf
    list_pNC[[4]] <- V

    # get theta for list of length 4
    theta <- get_theta(U, list_pNC)
    naclist[[2]] <- theta
    param_list <- c(naclist[[2]], param_list)

    result <- optim(par = param_list,
                    fn = loglikNC,
                    U = U,
                    naclist = naclist,
                    method ="BFGS",
                    control = list(fnscale = -1))
    naclist <- updatePars(result$par, naclist)
    param_list <- result$par
  }

  else if (length(naclist) == 3) {
    # check theta
    get_theta(U, naclist)
    naclist[[2]] <- theta

    result <- optim(par = naclist[[2]],
                    fn = loglikNC,
                    U = U,
                    naclist = naclist,
                    method ="BFGS",
                    control = list(fnscale = -1))
    naclist[[2]] <- result$par
    param_list <- result$par
    print(param_list)
  }
  return(list(naclist, param_list, result))
}

family <- "Gumbel" # copula family
tau <- c(0.2, 0.4, 0.6, 0.8) # Kendall’s tau
th <- iTau(archmCopula(family), tau = tau) # corresponding parameters
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

fitNCopula(U, new_list)[[3]]
fitNCopula(U, new_list_short)[[3]]

