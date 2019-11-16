#' nested copula with U and inner copulas.. get_theta for list with length of 4
#' @export
get_theta <- function(U, naclist) {
  theta <- 0
  family <- naclist[[1]]
  empirical_taus <- cor(U[,naclist[[3]]], method = "kendall")

  n <- length(naclist[[3]])
  corr_coef <- c() # vector of correlation coefficients
  for (index in 1:(n-1)) {
    corr_coef <- c(corr_coef, empirical_taus[index, (index+1):n])
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


#' fitNCopula function
#' @export
fitNCopula <- function(U, naclist) {
  param_list <- c()

  if (length(naclist) == 4) {
    subcopulas_list <- naclist[[4]]
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))

    # get initial parameter estimates for each descendant copula
    for (i in 1:length(subcopulas_list)) {
      list <- subcopulas_list[[i]]
      temp <- fitNCopula(U, list)
      naclist[[4]][[i]] <- temp[[1]] # update the naclist of the parent copula with new child-list
      list <- temp[[1]]
      param_list <- c(param_list, temp[[2]]) # save new param to be used for optim fcn later
      V[,i] <- pNC(U, list)  # estimate cdf using list with new param
    }

    # get initial parameter estimates for the current copula
    # get theta based on U and child copula's cdf,
    # then call optim to update current copula and all child copulas
    list_pNC <- naclist[-4]
    UV_concat <- cbind(U[, list_pNC[[3]]], V)
    list_pNC[[3]] <- seq_len(ncol(UV_concat))

    # get theta for list of length 4
    theta <- get_theta(UV_concat, list_pNC)
    naclist[[2]] <- theta
    param_list <- c(
      transform_par_unbounded(naclist[[2]], naclist)[[1]],
      param_list
    )

    result <- optim(par = param_list,
                    fn = loglikNC,
                    U = U,
                    naclist = naclist,
                    method ="BFGS",
                    control = list(fnscale = -1))
    naclist <- updatePars(
      transform_par(result$par, naclist)[[1]],
      naclist)[[1]]
    param_list <- result$par
  } else if (length(naclist) == 3) {
    # check theta
    theta <- get_theta(U, naclist)
    naclist[[2]] <- theta

    result <- optim(par = transform_par_unbounded(naclist[[2]], naclist)[[1]],
                    fn = loglikNC,
                    U = U,
                    naclist = naclist,
                    method ="BFGS",
                    control = list(fnscale = -1))
    naclist[[2]] <- transform_par(result$par, naclist)[[1]]
    param_list <- result$par
    print(param_list)
  }
  return(list(naclist, param_list, result))
}
