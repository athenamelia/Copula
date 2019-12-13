nac_Node <- structure(list(), class = "nac_node")
nac_Node <- list()
class(nac_Node) <- "nac_node"

is.nac_Node <- function(nac_Node) {
  if (inherits(nac_Node, "nac_node")) { return(TRUE) }
  else {return(FALSE)}
}

new_nac_node <- function(list(family = character(), theta = double(), U_indices = integer(), subcopula = NULL) {
  for (i in 1:length(subcopula)) {
    stopifnot(is.nac_Node(subcopula[[i]]))
  }
  stopifnot(is.character(family))
  stopifnot(is.double(theta))
  stopifnot(is.integer(U_indices))
  stopifnot(is.list(subcopula))
  
  nac_Node <- list(family, theta, U_indices, subcopula)
  structure(nac_Node, class = "nac_node")
}

get_family <- function(nac_Node) {
  return(nac_Node[[1]])
}

get_iniTheta <- function(nac_Node) {
  return(nac_Node[[2]])
}

get_U_indices <- function(nac_Node) {
  return(nac_Node[[3]])
}

get_subcopulas <- function(nac_Node) {
  return(nac_Node[[4]])
}

set_family <- function(nac_Node, family) {
  nac_Node[[1]] <- family
}

set_iniTheta <- function(nac_Node, theta) {
  nac_Node[[2]] <- theta
}

set_U_indices <- function(nac_Node, index) {
  nac_Node[[3]] <- index
}

set_subcopula <- function(nac_Node, subcopula) {
  nac_Node[[4]] <- subcopula
}

has_subcopulas <- function(nac_Node) {
  if (length(nac_Node) == 4) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

append_subcopula <- function(nac_Node, subcopula) {
  nac_Node <- append(nac_Node, subcopula)
  return(nac_Node)
}

set_theta_from_unbounded <- function(nac_Node, par, current_index = 1) {
  family <- get_family(nac_Node)
  
  if (family == "Clayton") {
    if (has_subcopulas(nac_Node) == FALSE) {
      if (length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- -1 + exp(par[current_index])
      } else {
        par[current_index] <- exp(par[current_index])
      }
      current_index <- current_index + 1
    }
    
    else if (has_subcopulas(nac_Node)) {
      if (length(get_subcopulas(nac_Node)) + length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- -1 + exp(par[current_index])
      } else {
        par[current_index] <- exp(par[current_index])
      }
      current_index <- current_index + 1
      
      subcopulas_list <- get_subcopulas(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  if (family == "Frank") {
    if (has_subcopulas(nac_Node) == FALSE) {
      if (par[current_index] == 0) {
        par[current_index] = 0.000001
      }
      current_index <- current_index + 1
    }
    
    else if (has_subcopulas(nac_Node)) {
      if (par[current_index] == 0) {
        par[current_index] = 0.000001
      }
      
      current_index <- current_index + 1
      subcopulas_list <- get_subcopulas(nac_Node)
      
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  if (family == "Gumbel" | family == "Joe"){
    if (has_subcopulas(nac_Node) == FALSE) {
      par[current_index] <- exp(par[current_index]) + 1
      current_index <- current_index + 1
    }
    else if (has_subcopulas(nac_Node)) {
      par[current_index] <- exp(par[current_index]) + 1
      current_index <- current_index + 1
      subcopulas_list <- get_subcopulas(nac_Node)
      
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  if (family == "Ali") {
    if (has_subcopulas(nac_Node) == FALSE) {
      par[current_index] <- 2 * (exp(par[current_index])/1+exp(par[current_index])) -1
      current_index <- current_index + 1
    }
    else if (has_subcopulas(nac_Node)) {
      par[current_index] <- 2 * (exp(par[current_index])/1+exp(par[current_index])) -1
      current_index <- current_index + 1
      subcopulas_list <- get_subcopulas(nac_Node)
      
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_from_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  return(list(par, current_index))
}


set_theta_unbounded <- function(nac_Node, par, current_index = 1) {
  family <- get_family(nac_Node)
  
  if (family == "Clayton") {
    if (has_subcopulas(nac_Node) == FALSE) {
      if (length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- log(par[current_index] + 1)
      } else {
        par[current_index] <- log(par[current_index])
      }
      current_index <- current_index + 1
    }
    
    else if (has_subcopulas(nac_Node)) {
      if (length(get_subcopulas(nac_Node)) + length(get_U_indices(nac_Node)) == 2) {
        par[current_index] <- log(par[current_index] + 1)
      } else {
        par[current_index] <- log(par[current_index])
      }
      current_index <- current_index + 1
      
      subcopulas_list <- get_subcopulas(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  if (family == "Frank") {
    # do nothing for frank family
  }
  
  if (family == "Gumbel" | family == "Joe"){
    if (has_subcopulas(nac_Node) == FALSE) {
      par[current_index] <- log(par[current_index] - 1)
      current_index <- current_index + 1
    }
    else if (has_subcopulas(nac_Node)) {
      par[current_index] <- log(par[current_index] - 1)
      current_index <- current_index + 1
      
      subcopulas_list <- get_subcopulas(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  
  if (family == "Ali") {
    if (has_subcopulas(nac_Node) == FALSE) {
      par[current_index] <- 1 - 1/ (1 - (par[current_index] + 1)/2)
      current_index <- current_index + 1
    }
    else if (has_subcopulas(nac_Node)) {
      par[current_index] <- 1 - 1/ (1 - (par[current_index] + 1)/2)
      current_index <- current_index + 1
      
      subcopulas_list <- get_subcopulas(nac_Node)
      for (i in 1:length(subcopulas_list)) {
        list <- subcopulas_list[[i]]
        temp <- set_theta_unbounded(nac_Node, par, current_index)
        par[current_index] <- temp[[1]][[current_index]]
        current_index <- temp[[2]]
      }
    }
  }
  return(list(par, current_index))
}

get_cdf <- function(nac_Node, U) {
  nestedCopula <- 0
  family <- get_family(nac_Node)
  theta <- get_iniTheta(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  U_indices <- get_U_indices(nac_Node)
  
  if (has_subcopulas(nac_Node)) {
    subcopulas_list <- get_subcopulas(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))
    
    for (v_ind in 1:length(subcopulas_list)) {
      V[,v_ind] <- get_cdf(subcopulas_list[[v_ind]], U)
    }

    # check family
    if (family == 'Clayton') {
      copula <- claytonCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Frank') {
      copula <- frankCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Gumbel') {
      copula <- gumbelCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Independence') {
      copula <- indepCopula(dim = ncol_U + ncol(V))
    } else if (family == 'Joe') {
      copula <- joeCopula(theta, dim = ncol_U + ncol(V))
    } else if (family == 'Ali') {
      copula <- amhCopula(theta, dim = ncol_U + ncol(V))
    }
    
    nestedCopula <- nestedCopula + pCopula(cbind(U[,U_indices],V), copula = copula)
  }
  
  else if (has_subcopulas(nac_Node) == FALSE) {
   if (family == 'Clayton') {
      copula <- claytonCopula(theta, dim = ncol_U)
    } else if (family == 'Frank') {
      copula <- frankCopula(theta, dim = ncol_U)
    } else if (family == 'Gumbel') {
      copula <- gumbelCopula(theta, dim = ncol_U)
    } else if (family == 'Independence') {
      copula <- indepCopula(dim = ncol_U)
    } else if (family == 'Joe') {
      copula <- joeCopula(theta, dim = ncol_U)
    } else if (family == 'Ali') {
      copula <- amhCopula(theta, dim = ncol_U)
    }
    
    U_indices <- get_U_indices(nac_Node) 
    nestedCopula <- nestedCopula + pCopula(U[,U_indices], copula = copula)
  }
  return(nestedCopula)
}

get_density <- function(nac_Node, U, log = TRUE) {
  density <- 0
  family <- get_family(nac_Node)
  theta <- get_iniTheta(nac_Node)
  ncol_U <- length(get_U_indices(nac_Node))
  U_indices <- get_U_indices(nac_Node)
  
  if (has_subcopulas(nac_Node) == FALSE) {
    if (family == 'Clayton') {
      density <- density + dCopula(U[,U_indices], claytonCopula(theta, dim = ncol_U, log = TRUE))
    } else if (family == 'Frank') {
      density <- density + dCopula(U[,U_indices], frankCopula(theta, dim = ncol_U, log = TRUE))
    } else if (family == 'Gumbel') {
      density <- density + dCopula(U[,U_indices], gumbelCopula(theta, dim = ncol_U, log = TRUE))
    } else if (family == 'Independence') {
      density <- density + dCopula(U[,U_indices], indepCopula(dim = ncol_U, log = TRUE))
    } else if (family == 'Joe') {
      density <- density + dCopula(U[,U_indices], joeCopula(theta, dim = ncol_U, log = TRUE))
    } else if (family == 'Ali') {
      density <- density + dCopula(U[,U_indices], amhCopula(theta, dim = ncol_U, log = TRUE))
    }
  }
  
  else if (has_subcopulas(nac_Node)) {
    subcopulas_list <- get_subcopulas(nac_Node)
    V <- matrix(NA, nrow = nrow(U), ncol = length(subcopulas_list))
    
    for (v_ind in 1:length(subcopulas_list)) {
      density <- density + get_density(subcopulas_list[[v_ind]], U, log = TRUE)
      V[,v_ind] <- get_cdf(subcopulas_list[[v_ind]], U)
    }
    
    X <- cbind(U[,U_indices], V)
    # check family
    if (family == 'Clayton') {
      density <- density + dCopula(X, claytonCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Frank') {
      density <- density + dCopula(X, frankCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Gumbel') {
      density <- density + dCopula(X, gumbelCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Independence') {
      density <- density + dCopula(X, indepCopula(dim = ncol(X)), log = TRUE)
    } else if (family == 'Joe') {
      density <- density + dCopula(X, joeCopula(theta, dim = ncol(X)), log = TRUE)
    } else if (family == 'Ali') {
      density <- density + dCopula(X, amhCopula(theta, dim = ncol(X)), log = TRUE)
    }
  }
  
  if (log) {
    return(density)
  } else {
    return(exp(density))
  }
}

get_iniCor <- function(nac_Node, U) {
  theta <- 0
  family <- get_family(nac_Node)
  U_indices <- get_U_indices(nac_Node)
  empirical_taus <- cor(U[,U_indices], method = "kendall")
  
  n <- length(U_indices)
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
      if (has_subcopulas(nac_Node) == FALSE && length(get_U_indices(nac_Node)) == 2) {
        theta <- -1
      } else if (has_subcopulas(nac_Node) && length(get_U_indices(nac_Node)) + length(get_subcopulas(nac_Node)) == 2) {
        theta <- -1
      }
    }
  }
  else if (theta < 0) {
    if (family == "Clayton") {
      if (has_subcopulas(nac_Node) == FALSE && length(get_U_indices(nac_Node)) >= 3) {
        theta <- 0
      } else if (has_subcopulas(nac_Node) && length(get_U_indices(nac_Node)) + length(get_subcopulas(nac_Node)) >= 3) {
        theta <- 0
      }
    }
  }
  else if (theta <= 1) {
    if (family == "Joe" | family == "Gumbel") {
      theta <- 1.000001
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

get_num_par <- function(nac_Node) {
  num_par <- 0
  if (has_subcopulas(nac_Node)) {
    num_par <- num_par + 1
    sub_list <- get_subcopulas(nac_Node)
    
    for (i in 1:length(sub_list)) {
      num_par <- num_par + num_par(sub_list[[i]])
    }
  }
  
  else if (has_subcopulas(nac_Node) == FALSE) {
    num_par <- num_par + 1
  }
  
  return(num_par)
}

update_par <- function(nac_Node, par, current_index = 1) {
  if (has_subcopulas(nac_Node)) {
    set_iniTheta(nac_Node, par[current_index])
    # if(is.na(naclist[[2]])) {
    #     cat("in updatePars")
    #     browser()
    # }
    current_index <- current_index + 1
    get_subcopulas()
    subcopulas <- get_subcopulas(nac_Node)
    for (child_index in 1:subcopulas) {
      temp <- update_par(subcopulas[[child_index]], par, current_index)
      set_subcopula(subcopulas[[child_index]], temp[[1]])
      current_index <- temp[[2]]
    }
  }
  
  else if (has_subcopulas(nac_Node) == FALSE) {
    set_iniTheta(nac_Node, par[current_index])
    # if(is.na(naclist[[2]])) {
    #     cat("in updatePars")
    #     browser()
    # }
    current_index <- current_index + 1
  }
  
  return(list(naclist, current_index))
}


get_loglik <- function(nac_Node, U, par) {
  # transform provided par (vector of real numbers) to a vector of parameter values
  # that are within the bounds for the copula we are estimating
  orig_par <- par
  par <- set_theta_from_unbounded(nac_Node, par, current_index = 1)[[1]]  
  
  # if(is.na(par)) {
  #   browser()
  # }
  
  # evaluate log-likelihood at provided parameter values
  new_nac_Node <- update_par(nac_Node, par, current_index = 1)[[1]]
  result <- sum(get_density(new_nac_Node, U, log = TRUE))
  return(result)
}

get_mix_density <- function(U, W, mixlist, log = TRUE) {
  log_sum <- vector()
  density <- 0
  num_of_copulas <- length(mixlist)
  
  for (i in 1:num_of_copulas) {
    nac_Node <- mixlist[[i]]
    log_sum <- append(log_sum, log(W[[i]]) + get_density(nac_Node, U, log = TRUE))
  }
  
  density <- logspace_sum(log_sum)
  
  if(log) {
    return(density)
  } else {
    return(exp(density))
  }
}

get_mix_loglik <- function(U, mixlist, parlist) {
  num_of_copulas <- length(mixlist)
  sum_params <- 0
  W <- vector()
  
  # softmax transformation
  # transform par into weights
  for (i in 1:num_of_copulas) {
    sum_params <- sum_params + exp(parlist[[i]])
  }
  
  for (i in 1:num_of_copulas) {
    W[[i]] <- exp(parlist[[i]]) / sum_params
  }
  
  # evaluate log-likelihood
  tryCatch(sum(get_mix_density(U, W, mixlist, log = TRUE)),
           error = function(e) {
             return(-1 * .Machine$double.xmax)
           })
}

get_xv <- function(x, nac_Node, k = NULL) {
  if (!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  
  d <- ncol(x)
  n <- nrow(x)
  stopifnot(is.numeric(x), d > 1, n > 0)
  
  if (is.null(k)) { k <- n}
  else {as.integer(k)}
  
  density <- 0
  folds <- createFolds(1:n, k, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:k) {
    test_ind <- c(folds[[i]])
    train_set <- pobs(x[-test_ind, , drop = FALSE], ties.method = 'average')
    optim_results <- optim(par = param,
                           fn = get_loglik,
                           U = train_set,
                           naclist = naclist,
                           method ="BFGS",
                           control = list(fnscale = -1))
    
    full_set <- pobs(x, ties.method = 'average')
    test_set <- full_set[test_ind, ]
    
    new_nac_Node <- update_par(nac_Node, exp(optim_results$par), current_index = 1)[[1]]
    density <- density + sum(get_density(new_nac_Node, test_set, log = TRUE))
  }
  density/k
}