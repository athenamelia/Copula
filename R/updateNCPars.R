# update nested copula parameters
updatePars <- function(par, naclist, current_index = 1) {
  if (length(naclist) == 3) {
    naclist[[1]] <- par[current_index]
    current_index <- current_index + 1

    for (child_index in 1:length(naclist[[3]])) {
      temp <- updatePars(par, naclist[[3]][[child_index]], current_index)
      naclist[[3]][[child_index]] <- temp[[1]]
      current_index <- temp[[2]]
    }
  }
  
  else if (length(naclist) == 2) {
    naclist[[1]] <- par[current_index]
    current_index <- current_index + 1
  }
  
  return(list(naclist, current_index))
}

param <- c(1, 2, 3)
nlist <- list(th[1], 1:2, list(list(th[2], 3:4),
                              list(th[3], 5:6)))
updatePars(param, nlist, current_index = 1)



# NAC structure
nlist_test <- list(th[1], 1, list(list(th[2], 2:3, list(list(th[3], 4, list(list(th[4], 5:6))))), 
                                  list(th[5], 7:10),
                                  list(th[6], 11, list(list(th[7], 12:15, list(list(th[8], 16)))))))
  
nlist_test
param_test <- c(1, 2, 3, 4, 5, 6, 7, 8)
updatePars(param_test, nlist_test, current_index = 1)
