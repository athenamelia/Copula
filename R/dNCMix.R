# Calculate density of mixture of copulas 

dNCMix <- function(U, W, mixlist, log = TRUE) {
  log_sum <- vector()
  density <- 0 
  num_of_copulas <- length(mixlist) 
  
  for (i in 1:num_of_copulas) {
    log_sum <- append(log_sum, log(W[[i]]) + dNC(U, mixlist[[i]], log = TRUE))
  }
  
  density <- logspace_sum(log_sum)
  
  if(log) {
    return(density)  
  } else {
    return(exp(density))
  }
}

