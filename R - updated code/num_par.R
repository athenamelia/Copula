# calculate number of parameters

# @param naclist nesting structure
# @return number of parameters

num_par <- function(naclist) {
  num_par <- 0
  if (length(naclist) == 4) {
    num_par <- num_par + 1
    sub_list <- naclist[[4]]

    for (i in 1:length(sub_list)) {
      num_par <- num_par + num_par(sub_list[[i]])
    }
  }

  else if (length(naclist) == 3) {
      num_par <- num_par + 1
  }

  return(num_par)
}
