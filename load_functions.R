
get_time <- function(x, data){
  t = data$time
  H2 = data$h2
  
  if (x < min(H2)) {
    t1 = t[1]
  } else {
    t1 = t[max(which(!x < H2))]
  }
  
  return(t1)
}

get_time_H0 <- function(x, data){
  t = data$time
  H = data$hazard
  
  if (x < min(H) | is.na(x)) {
    #t1 = t[1]
    t1 = min( t[t!=min(t)] )
  } else {
    t1 = t[max(which(!x < H))]
  }
  
  return(t1)
}

#install.packages("flexsurv")
library("flexsurv")
