
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

